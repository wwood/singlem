#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2023 Ben Woodcroft
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2023"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import logging
import sys
import os
import csv
import shutil
from multiprocessing import Pool
import tempfile
import polars as pl

import extern
from bird_tool_utils import iterable_chunks

from .biolib_lite.prodigal_biolib import Prodigal
from .biolib_lite.common import remove_extension
from tqdm.contrib.concurrent import process_map
from tqdm import tqdm

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '..')] + sys.path

from .metapackage import Metapackage
from .archive_otu_table import ArchiveOtuTable
from .sequence_classes import SeqReader
from .singlem import OrfMUtils
from .otu_table_collection import OtuTableCollection
from .regenerator import Regenerator
from .pipe import SearchPipe
from .sequence_database import SequenceDatabaseOtuTable, SequenceDatabase
from .checkm2 import CheckM2

taxonomy_prefixes = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
genome_file_suffixes = ['.fna', '.fa', '.fasta', '.fna.gz', '.fa.gz', '.fasta.gz']


def generate_taxonomy_for_new_genomes(**kwargs):
    threads = kwargs.pop('threads')
    new_genome_fasta_files = kwargs.pop('new_genome_fasta_files')
    pplacer_threads = kwargs.pop('pplacer_threads')
    working_directory = kwargs.pop('working_directory')
    gtdbtk_output_directory = kwargs.pop('gtdbtk_output_directory')
    output_taxonomies_file = kwargs.pop('output_taxonomies_file')
    excluded_genomes = kwargs.pop('excluded_genomes')
    if len(kwargs) > 0:
        raise Exception("Unexpected arguments detected: %s" % kwargs)

    # Create batchfile, tab separated in 2 columns (FASTA file, genome ID)
    name_to_genome_fasta = {}
    with tempfile.NamedTemporaryFile(mode='w') as batchfile:
        for genome_fasta in new_genome_fasta_files:
            if not os.path.isfile(genome_fasta):
                raise Exception("Specified genome does not appear to be a file: {}".format(genome_fasta))
            name = os.path.basename(genome_fasta)
            name_to_genome_fasta[name] = genome_fasta
            batchfile.write(genome_fasta + '\t' + name + '\n')
        batchfile.flush()

        if gtdbtk_output_directory:
            logging.info("Using pre-existing GTDBtk output directory: {}".format(gtdbtk_output_directory))
            gtdbtk_output = gtdbtk_output_directory
        else:
            # run gtdbtk to temporary output directory
            gtdbtk_output = os.path.join(working_directory, 'gtdbtk_output')
            logging.info("Running GTDBtk to generate taxonomy for new genomes ..")
            # logging.warning("mash_db used is specific to QUT's CMR cluster, will fix this in future")
            cmd = f'gtdbtk classify_wf --cpus {threads} --batchfile {batchfile.name} --out_dir {gtdbtk_output} --mash_db {working_directory}/gtdbtk_mash.msh'
            if pplacer_threads:
                cmd += f' --pplacer_cpus {pplacer_threads}'
            extern.run(cmd)
            logging.info("GTDBtk finished")

        # read in assigned taxonomy from gtdbtk output
        if not os.path.exists(gtdbtk_output):
            raise Exception("GTDBtk output directory not found: {}".format(gtdbtk_output))
        bacteria_taxonomy_path = os.path.join(gtdbtk_output, 'gtdbtk.bac120.summary.tsv')
        archaea_taxonomy_path = os.path.join(gtdbtk_output, 'gtdbtk.ar53.summary.tsv')

        if output_taxonomies_file:
            logging.info("Writing new genome taxonomies to {}".format(output_taxonomies_file))
            output_taxonomies_fh = open(output_taxonomies_file, 'w')
            output_taxonomies_fh.write('genome\ttaxonomy\n')

        new_taxonomies = {}
        excluded_genome_basenames = set([os.path.basename(x) for x in excluded_genomes])
        for tax_file in [bacteria_taxonomy_path, archaea_taxonomy_path]:
            if os.path.exists(tax_file):
                logging.info("Reading taxonomy from {}".format(tax_file))

                df = pl.read_csv(tax_file, separator='\t')
                for row in df.rows(named=True):
                    genome_name = row['user_genome']
                    taxonomy = row['classification'].split(';')
                    if genome_name in excluded_genome_basenames:
                        logging.debug("Ignoring genome {} because it is in the excluded_genomes list".format(genome_name))
                        continue
                    if len(taxonomy) != 7:
                        if taxonomy == ['Unclassified Bacteria'] or taxonomy == ['Unclassified Archaea']:
                            logging.warning(
                                "The genome {} was not given any classification, perhaps it is poor quality?".format(
                                    genome_name))
                            if output_taxonomies_file:
                                output_taxonomies_fh.write(
                                    '\t'.join([genome_name, 'Root; ' + '; '.join(row['classification'].split(';'))]) +
                                    '\n')
                        else:
                            raise Exception("Unexpected taxonomy length found in GTDBtk output file {}: {}".format(
                                tax_file, taxonomy))
                    elif taxonomy[6] != 's__':
                        logging.debug(
                            "Ignoring genome {} because it already has a species-level taxonomic assignment: {}".format(
                                genome_name, row['classification']))
                        if output_taxonomies_file:
                            output_taxonomies_fh.write(
                                '\t'.join([genome_name, 'Root; ' + '; '.join(row['classification'].split(';'))]) + '\n')
                    else:
                        placeholder_taxonomy = remove_file_extensions(genome_name)

                        for i, prefix in zip(range(len(taxonomy)), taxonomy_prefixes):
                            if taxonomy[i] == prefix:
                                taxonomy[i] = prefix + placeholder_taxonomy
                        logging.debug("Adding new species-level genome {} with taxonomy {}".format(
                            genome_name, taxonomy))
                        new_taxonomies[name_to_genome_fasta[genome_name]] = taxonomy
                        if output_taxonomies_file:
                            output_taxonomies_fh.write('\t'.join([genome_name, 'Root; ' + '; '.join(taxonomy)]) + '\n')

        if output_taxonomies_file:
            output_taxonomies_fh.close()

        if len(new_taxonomies) == 0:
            raise Exception(
                "No new species-level genomes found in GTDBtk output, so there are no genomes to add to the metapackage."
            )

        # Create taxonomy file
        logging.info("Creating taxonomy file for new genomes")
        taxonomy_path = os.path.join(working_directory, 'new_genomes_taxonomy.tsv')
        with open(taxonomy_path, 'w') as taxonomy_file:
            writer = csv.writer(taxonomy_file, delimiter='\t')
            for genome_fasta, taxonomy in new_taxonomies.items():
                writer.writerow([genome_fasta] + ['; '.join(taxonomy)])
        logging.debug("New genomes taxonomy file written to {}".format(taxonomy_path))

        # Remove fasta files where the species is not new
        initial_count = len(new_genome_fasta_files)
        new_genome_fasta_files = list([g for g in new_genome_fasta_files if g in new_taxonomies])
        logging.info(
            "Removed {} genomes from new genomes list because they already had a species-level taxonomic assignment".
            format(initial_count - len(new_genome_fasta_files)))
        logging.info("Number of new genomes to add to metapackage: {}".format(len(new_genome_fasta_files)))

        return taxonomy_path, new_genome_fasta_files

def remove_file_extensions(filename):
    for suffix in genome_file_suffixes:
        if filename.endswith(suffix):
            filename = filename[:-len(suffix)]
    return filename

def generate_new_singlem_package(myargs):
    (working_directory, old_spkg, matched_transcript_path, sequence_to_genome, genome_to_taxonomy) = myargs

    logging.info("Preparing GraftM create inputs for {}".format(old_spkg.graftm_package_basename()))

    # Change working dir so that each spkg is isolated
    working_directory = os.path.join(working_directory, 'create-' + old_spkg.graftm_package_basename())
    os.makedirs(working_directory)

    # Create new taxonomy
    # new_gpkg_path = os.path.join(working_directory, old_spkg.graftm_package_basename() + '.gpkg')
    new_spkg_path = os.path.join(working_directory, old_spkg.graftm_package_basename() + '.spkg')

    # Run graftM graft to get the protein sequences from the transcripts which matched the HMM
    graft_output_folder = os.path.join(working_directory, 'graft_output')
    extern.run('graftM graft --graftm_package {} --forward {} --output {} --search_only'.format(
        old_spkg.graftm_package_path(), matched_transcript_path, graft_output_folder))

    # Work out the sequence prefix
    with open(old_spkg.graftm_package().unaligned_sequence_database_path()) as f:
        first_line = f.readline()
        if first_line.startswith('>'):
            if '~' in first_line:
                sequence_prefix = first_line[1:].split('~')[0]
            else:
                logging.warning(
                    "Using an empty sequence prefix to delineate different spkgs in sequence names, because the old spkg  didn't specify one."
                )
                sequence_prefix = ''
        else:
            raise Exception("Unexpected first line in unaligned sequence database: {}".format(first_line))

    # Find protein sequences that matched from each genome
    matched_sequences = []
    orf_path = os.path.join(graft_output_folder, 'matched_transcripts/matched_transcripts_orf.fa')
    if not os.path.exists(orf_path):
        # No matches, so just copy the old spkg
        logging.info("No hits for {}, so reusing copying old spkg".format(old_spkg.graftm_package_basename()))
        shutil.copytree(old_spkg.base_directory(), new_spkg_path)

    else:
        with open(orf_path) as f:
            for name, seq, _ in SeqReader().readfq(f):
                matched_sequences.append((name, seq))
        logging.info("Matched {} sequences from new genomes".format(len(matched_sequences)))

        # Add the sequences from the old spkg, but only if they have taxonomy assigned
        new_sequences_fasta = os.path.join(working_directory, 'new_sequences.fasta')
        shutil.copyfile(old_spkg.graftm_package().unaligned_sequence_database_path(), new_sequences_fasta)

        # Append them to the new sequences fasta and taxonomy files
        taxonomy_file = os.path.join(working_directory, 'taxonomy.tsv')
        orfm_utils = OrfMUtils()
        with open(taxonomy_file, 'w') as f:
            old_taxhash = old_spkg.taxonomy_hash()
            for name, taxonomy in old_taxhash.items():
                f.write(name + '\t' + '; '.join(taxonomy) + '\n')

            with open(new_sequences_fasta, 'a') as new_sequences_fasta_fh:
                for seq_name, seq in matched_sequences:
                    genome = sequence_to_genome[orfm_utils.un_orfm_name(seq_name)]
                    new_sequences_fasta_fh.write('>' + sequence_prefix + '~' + seq_name + '\n' + seq + '\n')
                    f.write(sequence_prefix + '~' + seq_name + '\t' + genome_to_taxonomy[genome] + '\n')

        # Log how many of the genomes matched this marker
        logging.info("Added {} sequences from {} genomes".format(len(matched_sequences), len(genome_to_taxonomy)))

        # Run singlem regenerate to create a new spkg
        logging.info("Regenerating {} ..".format(old_spkg.graftm_package_basename()))

        Regenerator().regenerate(
            input_singlem_package=old_spkg.base_directory(),
            output_singlem_package=new_spkg_path,
            input_sequences=new_sequences_fasta,
            input_taxonomy=taxonomy_file,
            sequence_prefix=sequence_prefix,
            no_further_euks=True,
        )

        logging.info("Regenerated {} ..".format(old_spkg.graftm_package_basename()))

    return new_spkg_path


def run_pipe(params):
    transcript_paths, tf, old_metapackage_path, threads = params
    # Set it high as it is too verbose otherwise
    logging.basicConfig(level=logging.ERROR, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    SearchPipe().run(sequences=transcript_paths,
                     archive_otu_table=tf,
                     metapackage_path=old_metapackage_path,
                     threads=threads)


def generate_new_metapackage(num_threads, working_directory, old_metapackage_path, new_genome_fasta_files,
                             new_taxonomies_file, new_genome_transcripts_and_proteins):

    # Add the new genome data to each singlem package
    # For each package, the unaligned seqs are in the graftm package,
    # Taxonomy/seqinfo give the taxonomy of target and euk sequences
    old_metapackage = Metapackage.acquire(old_metapackage_path)

    genome_to_taxonomy = {}
    with open(new_taxonomies_file) as f:
        for row in csv.reader(f, delimiter='\t'):
            genome = row[0]
            if genome.endswith('.fasta'):
                genome = genome[:-6]
            elif genome.endswith('.fasta.gz'):
                genome = genome[:-9]
            elif genome.endswith('.fa'):
                genome = genome[:-3]
            elif genome.endswith('.fa.gz'):
                genome = genome[:-6]
            elif genome.endswith('.fna'):
                genome = genome[:-4]
            elif genome.endswith('.fna.gz'):
                genome = genome[:-7]
            else:
                raise Exception("Unexpected file extension for genome in new_taxonomies_file: {}".format(genome))
            genome_to_taxonomy[os.path.basename(genome)] = row[1]
    if len(genome_to_taxonomy) != len(new_genome_fasta_files) or len(genome_to_taxonomy) == 0:
        raise

    logging.info("Gathering OTUs from new genomes ..")
    # Multiprocess the singlem pipe commands to speed things up, since pipe is
    # slow for many genomes at the moment, especially since we need to assign
    # taxonomy to exclude off-target hits.
    old_metapackage_on_new_genomes_tempfiles = []
    # Create regular files in a tempdir so that we don't need so many file handles
    td = tempfile.TemporaryDirectory(prefix='supplement_metapackage_otu_table_genome_jsons')
    param_sets = []
    for i, chunk in enumerate(iterable_chunks(new_genome_fasta_files, 10)):
        chunk2 = [c for c in chunk if c is not None]
        tf = os.path.join(td.name, f"{i}.json")
        old_metapackage_on_new_genomes_tempfiles.append(tf)
        transcript_paths = [new_genome_transcripts_and_proteins[c].transcript_fasta for c in chunk2]
        param_sets.append((transcript_paths, tf, old_metapackage_path, 1))
    process_map(run_pipe, param_sets, max_workers=num_threads, desc="Running singlem pipe")

    # Create SDB for final metapackage
    logging.info("Creating sdb to include in metapackage ..")
    new_metapackage_sdb_path = os.path.join(working_directory, 'new_metapackage.sdb')
    all_new_otu_tables = OtuTableCollection()
    # Change the taxonomy to be correct
    # old_metapackage_on_new_genomes_no_off_target_collection = StreamingOtuTableCollection(
    # )  # Need streaming so each_sample_otus works
    for tf in old_metapackage_on_new_genomes_tempfiles:
        with open(tf) as g:
            ar = ArchiveOtuTable.read(g)
            # old_metapackage_on_new_genomes_no_off_target_collection.add_archive_otu_table_object(ar)

            collection = OtuTableCollection()
            collection.add_archive_otu_table_object(ar)
            all_new_otu_tables.add_archive_otu_table_object(
                collection.exclude_off_target_hits(old_metapackage.singlem_packages, return_archive_table=True))
            # Set the source data's taxonomy
            taxonomy_field = ArchiveOtuTable.TAXONOMY_FIELD_INDEX
            sample_name_field = ArchiveOtuTable.SAMPLE_ID_FIELD_INDEX
            for row in collection.archive_table_objects[0].data:
                row[taxonomy_field] = genome_to_taxonomy[row[sample_name_field]]

    # Add old metapackage's OTUs.
    db_otu_table = SequenceDatabaseOtuTable(old_metapackage.nucleotide_sdb())
    all_new_and_old_otu_tables = OtuTableCollection()
    all_new_and_old_otu_tables.add_otu_table_object(all_new_otu_tables)
    all_new_and_old_otu_tables.add_otu_table_object(db_otu_table)

    SequenceDatabase().create_from_otu_table(
        new_metapackage_sdb_path,
        all_new_and_old_otu_tables,
        num_threads=num_threads,
    )

    # Dump matched transcript sequences to a file for graftm graft
    matched_transcript_path = os.path.join(working_directory, 'matched_transcripts.fasta')
    # with open(matched_transcript_path, 'w') as matched_transcripts_fasta:
    #     Summariser().dump_raw_sequences_from_archive_otu_table(
    #         output_table_io=matched_transcripts_fasta,
    #         table_collection=old_metapackage_on_new_genomes_no_off_target_collection,
    #     )

    sequence_to_genome = {}
    with open(matched_transcript_path, 'w') as f:
        for otu in all_new_otu_tables:
            read_names = otu.read_names()
            seqs = otu.read_unaligned_sequences()
            for seq, read_name in zip(seqs, read_names):
                if read_name in sequence_to_genome:
                    raise Exception("Duplicate sequence name: {}".format(read_name))
                sequence_to_genome[read_name] = otu.sample_name
                f.write('>{}\n{}\n'.format(read_name, seq))

    new_metapackage_path = os.path.join(working_directory, 'new_metapackage.mpkg')

    # For each spkg in the old mpkg, create a new spkg
    to_process = [(working_directory, spkg, matched_transcript_path, sequence_to_genome, genome_to_taxonomy)
                  for spkg in old_metapackage.singlem_packages]
    if num_threads > 1:
        new_spkg_paths = Pool(num_threads).map(generate_new_singlem_package, to_process)
    else:
        new_spkg_paths = list([generate_new_singlem_package(x) for x in to_process])

    # Create a new metapackage from the singlem packages
    logging.info("Creating new metapackage ..")
    Metapackage.generate(singlem_packages=new_spkg_paths,
                         nucleotide_sdb=new_metapackage_sdb_path,
                         prefilter_diamond_db=old_metapackage.prefilter_db_path(),
                         output_path=new_metapackage_path,
                         threads=num_threads,
                         prefilter_clustering_threshold=None,
                         taxon_genome_lengths=None)
    logging.info("New metapackage created at {}".format(new_metapackage_path))

    td.cleanup()

    return new_metapackage_path


class TranscriptsAndProteins:

    def __init__(self):
        self.transcript_fasta = None
        self.protein_fasta = None
        self.translation_table = None


def run_prodigal_on_one_genome(params):
    genome_fasta, output_directory = params
    logging.debug("Running prodigal on {}, output {}".format(genome_fasta, output_directory))
    prodigal = Prodigal(cpus=1, verbose=False)
    prodigal_return = prodigal.run([genome_fasta], output_directory)

    genome_id = remove_extension(genome_fasta)
    # These must match result from biolib, which I modified
    aa_gene_file = os.path.join(output_directory, genome_id + '.faa')
    nt_gene_file = os.path.join(output_directory, genome_id + '.fna')

    to_return = TranscriptsAndProteins()
    to_return.protein_fasta = aa_gene_file
    to_return.transcript_fasta = nt_gene_file
    to_return.translation_table = prodigal_return[genome_id].best_translation_table
    return to_return


def generate_faa_and_transcript_fna_files_for_new_genomes(**kwargs):
    # working_directory=working_directory,
    # threads=args.threads,
    # new_genome_fasta_files=args.new_genome_fasta_files
    working_directory = kwargs.pop('working_directory')
    threads = kwargs.pop('threads')
    new_genome_fasta_files = kwargs.pop('new_genome_fasta_files')
    if len(kwargs) > 0:
        raise Exception("Unexpected arguments detected: %s" % kwargs)

    # Run prodigal on the new genomes, returning a translation table (4 or 11)
    # appropriate one.
    output_directory_base = os.path.join(working_directory, 'prodigal_output')
    if not os.path.exists(output_directory_base):
        os.makedirs(output_directory_base)
    genome_ids = set()
    param_sets = []
    for genome_fasta in new_genome_fasta_files:
        genome_id = remove_extension(genome_fasta)

        if genome_id in genome_ids:
            raise Exception("Duplicate genome ID: %s" % genome_id)
        genome_ids.add(genome_id)

        output_directory = os.path.join(output_directory_base, genome_id)
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)
        param_sets.append((genome_fasta, output_directory))

    logging.info("Running prodigal on {} new genomes ..".format(len(param_sets)))
    if threads > 1:
        transcripts_and_proteins = process_map(run_prodigal_on_one_genome, param_sets, max_workers=threads)
    else:
        transcripts_and_proteins = []
        for x in tqdm(param_sets):
            transcripts_and_proteins.append(run_prodigal_on_one_genome(x))

    num_translation_table4 = sum([1 for x in transcripts_and_proteins if x.translation_table == 4])
    num_translation_table11 = sum([1 for x in transcripts_and_proteins if x.translation_table == 11])
    logging.info("Prodigal finished. Found {} with translation table 4, {} with translation table 11".format(
        num_translation_table4, num_translation_table11))

    to_return = {}
    for genome_fasta, transcripts_and_proteins in zip(new_genome_fasta_files, transcripts_and_proteins):
        to_return[genome_fasta] = transcripts_and_proteins
    return to_return


class Supplementor:

    def supplement(self, **kwargs):
        new_genome_fasta_files = kwargs.pop('new_genome_fasta_files')
        new_taxonomies = kwargs.pop('new_taxonomies')
        input_metapackage = kwargs.pop('input_metapackage')
        output_metapackage = kwargs.pop('output_metapackage')
        threads = kwargs.pop('threads')
        pplacer_threads = kwargs.pop('pplacer_threads')
        working_directory = kwargs.pop('working_directory')
        gtdbtk_output_directory = kwargs.pop('gtdbtk_output_directory')
        output_taxonomies = kwargs.pop('output_taxonomies')
        # checkm2_quality_file=args.checkm2_quality_file,
        # no_quality_filter=args.no_quality_filter,
        checkm2_quality_file = kwargs.pop('checkm2_quality_file')
        checkm2_min_completeness = kwargs.pop('checkm2_min_completeness')
        checkm2_max_contamination = kwargs.pop('checkm2_max_contamination')
        no_quality_filter = kwargs.pop('no_quality_filter')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        with tempfile.TemporaryDirectory() as working_directory:
            if working_directory:
                working_directory = working_directory
            if not os.path.exists(working_directory):
                os.mkdir(working_directory)

            # Quality filter the new genomes
            num_before_quality_filter = len(new_genome_fasta_files)
            excluded_genomes = set()
            if no_quality_filter:
                logging.info("Skipping quality filtering of new genomes")
            elif checkm2_quality_file:
                logging.info("Quality filtering new genomes using checkm2 quality file: %s" % checkm2_quality_file)
                checkm = CheckM2(checkm2_quality_file)
                new_genome_fasta_files_without_file_extension = set(checkm.genomes_of_sufficient_quality(checkm2_min_completeness, checkm2_max_contamination))
                # Check the new genomes are in the checkm2 quality file
                for g in new_genome_fasta_files:
                    if remove_extension(g) not in checkm:
                        raise Exception("Genome %s not found in checkm2 quality file" % g)
                new_genome_fasta_files2 = []
                for x in new_genome_fasta_files:
                    if remove_extension(x) in new_genome_fasta_files_without_file_extension:
                        new_genome_fasta_files2.append(x)
                    else:
                        excluded_genomes.add(x)
                new_genome_fasta_files = new_genome_fasta_files2
                logging.info("Removed {} genomes of insufficient quality, leaving {} genomes to be added to the metapackage".format(
                    num_before_quality_filter - len(new_genome_fasta_files), len(new_genome_fasta_files)))
            else:
                raise Exception("Must provide either --no-quality-filter or --checkm2-quality-file")
            

            # Generate faa and transcript fna files for the new genomes
            logging.info("Generating faa and transcript fna files for the new genomes ..")
            genome_transcripts_and_proteins = generate_faa_and_transcript_fna_files_for_new_genomes(
                working_directory=working_directory, threads=threads, new_genome_fasta_files=new_genome_fasta_files)

            if new_taxonomies:
                taxonomy_file = new_taxonomies
                new_genome_fasta_files = new_genome_fasta_files
            else:
                logging.info("Generating taxonomy for new genomes ..")
                taxonomy_file, new_genome_fasta_files = generate_taxonomy_for_new_genomes(
                    working_directory=working_directory,
                    threads=threads,
                    new_genome_fasta_files=new_genome_fasta_files,
                    gtdbtk_output_directory=gtdbtk_output_directory,
                    pplacer_threads=pplacer_threads,
                    output_taxonomies_file=output_taxonomies,
                    excluded_genomes=excluded_genomes)

            # Run the genomes through pipe with genome fasta input to identify the new sequences
            logging.info("Generating new SingleM packages and metapackage ..")
            new_metapackage = generate_new_metapackage(threads, working_directory, input_metapackage,
                                                       new_genome_fasta_files, taxonomy_file,
                                                       genome_transcripts_and_proteins)

            logging.info("Copying generated metapackage to {}".format(output_metapackage))
            shutil.copytree(new_metapackage, output_metapackage)
            logging.info("Finished suppplement.")
