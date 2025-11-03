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
from multiprocessing import Pool, Manager, get_context
import tempfile
import re
os.environ["POLARS_MAX_THREADS"] = str(1)
os.environ["OPENBLAS_NUM_THREADS"] = str(1)
import polars as pl
from Bio import SearchIO

import extern
import subprocess

from .biolib_lite.prodigal_biolib import Prodigal
from .biolib_lite.common import remove_extension
from tqdm.contrib.concurrent import process_map
from tqdm import tqdm

# sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__))), '..'] + sys.path

from .metapackage import Metapackage, CUSTOM_TAXONOMY_DATABASE_NAME, GTDB_DATABASE_NAME
from .archive_otu_table import ArchiveOtuTable
from .sequence_classes import SeqReader
from .utils import FastaNameToSampleName, OrfMUtils
from .otu_table_collection import OtuTableCollection
from .regenerator import Regenerator
from .pipe import SearchPipe
from .sequence_database import SequenceDatabaseOtuTable, SequenceDatabase
from .checkm2 import CheckM2
from .otu_table_entry import OtuTableEntry
from .genome_size import GenomeSizes

taxonomy_prefixes = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
genome_file_suffixes = ['.fna', '.fa', '.fasta', '.fna.gz', '.fa.gz', '.fasta.gz']


def generate_taxonomy_for_new_genomes(**kwargs):
    threads = kwargs.pop('threads')
    new_genome_fasta_files = kwargs.pop('new_genome_fasta_files')
    pplacer_threads = kwargs.pop('pplacer_threads')
    working_directory = kwargs.pop('working_directory')
    gtdbtk_output_directory = kwargs.pop('gtdbtk_output_directory')
    taxonomy_file = kwargs.pop('taxonomy_file')
    output_taxonomies_file = kwargs.pop('output_taxonomies_file')
    excluded_genomes = kwargs.pop('excluded_genomes')
    old_metapackage = kwargs.pop('old_metapackage')
    skip_taxonomy_check = kwargs.pop('skip_taxonomy_check')
    if len(kwargs) > 0:
        raise Exception("Unexpected arguments detected: %s" % kwargs)

    if not skip_taxonomy_check:
        logging.info("Reading old taxonomies for checking concordance with new taxonomy strings ..")
        known_taxons = set()
        for spkg in old_metapackage.singlem_packages:
            for taxon in spkg.taxonomy_hash().values():
                taxon_building = None
                for t in taxon:
                    if taxon_building is None:
                        taxon_building = t
                    else:
                        taxon_building += ';' + t
                    known_taxons.add(taxon_building)
        logging.info("Read in {} known taxons".format(len(known_taxons)))

    # Create batchfile, tab separated in 2 columns (FASTA file, genome ID)
    name_to_genome_fasta = {}
    for genome_fasta in new_genome_fasta_files:
        if not os.path.isfile(genome_fasta):
            raise Exception("Specified genome does not appear to be a file: {}".format(genome_fasta))
        name = os.path.basename(genome_fasta)
        name_to_genome_fasta[name] = genome_fasta
        name_to_genome_fasta[remove_file_extensions(name)] = genome_fasta

    taxonomies_to_process = []
    if taxonomy_file:
        logging.info("Using pre-existing taxonomy file: {}".format(taxonomy_file))
        # Read in taxonomy file as 2col TSV
        with open(taxonomy_file) as f:
            for row in csv.reader(f, delimiter='\t'):
                if len(row) != 2:
                    raise Exception("Unexpected number of columns (expected 2) in taxonomy file: {}".format(row))
                taxonomies_to_process.append((row[0], row[1]))

    else:
        with tempfile.NamedTemporaryFile(mode='w') as batchfile:
            for genome_fasta in new_genome_fasta_files:
                name = os.path.basename(genome_fasta)
                batchfile.write(genome_fasta + '\t' + name + '\n')
            batchfile.flush()

            if gtdbtk_output_directory:
                logging.info("Using pre-existing GTDBtk output directory: {}".format(gtdbtk_output_directory))
                gtdbtk_output = gtdbtk_output_directory
            else:
                # run gtdbtk to temporary output directory
                gtdbtk_output = os.path.join(working_directory, 'gtdbtk_output')
                logging.info("Running GTDBtk to generate taxonomy for new genomes ..")
                if not shutil.which("gtdbtk"):
                    raise Exception("gtdbtk is not installed by default; install it (e.g. 'conda install gtdbtk').")
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

            for tax_file in [bacteria_taxonomy_path, archaea_taxonomy_path]:
                if os.path.exists(tax_file):
                    logging.info("Reading taxonomy from {}".format(tax_file))

                    df = pl.read_csv(tax_file, separator='\t', null_values = "N/A")
                    for row in df.rows(named=True):
                        genome_name = row['user_genome']
                        taxonomies_to_process.append((genome_name, row['classification']))

    # For loop setup
    excluded_genome_basenames = set([os.path.basename(x) for x in excluded_genomes])
    new_taxonomies = {}
    num_genomes_without_fasta = 0
    if output_taxonomies_file:
        logging.info("Writing new genome taxonomies to {}".format(output_taxonomies_file))
        output_taxonomies_fh = open(output_taxonomies_file, 'w')
        output_taxonomies_fh.write('genome\ttaxonomy\n')

    # For loop
    for genome_name, taxonomy_str in taxonomies_to_process:
        taxonomy = list([s.strip() for s in taxonomy_str.split(';')])

        if genome_name in excluded_genome_basenames:
            logging.debug(
                "Ignoring genome {} because it is in the excluded_genomes list".format(genome_name))
            continue
        if genome_name not in name_to_genome_fasta:
            logging.debug(
                "Genome {} was not found in the list of genomes to be included".format(genome_name))
            num_genomes_without_fasta += 1
            continue

        if len(taxonomy) != 7:
            if taxonomy == ['Unclassified Bacteria'] or taxonomy == ['Unclassified Archaea']:
                logging.warning(
                    "The genome {} was not given any classification, perhaps it is poor quality?".format(
                        genome_name))
                if output_taxonomies_file:
                    output_taxonomies_fh.write(
                        '\t'.join([genome_name, 'Root; ' + '; '.join(taxonomy)]) + '\n')
            else:
                raise Exception("Unexpected taxonomy length found: {}. An example taxonomy is {}".format(
                    taxonomy, list(known_taxons)[0]))
        else:
            # Check that the taxonomy provided is a known taxonomy in the current metapackage.
            if not skip_taxonomy_check:
                taxonomy_to_check = re.sub(r'; {0,1}.__$', '', re.sub(r'; {0,1}.__;.*', '', ';'.join(taxonomy)))
                if taxonomy_to_check not in known_taxons:
                    raise Exception(
                        "The taxonomy {} for genome {} (originally {}) is not a known taxonomy in the current metapackage. This error might be due to a mismatch between the GTDB version used to generate the metapackage to be supplemented, and the GTDB version backing GTDBtk. For instance, are they both R214 or both R220?".format(
                            taxonomy_to_check, genome_name, taxonomy_str))
            if taxonomy[6] != 's__':
                logging.debug(
                    "Ignoring genome {} because it already has a species-level taxonomic assignment: {}".format(
                        genome_name, taxonomy))
                if output_taxonomies_file:
                    output_taxonomies_fh.write(
                        '\t'.join([genome_name, 'Root; ' + '; '.join(taxonomy)]) + '\n')
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

    if num_genomes_without_fasta > 0:
        if taxonomy_file:
            source = 'taxonomy file'
        else:
            source = 'GTDBtk output'
        logging.warning(
            "There were {} genomes in the {} that were not found in the list of genomes to be included. Ignoring these."
            .format(num_genomes_without_fasta, source))

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
                    key = orfm_utils.un_orfm_name(seq_name)
                    if key not in sequence_to_genome:
                        logging.debug(
                            "Skipping sequence {} because it is not in the sequence_to_genome map, interpreted to mean it matched the HMM but not the window specifically."
                            .format(seq_name))
                    else:
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
            no_further_euks=True,
        )

        logging.info("Regenerated {} ..".format(old_spkg.graftm_package_basename()))

    return new_spkg_path


def run_hmmsearch_on_one_genome(lock, data, matched_transcripts_fna, working_directory, hmmsearch_evalue, concatenated_hmms):
    total_num_transcripts = 0
    failure_genomes = 0
    num_transcriptomes = 0
    num_found_transcripts = 0

    (genome, tranp) = data

    hmmsearch_output = os.path.join(working_directory, 'hmmsearch_output',
                                    os.path.basename(genome) + '.hmmsearch')
    # hmmsearch_cmd = ['hmmsearch', '-o /dev/null', '-E 1e-20', '--tblout', '>(sed "s/  */\t/g" >', hmmsearch_output, ') --cpu', str(num_threads), concatenated_hmms, tranp.protein_fasta]
    hmmsearch_cmd = [
        'hmmsearch', '-o /dev/null', '-E ', hmmsearch_evalue, '--tblout', hmmsearch_output, concatenated_hmms, tranp.protein_fasta
    ]
    logging.debug("Running hmmsearch on {} ..".format(genome))
    # In the past have got strange errors about not extern.run not existing, so
    # using subprocess instead.
    subprocess.check_call(['bash','-c',' '.join(hmmsearch_cmd)])
    logging.debug("Ran hmmsearch on {} ..".format(genome))
    num_transcriptomes += 1

    qresults = list(SearchIO.parse(hmmsearch_output, 'hmmer3-tab'))
    found_hits_hmm = []
    found_hits_transcript = []
    found_hits_bitscore = []
    for qresult in qresults:
        for hit in qresult.hits:
            found_hits_hmm.append(qresult.id)
            found_hits_transcript.append(hit.id)
            found_hits_bitscore.append(hit.bitscore)
            total_num_transcripts += 1

    if total_num_transcripts == 0:
        logging.debug("No hits found for {}".format(genome))
        failure_genomes += 1
        return (total_num_transcripts, failure_genomes, num_transcriptomes, num_found_transcripts)

    df = pl.DataFrame({
        'hmm': found_hits_hmm,
        'transcript': found_hits_transcript,
        'bitscore': found_hits_bitscore
    })
    # Take best hit for each transcript
    df2 = df.sort('bitscore', descending=True).group_by('transcript').first()
    # Remove genes where there is >1 hit
    # If any HMM is associated with multiple genes, remove it as being potentially dubious
    ok_genes = list(df2.group_by('hmm').count().filter(pl.col('count') == 1).select('hmm'))[0]
    df3 = df2.filter(pl.col('hmm').is_in(ok_genes))
    logging.debug("After filtering, {} hits remain, from {} total hits".format(len(df3), len(df2)))

    # Create a set of the transcripts which hit. Lazily assumes that the
    # best hit in the hmmsearch output is the best hit by singlem pipe
    matched_transcript_ids = set(df3['transcript'])

    # Extract the transcripts belonging to these hits
    # Rename the transcripts, with genome and original ID separated by a ~
    num_printed = 0
    with open(tranp.transcript_fasta) as g:
        with lock:
            with open(matched_transcripts_fna, 'a') as f:
                for name, seq, _ in SeqReader().readfq(g):
                    if name in matched_transcript_ids:
                        genome_basename = remove_extension(os.path.basename(genome))
                        new_name = genome_basename + '‡' + name # Use ‡ to separate genome and original ID, must be kept in check with elsewhere in the code
                        print('>' + new_name + '\n' + seq + '\n', file=f)
                        num_printed += 1
    logging.debug("Printed {} transcripts for {}".format(num_printed, genome))
    if num_printed != len(matched_transcript_ids):
        logging.error("Expected to print {} transcripts for {}, but only printed {}".format(
            len(matched_transcript_ids), genome, num_printed))
    num_found_transcripts += num_printed

    return (total_num_transcripts, failure_genomes, num_transcriptomes, num_found_transcripts)


def gather_hmmsearch_results(num_threads, working_directory, old_metapackage, new_genome_transcripts_and_proteins,
                             hmmsearch_evalue):
    # Run hmmsearch using a concatenated set of HMMs from each graftm package in the metapackage
    # Create concatenated HMMs in working_directory/concatenated_alignment_hmms.hmm
    concatenated_hmms = os.path.join(working_directory, 'concatenated_alignment_hmms.hmm')
    num_hmms = 0
    with open(concatenated_hmms, 'w') as f:
        for spkg in old_metapackage.singlem_packages:
            with open(spkg.graftm_package().alignment_hmm_path()) as g:
                f.write(g.read())
                num_hmms += 1

    logging.info("Running HMMSEARCH using a concatenated set of {} HMMs".format(num_hmms))

    # Take the best hitting HMM for each transcript
    # mkdir working_directory/hmmsearch_output
    os.mkdir(os.path.join(working_directory, 'hmmsearch_output'))

    # Create a new file which is a concatenation of the transcripts we want to include
    # Use a lock to prevent race conditions since each worker writes this this
    matched_transcripts_fna = os.path.join(working_directory, 'matched_transcripts.fna')

    total_num_transcripts = 0
    total_failure_genomes = 0
    total_num_transcriptomes = 0
    total_num_found_transcripts = 0

    with Manager() as manager:
        lock = manager.Lock()

        # Since the multispawned function uses polars, we must use spawn
        # context, otherwise we get deadlock. See
        # https://pola-rs.github.io/polars/user-guide/misc/multiprocessing/#example
        with get_context('spawn').Pool(num_threads) as pool:
            map_result = pool.starmap(
                run_hmmsearch_on_one_genome,
                [(lock, data, matched_transcripts_fna, working_directory, hmmsearch_evalue, concatenated_hmms) for data in new_genome_transcripts_and_proteins.items()],
                chunksize=1)

            for (num_transcripts, failure_genomes, num_transcriptomes, num_found_transcripts) in map_result:
                total_num_transcripts += num_transcripts
                total_failure_genomes += failure_genomes
                total_num_transcriptomes += num_transcriptomes
                total_num_found_transcripts += num_found_transcripts

    logging.info(
        "Ran hmmsearch on {} genomes, finding {} marker genes. {} were excluded based on 2+ copy number.".format(
            total_num_transcriptomes, total_num_found_transcripts, total_num_transcripts - total_num_found_transcripts))
    if failure_genomes > 0:
        logging.warning("hmmsearch failed to find any marker genes for {} genomes".format(failure_genomes))

    if total_num_found_transcripts == 0:
        logging.error("Unable to find any marker genes from the provided genomes, so cannot generate a new metapackage")
        return None
    else:
        return matched_transcripts_fna


def generate_new_metapackage(num_threads, working_directory, old_metapackage, new_genome_fasta_files,
                             new_taxonomies_file, new_genome_transcripts_and_proteins, hmmsearch_evalue,
                             checkm2_quality_file, no_taxon_genome_lengths, new_taxonomy_database_name,
                             new_taxonomy_database_version):

    # Add the new genome data to each singlem package
    # For each package, the unaligned seqs are in the graftm package,
    # Taxonomy/seqinfo give the taxonomy of target and euk sequences

    genome_to_taxonomy = {}
    with open(new_taxonomies_file) as f:
        for row in csv.reader(f, delimiter='\t'):
            genome = FastaNameToSampleName.fasta_to_name(row[0])
            genome_to_taxonomy[os.path.basename(genome)] = row[1]
    if len(genome_to_taxonomy) != len(new_genome_fasta_files) or len(genome_to_taxonomy) == 0:
        raise Exception("Unexpected number of genomes with taxonomy vs. fasta files")

    # Recalculate genome sizes for each species and higher level taxon
    if not no_taxon_genome_lengths:
        logging.info("Recalculating genome sizes for each species and higher level taxon ..")
        new_genome_sizes = recalculate_genome_sizes(
            old_metapackage,
            genome_to_taxonomy,
            checkm2_quality_file,
            new_genome_fasta_files)
        taxon_genome_lengths_tmpfile = tempfile.NamedTemporaryFile(mode='w', prefix='taxon_genome_lengths_', suffix='.tsv')
        new_genome_sizes.write_csv(taxon_genome_lengths_tmpfile.name, separator='\t')
        logging.info("Finished recalculating genome sizes")

    logging.info("Gathering OTUs from new genomes ..")
    matched_transcripts_fna = gather_hmmsearch_results(num_threads, working_directory, old_metapackage,
                                                       new_genome_transcripts_and_proteins, hmmsearch_evalue)
    if matched_transcripts_fna is None:
        # No transcripts were found, so no new OTUs to add
        return None

    # Run singlem pipe --forward on that, not assigning taxonomy - we know taxonomy from gtdbtk
    logging.info("Running singlem pipe --forward on gathered transcripts ..")
    # pipe messes up the tempdir, so cache and set back
    original_tmpdir = tempfile.gettempdir()
    new_genomes_otu_table_unassigned = SearchPipe().run_to_otu_table(
        sequences=[matched_transcripts_fna],
        metapackage_object=old_metapackage,
        threads=num_threads,
        diamond_prefilter=True,
        assign_taxonomy=False)  # Maybe could remove off-targets here?
    tempfile.tempdir = original_tmpdir

    # Assign taxonomy. When multiple genomes have the same OTU, add each genome
    # as a separate OTU.
    taxonomy_field = ArchiveOtuTable.TAXONOMY_FIELD_INDEX
    sample_name_field = ArchiveOtuTable.SAMPLE_ID_FIELD_INDEX
    read_name_field = ArchiveOtuTable.READ_NAME_FIELD_INDEX
    new_genomes_otu_table = ArchiveOtuTable()
    sequence_to_genome = {}

    bad_chars = ['‡']
    # raise if there are bad chars in genome names
    for g in genome_to_taxonomy.keys():
        if any(char in g for char in bad_chars):
            raise ValueError("Because they are used as delimiters within SingleM, the '‡' characters cannot be used in genome names, or transcript names. The offending genome name is: {}".format(g))

    for otu in new_genomes_otu_table_unassigned:
        for read_name in otu.data[6]:
            # Raise it there are >1 of any single bad char in the read name
            for bad in bad_chars:
                if read_name.count(bad) > 1:
                    raise ValueError("Because they are used as delimiters within SingleM, the '‡' characters cannot be used in genome names, or transcript names. The offending read name is: {}".format(read_name))
            genome, _ = read_name.split('‡')
            sequence_to_genome[read_name] = FastaNameToSampleName.fasta_to_name(genome)

            data2 = otu.data.copy()
            data2[sample_name_field] = genome
            data2[taxonomy_field] = genome_to_taxonomy[FastaNameToSampleName.fasta_to_name(genome)]
            data2[read_name_field] = [read_name]
            entry = OtuTableEntry()
            entry.data = data2
            new_genomes_otu_table.add([entry])

    # Create SDB for final metapackage
    logging.info("Creating sdb to include in metapackage ..")
    new_metapackage_sdb_path = os.path.join(working_directory, 'new_metapackage.sdb')

    # Add old metapackage's OTUs.
    db_otu_table = SequenceDatabaseOtuTable(old_metapackage.nucleotide_sdb())
    all_new_and_old_otu_tables = OtuTableCollection()
    all_new_and_old_otu_tables.add_otu_table_object(new_genomes_otu_table)
    all_new_and_old_otu_tables.add_otu_table_object(db_otu_table)

    SequenceDatabase().create_from_otu_table(
        new_metapackage_sdb_path,
        all_new_and_old_otu_tables,
        num_threads=num_threads,
    )

    new_metapackage_path = os.path.join(working_directory, 'new_metapackage.mpkg')

    # For each spkg in the old mpkg, create a new spkg
    to_process = [(working_directory, spkg, matched_transcripts_fna, sequence_to_genome, genome_to_taxonomy)
                  for spkg in old_metapackage.singlem_packages]
    if num_threads > 1:
        new_spkg_paths = Pool(num_threads).map(generate_new_singlem_package, to_process)
    else:
        new_spkg_paths = list([generate_new_singlem_package(x) for x in to_process])
    
    # Check if old_metapackage uses avg_num_genes_per_species
    if old_metapackage.avg_num_genes_per_species() is not None:
        # Will need to recalculate this for the new metapackage
        calculate_average_num_genes_per_species = True
    else:
        calculate_average_num_genes_per_species = False

    # Create a new metapackage from the singlem packages
    logging.info("Creating new metapackage ..")
    Metapackage.generate(singlem_packages=new_spkg_paths,
                         nucleotide_sdb=new_metapackage_sdb_path,
                         prefilter_diamond_db=old_metapackage.prefilter_db_path(),
                         output_path=new_metapackage_path,
                         threads=num_threads,
                         prefilter_clustering_threshold=None,
                         taxon_genome_lengths=taxon_genome_lengths_tmpfile.name if not no_taxon_genome_lengths else None,
                         taxonomy_database_name=new_taxonomy_database_name,
                         taxonomy_database_version=new_taxonomy_database_version,
                         diamond_prefilter_performance_parameters=old_metapackage.diamond_prefilter_performance_parameters(),
                         diamond_taxonomy_assignment_performance_parameters=old_metapackage.diamond_taxonomy_assignment_performance_parameters(),
                         makeidx_sensitivity_params=old_metapackage.makeidx_sensitivity_params(),
                         calculate_average_num_genes_per_species=calculate_average_num_genes_per_species)
    logging.info("New metapackage created at {}".format(new_metapackage_path))

    if not no_taxon_genome_lengths:
        taxon_genome_lengths_tmpfile.close()

    return new_metapackage_path


class TranscriptsAndProteins:
    def __init__(self):
        self.transcript_fasta = None
        self.protein_fasta = None
        self.translation_table = None


def recalculate_genome_sizes(
        old_metapackage,
        genome_to_taxonomy,
        checkm2_quality_file,
        new_genome_fasta_files):

    # Read old taxon lengths, just recording the species-level ones
    old_taxon_lengths = {}
    old_metapackage_taxon_lengths = old_metapackage.taxon_genome_lengths()
    if old_metapackage_taxon_lengths is None:
        raise Exception("Supplement mode requires an input metapackage which has included taxon genome lengths")

    for row in pl.DataFrame(old_metapackage_taxon_lengths).rows(named=True):
        # rank    genome_size
        rank = row['rank']
        if rank.startswith('s__'):
            old_taxon_lengths[row['rank']] = row['genome_size']
    logging.info("Read {} species-level taxon lengths from old metapackage".format(len(old_taxon_lengths)))

    # Map previous rank IDs, which are in the form of "s__xxx", to fully defined taxons
    logging.info("Reading taxonomies from old metapackage ..")
    old_taxonomies = old_metapackage.get_all_taxonomy_strings()
    logging.debug("Read {} full taxonomy strings from old metapackage".format(len(old_taxonomies)))
    species_to_full = {}
    logging.info("Creating mapping of species-level taxon IDs to full taxonomies ..")
    for taxonomy in old_taxonomies:
        last = taxonomy.split(';')[-1].strip()
        species_to_full[last] = taxonomy

    # In really rare instances (e.g. spire_mag_01111928) there are no reads
    # included in the metapackage because all are filtered as being double. But
    # those ones are recorded in the old_taxon_lengths. Avoid these ones.
    logging.debug("Filtering out species from the old metapackage which have no entries in singlem profile ..")
    species_to_add = []
    num_excluded = 0
    for s in old_taxon_lengths.keys():
        if s in species_to_full:
            species_to_add.append(s)
        else:
            logging.debug("Species {} was not found in the old metapackage's full taxonomies, suggesting potentially poor quality".format(s))
            num_excluded += 1
    if num_excluded > 0:
        logging.info("Excluded {} species from the old metapackage because they were not found in the full taxonomies".format(num_excluded))
    logging.debug("Creating gc dataframe, which has been known to cause issues with polars ..")
    gc = pl.DataFrame(
        {
            'species': species_to_add,
            'gtdb_taxonomy': [species_to_full[s] for s in species_to_add],
            'genome_size': [old_taxon_lengths[s] for s in species_to_add]
        }
    )

    # Add new genomes to the species-level taxon lengths
    # First read checkm2 estimate
    logging.info("Reading checkm2 quality file ..")
    checkm2 = CheckM2(checkm2_quality_file)

    logging.info("Reading genome FASTA files to calculate genome lengths ..")
    new_taxon_lengths = {}
    all_stats = checkm2.get_all_stats()
    for genome_fasta in new_genome_fasta_files:
        checkm_stats = all_stats[FastaNameToSampleName.fasta_to_name(genome_fasta)]

        corrected_length = GenomeSizes.corrected_genome_size(checkm_stats.genome_size, checkm_stats.completeness, checkm_stats.contamination)
        taxonomy = genome_to_taxonomy[FastaNameToSampleName.fasta_to_name(genome_fasta)]

        new_taxon_lengths[taxonomy] = corrected_length

    gtdb_taxonomies = new_taxon_lengths.keys()
    news = pl.DataFrame({
        'gtdb_taxonomy': gtdb_taxonomies,
        'genome_size': [new_taxon_lengths[n] for n in gtdb_taxonomies]
    })

    gc = pl.concat([
        gc.select(['gtdb_taxonomy', 'genome_size']),
        news
    ])

    # Calculate the genome lengths for each taxon level
    logging.info("Calculating inferred genome lengths for each taxon level ..")
    all_rank_genome_sizes = GenomeSizes.calculate_rank_genome_sizes(gc['gtdb_taxonomy'], gc['genome_size'])

    # Return dataframe of taxon lengths in same shape as metapackage.generate expects
    return all_rank_genome_sizes.select(['rank', 'genome_size'])


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
    for x in transcripts_and_proteins:
        logging.debug("Translation table for {} was {}".format(x.protein_fasta, x.translation_table))
    logging.info("Prodigal finished. Found {} with translation table 4, {} with translation table 11".format(
        num_translation_table4, num_translation_table11))

    to_return = {}
    for genome_fasta, transcripts_and_proteins in zip(new_genome_fasta_files, transcripts_and_proteins):
        to_return[genome_fasta] = transcripts_and_proteins
    return to_return


def read_gene_definitions(gene_definitions):
    df = pl.read_csv(gene_definitions, separator='\t')
    if df.columns != ['genome_fasta', 'transcript_fasta', 'protein_fasta']:
        raise Exception("Gene definitions file must have headings 'genome_fasta', 'transcript_fasta', 'protein_fasta'. Found: {}".format(df.columns))

    to_return = {}
    for row in df.rows(named=True):
        genome = os.path.abspath(row['genome_fasta'])
        to_return[genome] = TranscriptsAndProteins()
        to_return[genome].transcript_fasta = os.path.abspath(row['transcript_fasta'])
        to_return[genome].protein_fasta = os.path.abspath(row['protein_fasta'])
        # Test they exist here, because otherwise they silently fail later
        if not os.path.exists(to_return[genome].transcript_fasta):
            raise Exception("Transcript fasta file {} does not exist".format(to_return[genome].transcript_fasta))
        if not os.path.exists(to_return[genome].protein_fasta):
            raise Exception("Protein fasta file {} does not exist".format(to_return[genome].protein_fasta))
        if not os.path.exists(genome):
            raise Exception("Genome fasta file {} does not exist".format(genome))
    logging.info("Read {} gene definitions from {}".format(len(to_return), gene_definitions))
    return to_return


def dereplicate_genomes_with_galah(**kwargs):
    threads = kwargs.pop('threads')
    genomes_to_dereplicate = kwargs.pop('genomes_to_dereplicate')
    checkm2_quality_file = kwargs.pop('checkm2_quality_file')
    if len(kwargs) > 0:
        raise Exception("Unexpected arguments detected: %s" % kwargs)

    quality_arg = ""
    if checkm2_quality_file is None:
        logging.warning("Galah is being run without a checkm2 quality file. This is not recommended, because then galah may choose non-optimal genomes as representative.")
    else:
        quality_arg = "--checkm2-quality-report {}".format(checkm2_quality_file)
    logging.info("Running galah to dereplicate {} genomes ..".format(len(genomes_to_dereplicate)))
    new_genome_paths = extern.run("galah cluster {} --threads {} --genome-fasta-list /dev/stdin --ani 95 --output-representative-list /dev/stdout".format(
        quality_arg, threads), stdin="\n".join(genomes_to_dereplicate)).splitlines()
    logging.info("After dereplication, {} genomes remain".format(len(new_genome_paths)))
    return new_genome_paths


class Supplementor:

    def supplement(self, **kwargs):
        new_genome_fasta_files = kwargs.pop('new_genome_fasta_files')
        new_genome_fasta_files_list = kwargs.pop('new_genome_fasta_files_list')
        new_taxonomies = kwargs.pop('new_taxonomies')
        input_metapackage = kwargs.pop('input_metapackage')
        output_metapackage = kwargs.pop('output_metapackage')
        threads = kwargs.pop('threads')
        pplacer_threads = kwargs.pop('pplacer_threads')
        predefined_working_directory = kwargs.pop('working_directory')
        gtdbtk_output_directory = kwargs.pop('gtdbtk_output_directory')
        taxonomy_file = kwargs.pop('taxonomy_file')
        output_taxonomies = kwargs.pop('output_taxonomies')
        checkm2_quality_file = kwargs.pop('checkm2_quality_file')
        checkm2_min_completeness = kwargs.pop('checkm2_min_completeness')
        checkm2_max_contamination = kwargs.pop('checkm2_max_contamination')
        no_quality_filter = kwargs.pop('no_quality_filter')
        hmmsearch_evalue = kwargs.pop('hmmsearch_evalue')
        # no_dereplication=args.no_dereplication,
        # dereplicate_with_galah=args.dereplicate_with_galah,
        no_dereplication = kwargs.pop('no_dereplication')
        dereplicate_with_galah = kwargs.pop('dereplicate_with_galah')
        skip_taxonomy_check = kwargs.pop('skip_taxonomy_check')
        no_taxon_genome_lengths = kwargs.pop('no_taxon_genome_lengths')
        gene_definitions = kwargs.pop('gene_definitions')
        ignore_taxonomy_database_incompatibility = kwargs.pop('ignore_taxonomy_database_incompatibility')
        new_taxonomy_database_name = kwargs.pop('new_taxonomy_database_name')
        new_taxonomy_database_version = kwargs.pop('new_taxonomy_database_version')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)

        os.environ["POLARS_MAX_THREADS"] = str(1)
        import polars as pl

        if new_genome_fasta_files is None and new_genome_fasta_files_list is None:
            raise Exception("Must provide either --new-genome-fasta-files or --new-genome-fasta-files-list")
        elif new_genome_fasta_files is not None and new_genome_fasta_files_list is not None:
            raise Exception("Must provide either --new-genome-fasta-files or --new-genome-fasta-files-list, not both")
        elif new_genome_fasta_files_list is not None:
            for file_list in new_genome_fasta_files_list:
                with open(file_list) as f:
                    new_genome_fasta_files = [os.path.abspath(x.strip()) for x in f.readlines()]
        else:
            new_genome_fasta_files = [os.path.abspath(x) for x in new_genome_fasta_files]

        if input_metapackage is None:
            old_metapackage = Metapackage.acquire_default()
        else:
            old_metapackage = Metapackage.acquire(input_metapackage)
        if ignore_taxonomy_database_incompatibility:
            logging.warning("Ignoring potential taxonomy database incompatibility as specified. You are responsible for ensuring that any new genomes added to the metapackage are from species not currently in the input metapackage.")
        elif old_metapackage.version < 5:
            raise Exception("The input metapackage was generated using an older version of SingleM. You are responsible for ensuring that any new genomes added to the metapackage are from species not currently in the input metapackage, since these old versions did not track which genome database was used.")
        elif old_metapackage.taxonomy_database_name() != GTDB_DATABASE_NAME:
            raise Exception("The input metapackage was not generated using the GTDB taxonomy database, but was recorded as '{}' version '{}'. We are halting here to alert you to the fact that supplementing a metapackage that is already supplemented makes you responsible for ensuring that there aren't any species which have multiple representatives e.g. one from the first supplement and one from the second. If you are confident there are no species like this, rerun supplement mode with --ignore-taxonomy-database-incompatibility.".format(
                old_metapackage.taxonomy_database_name(),
                old_metapackage.taxonomy_database_version()))
        
        if new_taxonomy_database_name is None:
            new_taxonomy_database_name = CUSTOM_TAXONOMY_DATABASE_NAME

        with tempfile.TemporaryDirectory() as working_directory:
            if predefined_working_directory:
                working_directory = predefined_working_directory
            if not os.path.exists(working_directory):
                os.mkdir(working_directory)
            # Set tmpdir to the working directory/tmp so that we don't fill
            # up /tmp when working directory is specified, and we keep
            # things cleaner. See https://github.com/wwood/singlem/issues/168
            tempfile_directory = os.path.join(working_directory, 'tmp')
            os.mkdir(tempfile_directory)
            tempfile.tempdir = tempfile_directory

            # Quality filter the new genomes
            num_before_quality_filter = len(new_genome_fasta_files)
            excluded_genomes = set()
            if no_quality_filter:
                logging.info("Skipping quality filtering of new genomes")
            elif checkm2_quality_file:
                logging.info("Quality filtering new genomes using checkm2 quality file: %s" % checkm2_quality_file)
                checkm = CheckM2(checkm2_quality_file)
                new_genome_fasta_files_without_file_extension = set(
                    checkm.genomes_of_sufficient_quality(checkm2_min_completeness, checkm2_max_contamination))
                # Check the new genomes are in the checkm2 quality file. Cache the names as a set for fast lookup
                checkm_names = set(checkm.names())
                for g in new_genome_fasta_files:
                    if remove_extension(g) not in checkm_names:
                        raise Exception("Genome %s not found in checkm2 quality file" % g)
                new_genome_fasta_files2 = []
                for x in new_genome_fasta_files:
                    if remove_extension(x) in new_genome_fasta_files_without_file_extension:
                        new_genome_fasta_files2.append(x)
                    else:
                        excluded_genomes.add(x)
                new_genome_fasta_files = new_genome_fasta_files2
                logging.info(
                    "Removed {} genomes of insufficient quality, leaving {} genomes to be added to the metapackage".
                    format(num_before_quality_filter - len(new_genome_fasta_files), len(new_genome_fasta_files)))
            else:
                raise Exception("Must provide either --no-quality-filter or --checkm2-quality-file")

            # Dereplicate the new genomes if required
            if no_dereplication:
                pass
            elif dereplicate_with_galah:
                logging.info("Dereplicating new genomes using Galah ..")
                new_genome_fasta_files = dereplicate_genomes_with_galah(
                    threads=threads, genomes_to_dereplicate=new_genome_fasta_files, checkm2_quality_file=checkm2_quality_file)
            else:
                raise Exception("Must provide either --no-dereplication or --dereplicate-with-galah")

            # Generate faa and transcript fna files for the new genomes
            if gene_definitions:
                logging.info("Reading transcript and protein files from gene definitions file ..")
                genome_transcripts_and_proteins = read_gene_definitions(gene_definitions)
            else:
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
                    taxonomy_file=taxonomy_file,
                    pplacer_threads=pplacer_threads,
                    output_taxonomies_file=output_taxonomies,
                    excluded_genomes=excluded_genomes,
                    old_metapackage=old_metapackage,
                    skip_taxonomy_check=skip_taxonomy_check)
                # Remove genomes that were excluded by not being novel at the species level
                # and check that we have all the genomes we expect in the gene calls
                genome_transcripts_and_proteins1 = {}
                for ng in new_genome_fasta_files:
                    if ng in genome_transcripts_and_proteins:
                        genome_transcripts_and_proteins1[ng] = genome_transcripts_and_proteins[ng]
                    else:
                        raise Exception("Transcript and protein sequences not found for genome {} ".format(ng))
                genome_transcripts_and_proteins = genome_transcripts_and_proteins1

            # Run the genomes through pipe with genome fasta input to identify the new sequences
            logging.info("Generating new SingleM packages and metapackage ..")
            new_metapackage = generate_new_metapackage(threads, working_directory, old_metapackage,
                                                       new_genome_fasta_files, taxonomy_file,
                                                       genome_transcripts_and_proteins, hmmsearch_evalue,
                                                       checkm2_quality_file, no_taxon_genome_lengths,
                                                       new_taxonomy_database_name, new_taxonomy_database_version)

            if new_metapackage is not None:
                logging.info("Copying generated metapackage to {}".format(output_metapackage))
                shutil.copytree(new_metapackage, output_metapackage)
                logging.info("Finished supplement.")

# if __name__ == "__main__":
#     import logging
#     logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    
#     old_metapackage = Metapackage.acquire('/home/woodcrob/git/singlem/db/S4.2.2.GTDB_r214.metapackage_20240502.smpkg.zb/payload_directory')
#     recalculate_genome_sizes(
#         old_metapackage,
#         None,
#         None,
#         None)
