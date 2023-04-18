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

import argparse
import logging
import sys
import os
import csv
import shutil
from multiprocessing import Pool
import tempfile
from random import shuffle

import extern

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

from singlem.metapackage import Metapackage
from singlem.archive_otu_table import ArchiveOtuTable
from singlem.sequence_classes import SeqReader
from singlem.singlem import OrfMUtils
from singlem.condense import CondensedCommunityProfile
from singlem.otu_table_collection import OtuTableCollection


def generate_new_singlem_package(myargs):
    (working_directory,
    old_spkg,
    matched_transcript_path,
    sequence_to_genome,
    genome_to_taxonomy) = myargs

    logging.info("Preparing GraftM create inputs for {}".format(old_spkg.graftm_package_basename()))

    # Change working dir so that each spkg is isolated
    working_directory = os.path.join(working_directory, 'create-'+old_spkg.graftm_package_basename())
    os.makedirs(working_directory)

    # Create new taxonomy
    # new_gpkg_path = os.path.join(working_directory, old_spkg.graftm_package_basename() + '.gpkg')
    new_spkg_path = os.path.join(working_directory, old_spkg.graftm_package_basename() + '.spkg')

    # Run graftM graft to get the protein sequences from the transcripts which matched the HMM
    graft_output_folder = os.path.join(working_directory, 'graft_output')
    extern.run('graftM graft --graftm_package {} --forward {} --output {} --search_only'.format(
        old_spkg.graftm_package_path(),
        matched_transcript_path,
        graft_output_folder))

    # Work out the sequence prefix
    with open(old_spkg.graftm_package().unaligned_sequence_database_path()) as f:
        first_line = f.readline()
        if first_line.startswith('>'):
            sequence_prefix = first_line[1:].split('~')[0]
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
                    new_sequences_fasta_fh.write('>' + sequence_prefix+'~'+seq_name + '\n' + seq + '\n')
                    f.write(sequence_prefix+'~'+seq_name + '\t' + genome_to_taxonomy[genome] + '\n')

        # Log how many of the genomes matched this marker
        logging.info("Added {} sequences from {} genomes".format(
            len(matched_sequences),
            len(genome_to_taxonomy)))

        # Run singlem regenerate to create a new spkg
        logging.info("Regenerating {} ..".format(old_spkg.graftm_package_basename()))
        cmd = 'singlem regenerate --no-further-euks --input-singlem-package {} --output-singlem-package {} --input-sequences {} --input-taxonomy {}'.format(
            old_spkg.base_directory(),
            new_spkg_path,
            new_sequences_fasta,
            taxonomy_file)
        extern.run(cmd)

        logging.info("Regenerated {} ..".format(old_spkg.graftm_package_basename()))

    return new_spkg_path


def generate_new_metapackage(
    num_threads,
    working_directory,
    old_metapackage_path,
    new_genome_fasta_files,
    new_taxonomies_file):

    # Add the new genome data to each singlem package
    # For each package, the unaligned seqs are in the graftm package,
    # Taxonomy/seqinfo give the taxonomy of target and euk sequences
    old_metapackage = Metapackage.acquire(old_metapackage_path)

    genome_otu_table = os.path.join(working_directory, 'genome_otu_table.json')

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
    extern.run('singlem pipe --genome-fasta-files {} --archive-otu-table {} --metapackage {} --no-assign-taxonomy'.format(
        ' '.join(new_genome_fasta_files),
        genome_otu_table,
        old_metapackage_path))

    # Create SDB for final metapackage
    logging.info("Creating sdb to include in metapackage ..")
    new_metapackage_sdb_path = os.path.join(working_directory, 'new_metapackage.sdb')
    with tempfile.NamedTemporaryFile(prefix='supplement_metapackage_otu_table_renamed') as f:
        # Change the taxonomy to be correct
        with open(genome_otu_table) as g:
            otus = ArchiveOtuTable.read(g)

            for otu in otus:
                otu.taxonomy = genome_to_taxonomy[otu.sample_name]
                f.write((str(otu)+"\n").encode())

        f.flush()
        extern.run('singlem makedb --otu-tables {} <(singlem query --dump --db {}) --db {} --sequence-database-methods smafa-naive --threads {}'.format(
            f.name, old_metapackage.nucleotide_sdb().base_directory, new_metapackage_sdb_path, num_threads
        ))

    # Dump matched transcript sequences to a file for graftm graft
    matched_transcript_path = os.path.join(working_directory, 'matched_transcripts.fasta')
    extern.run('singlem summarise --input-archive-otu-table {} --unaligned-sequences-dump-file {}'.format(
        genome_otu_table,
        matched_transcript_path))
    with open(genome_otu_table) as f:
        genome_archive_otu = ArchiveOtuTable.read(f)
    sequence_to_genome = {}
    with open(matched_transcript_path, 'w') as f:
        for otu in genome_archive_otu:
            read_names = otu.read_names()
            seqs = otu.read_unaligned_sequences()
            for seq, read_name in zip(seqs, read_names):
                if read_name in sequence_to_genome:
                    raise Exception("Duplicate sequence name: {}".format(read_name))
                sequence_to_genome[read_name] = otu.sample_name
                f.write('>{}\n{}\n'.format(read_name, seq))


    new_metapackage_path = os.path.join(working_directory, 'new_metapackage.mpkg')
    
    # For each spkg in the old mpkg, create a new spkg
    to_process = [(working_directory, spkg, matched_transcript_path, sequence_to_genome, genome_to_taxonomy) for spkg in old_metapackage.singlem_packages]
    new_spkg_paths = Pool(num_threads).map(generate_new_singlem_package, to_process)

    # Create a new metapackage from the singlem packages
    logging.info("Creating new metapackage ..")
    extern.run(
        'singlem metapackage --nucleotide-sdb {} --singlem-packages {} --metapackage {} --prefilter-diamond-db {} --no-taxon-genome-lengths'.format(
            new_metapackage_sdb_path, ' '.join(new_spkg_paths), new_metapackage_path, old_metapackage.prefilter_db_path()))
    logging.info("New metapackage created at {}".format(new_metapackage_path))
    return new_metapackage_path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser = argparse.ArgumentParser(parents=[parent_parser])
    subparsers = parser.add_subparsers(title="Sub-commands", dest='subparser_name')

    metapackage_description = 'Create a new metapackage from a vanilla one plus new genomes'
    metapackage_parser = subparsers.add_parser('metapackage')
    metapackage_parser.add_argument('--new-genome-fasta-files', help='FASTA files of new genomes', nargs='+', required=True)
    metapackage_parser.add_argument('--new-taxonomies', help='newline separated file containing taxonomies of new genomes (path<TAB>taxonomy). Must be fully specified to species level.', required=True)
    metapackage_parser.add_argument('--input-metapackage', help='metapackage to build upon', required=True)
    metapackage_parser.add_argument('--output-metapackage', help='output metapackage', required=True)
    metapackage_parser.add_argument('--threads', help='parallelisation', type=int, default=1)
    args = parser.parse_args()

    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        logging.getLogger('nmslib').setLevel(logging.ERROR)
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
        logging.getLogger('nmslib').setLevel(logging.WARN)
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

    if args.subparser_name == 'metapackage':
        with tempfile.TemporaryDirectory() as working_directory:
            # working_directory = '/tmp/ben_find_tmp'
            # os.mkdir(working_directory)

            # Run the genomes through pipe with genome fasta input to identify the new sequences
            logging.info("Generating new SingleM packages and metapackage ..")
            new_metapackage = generate_new_metapackage(
                args.threads,
                working_directory,
                args.input_metapackage,
                args.new_genome_fasta_files,
                args.new_taxonomies)

            logging.info("Copying generated metapackage to {}".format(args.output_metapackage))
            shutil.copytree(new_metapackage, args.output_metapackage)

    else:
        raise Exception("Unexpected subparser name: {}".format(args.subparser_name))

    logging.info("Done!")