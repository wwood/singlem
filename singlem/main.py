#!/usr/bin/env python3

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2015-2024"
__credits__ = ["Ben Woodcroft", "Samuel Aroney", "Rossen Zhao"]
__license__ = "GPL3+"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import logging
import sys
import os
import gzip
import tempfile
import json

from bird_tool_utils import *
from bird_tool_utils.people import *

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

import singlem
import singlem.pipe as pipe
from singlem.pipe import SearchPipe
from singlem.condense import Condenser
from singlem.metapackage import DATA_ENVIRONMENT_VARIABLE, CUSTOM_TAXONOMY_DATABASE_NAME
from singlem import OTU_TABLE_OUTPUT_FORMAT, ARCHIVE_TABLE_OUTPUT_FORMAT

from singlem.condense import DEFAULT_MIN_TAXON_COVERAGE as CONDENSE_DEFAULT_MIN_TAXON_COVERAGE
from singlem.condense import DEFAULT_GENOME_MIN_TAXON_COVERAGE as CONDENSE_DEFAULT_GENOME_MIN_TAXON_COVERAGE
from singlem.condense import DEFAULT_TRIM_PERCENT as CONDENSE_DEFAULT_TRIM_PERCENT

DEFAULT_WINDOW_SIZE = 60
SPECIES_LEVEL_AVERAGE_IDENTITY = float(DEFAULT_WINDOW_SIZE - SearchPipe.DEFAULT_MAX_SPECIES_DIVERGENCE) / DEFAULT_WINDOW_SIZE


def seqs(args):
    from singlem.sequence_classes import SeqReader as SingleMSeqReader
    from singlem.metagenome_otu_finder import MetagenomeOtuFinder

    if args.alignment_type == 'aa':
        is_protein_alignment = True
    elif args.alignment_type == 'dna':
        is_protein_alignment = False
    else:
        raise Exception("Unexpected alignment type '%s'" % args.alignment_type)

    # Read in the fasta Alignment
    protein_alignment = SingleMSeqReader().alignment_from_alignment_file(args.alignment)
    logging.info("Read in %i aligned protein sequences e.g. %s %s" % (
        len(protein_alignment),
        protein_alignment[0].name,
        protein_alignment[0].seq))

    best_position = MetagenomeOtuFinder().find_best_window(
        protein_alignment,
        args.window_size,
        is_protein_alignment,
        args.hmm)
    # Report the best position as a 1-based index
    best_position += 1

    logging.info("Found best start position %i" % best_position)
    print(best_position)

# Make pipe argument functions here so the code can be re-used between pipe and renew
def add_common_pipe_arguments(argument_group, extra_args=False):
    if extra_args:
        sequence_input_group = argument_group.add_mutually_exclusive_group(required=True)
        # Keep parity of these arguments with the 'read_fraction' command
        sequence_input_group.add_argument('-1','--forward','--reads','--sequences',
                                    nargs='+',
                                    metavar='sequence_file',
                                    help='nucleotide read sequence(s) (forward or unpaired) to be searched. Can be FASTA or FASTQ format, GZIP-compressed or not, short or long (but Nanopore >=10.4.1 or PacBio HiFi reads recommended).')
        argument_group.add_argument('-2', '--reverse',
                                    nargs='+',
                                    metavar='sequence_file',
                                    help='reverse reads to be searched. Can be FASTA or FASTQ format, GZIP-compressed or not.')
        sequence_input_group.add_argument('-f', '--genome-fasta-files',
                                    nargs='+',
                                    metavar='PATH',
                                    help='Path(s) to genome FASTA files. These are processed like input given with --forward, but use higher default values for --min-taxon-coverage and --min-orf-length.')
        sequence_input_group.add_argument('-d', '--genome-fasta-directory',
                                    metavar='PATH',
                                    help='Directory containing genome FASTA files. Treated identically to --forward input with higher default values for --min-taxon-coverage and --min-orf-length.')
        sequence_input_group.add_argument('--genome-fasta-list',
                                    metavar='PATH',
                                    help='File containing genome FASTA paths, one per line. Behaviour matches --forward with higher default values for --min-taxon-coverage and --min-orf-length.')
        argument_group.add_argument('-x', '--genome-fasta-extension',
                                    metavar='EXT',
                                    help='File extension of genomes in the directory specified with -d/--genome-fasta-directory. [default: fna]',
                                    default='fna')
    argument_group.add_argument('-p', '--taxonomic-profile', metavar='FILE', help="output a 'condensed' taxonomic profile for each sample based on the OTU table. Taxonomic profiles output can be further converted to other formats using singlem summarise.")
    argument_group.add_argument('--taxonomic-profile-krona', metavar='FILE', help="output a 'condensed' taxonomic profile for each sample based on the OTU table")
    argument_group.add_argument('--otu-table', metavar='filename', help='output OTU table')
    current_default = pipe.DEFAULT_THREADS
    argument_group.add_argument('--threads', type=int, metavar='num_threads', help='number of CPUS to use [default: %i]' % current_default, default=current_default)
    current_default = SearchPipe.DEFAULT_TAXONOMY_ASSIGNMENT_METHOD
    argument_group.add_argument(
        '--assignment-method', '--assignment_method',
        choices=(
                pipe.SMAFA_NAIVE_THEN_DIAMOND_ASSIGNMENT_METHOD,
                pipe.SCANN_NAIVE_THEN_DIAMOND_ASSIGNMENT_METHOD,
                pipe.ANNOY_THEN_DIAMOND_ASSIGNMENT_METHOD,
                pipe.SCANN_THEN_DIAMOND_ASSIGNMENT_METHOD,
                pipe.DIAMOND_ASSIGNMENT_METHOD,
                pipe.DIAMOND_EXAMPLE_BEST_HIT_ASSIGNMENT_METHOD,
                pipe.ANNOY_ASSIGNMENT_METHOD,
                pipe.PPLACER_ASSIGNMENT_METHOD),
        help='Method of assigning taxonomy to OTUs and taxonomic profiles [default: %s]\n\n' % (current_default) +
            table_roff([
                ["Method", "Description"],
                [pipe.SMAFA_NAIVE_THEN_DIAMOND_ASSIGNMENT_METHOD, "Search for the most similar window sequences <= 3bp different using a brute force algorithm (using the smafa implementation) over all window sequences in the database, and if none are found use DIAMOND blastx of all reads from each OTU."],
                [pipe.SCANN_NAIVE_THEN_DIAMOND_ASSIGNMENT_METHOD, "Search for the most similar window sequences <= 3bp different using a brute force algorithm over all window sequences in the database, and if none are found use DIAMOND blastx of all reads from each OTU."],
                [pipe.ANNOY_THEN_DIAMOND_ASSIGNMENT_METHOD, "Same as {}, except search using ANNOY rather than using brute force. Requires a non-standard metapackage.".format(pipe.SCANN_NAIVE_THEN_DIAMOND_ASSIGNMENT_METHOD)],
                [pipe.SCANN_THEN_DIAMOND_ASSIGNMENT_METHOD, "Same as {}, except search using SCANN rather than using brute force. Requires a non-standard metapackage.".format(pipe.SCANN_NAIVE_THEN_DIAMOND_ASSIGNMENT_METHOD)],
                [pipe.DIAMOND_ASSIGNMENT_METHOD, "DIAMOND blastx best hit(s) of all reads from each OTU."],
                [pipe.DIAMOND_EXAMPLE_BEST_HIT_ASSIGNMENT_METHOD, "DIAMOND blastx best hit(s) of all reads from each OTU, but report the best hit as a sequence ID instead of a taxonomy."],
                [pipe.ANNOY_ASSIGNMENT_METHOD, "Search for the most similar window sequences <= 3bp different using ANNOY, otherwise no taxonomy is assigned. Requires a non-standard metapackage."],
                [pipe.PPLACER_ASSIGNMENT_METHOD, "Use pplacer to assign taxonomy of each read in each OTU. Requires a non-standard metapackage."]
            ]),
        default=current_default)

    argument_group.add_argument('--output-extras', action='store_true',
        help='give extra output for each sequence identified (e.g. the read(s) each OTU was generated from) in the output OTU table [default: not set]',
        default=False)

def add_less_common_pipe_arguments(argument_group, extra_args=False):
    argument_group.add_argument('--archive-otu-table', metavar='filename', help='output OTU table in archive format for making DBs etc. [default: unused]')
    argument_group.add_argument('--metapackage', help='Set of SingleM packages to use [default: use the default set]')
    argument_group.add_argument('--sra-files',
            nargs='+',
            metavar='sra_file',
            help='"sra" format files (usually from NCBI SRA) to be searched')
    argument_group.add_argument('--read-chunk-size',
            type=int,
            metavar='num_reads',
            help='Size chunk to process at a time (in number of reads). Requires --sra-files.')
    argument_group.add_argument('--read-chunk-number',
            type=int,
            metavar='chunk_number',
            help='Process only this specific chunk number (1-based index). Requires --sra-files.')
    argument_group.add_argument('--output-jplace', metavar='filename', help='Output a jplace format file for each singlem package to a file starting with this string, each with one entry per OTU. Requires \'%s\' as the --assignment_method [default: unused]' % pipe.PPLACER_ASSIGNMENT_METHOD)
    argument_group.add_argument('--singlem-packages', nargs='+', help='SingleM packages to use [default: use the set from the default metapackage]')
    argument_group.add_argument('--assignment-singlem-db', '--assignment_singlem_db', help='Use this SingleM DB when assigning taxonomy [default: not set, use the default]')
    argument_group.add_argument('--diamond-taxonomy-assignment-performance-parameters',
                                help='Performance-type arguments to use when calling \'diamond blastx\' during the taxonomy assignment step. [default: use setting defined in metapackage when set, otherwise use \'%s\']' % SearchPipe.DEFAULT_DIAMOND_ASSIGN_TAXONOMY_PERFORMANCE_PARAMETERS,
                                default=None)
    argument_group.add_argument('--evalue',
                                help='HMMSEARCH e-value cutoff to use for sequence gathering [default: %s]' % SearchPipe.DEFAULT_HMMSEARCH_EVALUE, default=SearchPipe.DEFAULT_HMMSEARCH_EVALUE)
    argument_group.add_argument('--min-orf-length',
                                metavar='length',
                                help='When predicting ORFs require this many base pairs uninterrupted by a stop codon [default: %i for reads, %i for genomes]' % (SearchPipe.DEFAULT_MIN_ORF_LENGTH, SearchPipe.DEFAULT_GENOME_MIN_ORF_LENGTH),
                                type=int)
    argument_group.add_argument('--restrict-read-length',
                                metavar='length',
                                help='Only use this many base pairs at the start of each sequence searched [default: no restriction]',
                                type=int)
    argument_group.add_argument('--translation-table',
                                metavar='number',
                                type=int,
                                help='Codon table for translation. By default, translation table 4 is used, which is the same as translation table 11 (the usual bacterial/archaeal one), except that the TGA codon is translated as tryptophan, not as a stop codon. Using table 4 means that the minority of organisms which use table 4 are not biased against, without a significant effect on the majority of bacteria and archaea that use table 11. See http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes for details on specific tables. [default: %i]' % SearchPipe.DEFAULT_TRANSLATION_TABLE,
                                default=SearchPipe.DEFAULT_TRANSLATION_TABLE)
    argument_group.add_argument('--filter-minimum-protein',
                                metavar='length',
                                help='Ignore reads aligning in less than this many positions to each protein HMM when using --no-diamond-prefilter [default: %i]' % SearchPipe.DEFAULT_FILTER_MINIMUM_PROTEIN,
                                type=int, default=SearchPipe.DEFAULT_FILTER_MINIMUM_PROTEIN)
    argument_group.add_argument('--max-species-divergence', metavar='INT',
                                help='Maximum number of different bases acids to allow between a sequence and the best hit in the database so that it is assigned to the species level. [default: %i]' % SearchPipe.DEFAULT_MAX_SPECIES_DIVERGENCE,
                                type=int, default=SearchPipe.DEFAULT_MAX_SPECIES_DIVERGENCE)
    argument_group.add_argument('--exclude-off-target-hits', action='store_true', help="Exclude hits that are not in the target domain of each SingleM package")
    argument_group.add_argument('--min-taxon-coverage',
                                metavar='FLOAT',
                                help='Minimum coverage to report in a taxonomic profile. [default: {} for reads, {} for genomes]'.format(CONDENSE_DEFAULT_MIN_TAXON_COVERAGE, CONDENSE_DEFAULT_GENOME_MIN_TAXON_COVERAGE),
                                type=float)
    
    if extra_args:
        argument_group.add_argument('--working-directory', metavar='directory', help='use intermediate working directory at a specified location, and do not delete it upon completion [default: not set, use a temporary directory]')
        argument_group.add_argument('--working-directory-dev-shm', default=False, action='store_true', help='use an intermediate results temporary working directory in /dev/shm rather than the default [default: the usual temporary working directory, currently {}]'.format(
            tempfile.gettempdir()
        ))
        argument_group.add_argument('--force', action='store_true', help='overwrite working directory if required [default: not set]')
        argument_group.add_argument('--filter-minimum-nucleotide',
                                    metavar='length',
                                    help='Ignore reads aligning in less than this many positions to each nucleotide HMM [default: %i]' % SearchPipe.DEFAULT_FILTER_MINIMUM_NUCLEOTIDE,
                                    type=int, default=SearchPipe.DEFAULT_FILTER_MINIMUM_NUCLEOTIDE)
        argument_group.add_argument('--include-inserts', action='store_true',
                                    help='print the entirety of the sequences in the OTU table, not just the aligned nucleotides [default: not set]', default=False)
        argument_group.add_argument('--known-otu-tables', nargs='+',
                                    help='OTU tables previously generated with trusted taxonomies for each sequence [default: unused]')
        argument_group.add_argument('--no-assign-taxonomy', action='store_true',
                                    help='Do not assign any taxonomy except for those already known [default: not set]',
                                    default=False)
        argument_group.add_argument('--known-sequence-taxonomy', metavar='FILE',
                                    help='A 2-column "sequence<tab>taxonomy" file specifying some sequences that have known taxonomy [default: unused]')
        argument_group.add_argument('--no-diamond-prefilter', action='store_true',
                                    help='Do not parse sequence data through DIAMOND blastx using a database constructed from the set of singlem packages. Should be used with --hmmsearch-package-assignment. NOTE: ignored for nucleotide packages [default: protein packages: use the prefilter, nucleotide packages: do not use the prefilter]',
                                    default=False)
        argument_group.add_argument('--diamond-prefilter-performance-parameters',
                                    help='Performance-type arguments to use when calling \'diamond blastx\' during the prefiltering. By default, SingleM should run in <4GB of RAM except in very large (>100Gbp) metagenomes. [default: use setting defined in metapackage when set, otherwise use \'%s\']' % SearchPipe.DEFAULT_PREFILTER_PERFORMANCE_PARAMETERS,
                                    default=None)
        argument_group.add_argument('--hmmsearch-package-assignment', '--hmmsearch_package_assignment', action='store_true',
                                    help='Assign each sequence to a SingleM package using HMMSEARCH, and a sequence may then be assigned to multiple packages. [default: not set]',
                                    default=False)
        argument_group.add_argument('--diamond-prefilter-db',
                                    help='Use this DB when running DIAMOND prefilter [default: use the one in the metapackage, or generate one from the SingleM packages]')
        argument_group.add_argument('--assignment-threads',type=int,
                                    help='Use this many processes in parallel while assigning taxonomy [default: %i]' % SearchPipe.DEFAULT_ASSIGNMENT_THREADS,
                                    default=SearchPipe.DEFAULT_ASSIGNMENT_THREADS)
        argument_group.add_argument('--sleep-after-mkfifo', type=int,
                                    help='Sleep for this many seconds after running os.mkfifo [default: None]')

def validate_pipe_args(args, subparser='pipe'):
    if not args.otu_table and not args.archive_otu_table and not args.taxonomic_profile and not args.taxonomic_profile_krona:
        raise Exception("At least one of --output-taxonomic-profile, --output-taxonomic-profile-krona, --otu-table, or --archive-otu-table must be specified")
    if args.output_jplace and args.assignment_method != pipe.PPLACER_ASSIGNMENT_METHOD:
        raise Exception("If --output-jplace is specified, then --assignment-method must be set to %s" % pipe.PPLACER_ASSIGNMENT_METHOD)
    if args.metapackage and args.singlem_packages:
        raise Exception("Can only specify a metapackage or a singlem package set, not both")
    if args.output_extras and not args.otu_table:
        raise Exception("Can't use --output-extras without --otu-table")
    if subparser == 'pipe':
        if args.include_inserts and not args.otu_table and not args.archive_otu_table:
            raise Exception("Can't use --include-inserts without --otu-table or --archive-otu-table")
        if args.metapackage and args.diamond_prefilter_db:
            raise Exception("Can't use a metapackage with --diamond-prefilter-db")
        if args.output_jplace and args.known_otu_tables:
            raise Exception("Currently --output-jplace and --known-otu-tables are incompatible")
        if args.output_jplace and args.no_assign_taxonomy:
            raise Exception("Currently --output-jplace and --no-assign-taxonomy are incompatible")
        if args.known_sequence_taxonomy and not args.no_assign_taxonomy:
            raise Exception(
                "Currently --known-sequence-taxonomy requires --no-assign-taxonomy to be set also")
        if args.reverse and args.output_jplace:
            raise Exception("Currently --jplace-output cannot be used with --reverse")
        if args.working_directory and args.working_directory_dev_shm:
            raise Exception("Cannot specify both --working-directory and --working-directory-dev-shm")
        if args.sra_files and args.no_diamond_prefilter:
            raise Exception("SRA input data requires a DIAMOND prefilter step, currently")
        if args.no_assign_taxonomy and (args.taxonomic_profile or args.taxonomic_profile_krona):
            raise Exception("Can't use --no-assign-taxonomy with --output-taxonomic-profile or --output-taxonomic-profile-krona")
        if args.read_chunk_size and not args.sra_files:
            raise Exception("Can't use --read-chunk-size without --sra-files")
        if args.read_chunk_number and not args.sra_files:
            raise Exception("Can't use --read-chunk-number without --sra-files")
        if bool(args.read_chunk_size) != bool(args.read_chunk_number):
            raise Exception("Either none or both of --read-chunk-size and --read-chunk-number should be set")
        if args.read_chunk_size and len(args.sra_files) > 1:
            raise Exception("Can't use --read-chunk-size with more than one --sra-file")

def add_condense_arguments(parser):
    input_condense_arguments = parser.add_argument_group("Input arguments (1+ required)")
    input_condense_arguments.add_argument('--input-archive-otu-tables', '--input-archive-otu-table', nargs='+', help="Condense from these archive tables")
    input_condense_arguments.add_argument('--input-archive-otu-table-list',
        help="Condense from the archive tables newline separated in this file")
    input_condense_arguments.add_argument('--input-gzip-archive-otu-table-list',
        help="Condense from the gzip'd archive tables newline separated in this file")

    output_condense_arguments = parser.add_argument_group("Output arguments (1+ required)")
    output_condense_arguments.add_argument('-p', '--taxonomic-profile', metavar='filename', help="output OTU table")
    output_condense_arguments.add_argument('--taxonomic-profile-krona', metavar='filename', help='name of krona file to generate.')
    output_condense_arguments.add_argument('--output-after-em-otu-table', metavar='filename', help="output OTU table after expectation maximisation has been applied. Note that this table usually contains multiple rows with the same window sequence.")

    optional_condense_arguments = parser.add_argument_group("Other options")
    optional_condense_arguments.add_argument('--metapackage', help='Set of SingleM packages to use [default: use the default set]')
    current_default = CONDENSE_DEFAULT_MIN_TAXON_COVERAGE
    optional_condense_arguments.add_argument('--min-taxon-coverage',metavar='FRACTION',
        help='Set taxons with less coverage to coverage=0. [default: {}]'.format(current_default), default=current_default, type=float)
    current_default = CONDENSE_DEFAULT_TRIM_PERCENT
    optional_condense_arguments.add_argument('--trim-percent', type=float, default=current_default, help="percentage of markers to be trimmed for each taxonomy [default: {}]".format(current_default))

def generate_streaming_otu_table_from_args(args,
    input_prefix=False, query_prefix=False, archive_only=False, min_archive_otu_table_version=None):

    if archive_only:
        otu_tables = False
        otu_tables_list = False
    if input_prefix:
        if not archive_only:
            otu_tables = args.input_otu_tables
            otu_tables_list = args.input_otu_tables_list
        archive_otu_tables = args.input_archive_otu_tables
        archive_otu_table_list = args.input_archive_otu_table_list
        gzip_archive_otu_table_list = args.input_gzip_archive_otu_table_list
    elif query_prefix:
        otu_tables = args.query_otu_table
        otu_tables_list = args.query_otu_tables_list
        archive_otu_tables = args.query_archive_otu_tables
        archive_otu_table_list = args.query_archive_otu_table_list
        gzip_archive_otu_table_list = args.query_gzip_archive_otu_table_list
    else:
        if not archive_only:
            otu_tables = args.otu_tables
            otu_tables_list = args.otu_tables_list
        archive_otu_tables = args.archive_otu_tables
        archive_otu_table_list = args.archive_otu_table_list
        gzip_archive_otu_table_list = args.gzip_archive_otu_table_list

    if archive_only:
        if not archive_otu_tables and not archive_otu_table_list and not gzip_archive_otu_table_list:
            raise Exception("{} requires input archive OTU tables".format(args.subparser_name))
    else:
        if not otu_tables and not otu_tables_list and not archive_otu_tables and \
            not archive_otu_table_list and not gzip_archive_otu_table_list:
            raise Exception("{} requires input OTU tables or archive OTU tables".format(args.subparser_name))

    from singlem.otu_table_collection import StreamingOtuTableCollection

    otus = StreamingOtuTableCollection()
    if min_archive_otu_table_version:
        otus.min_archive_otu_table_version = min_archive_otu_table_version
    if otu_tables:
        for o in otu_tables:
            otus.add_otu_table_file(o)
    if otu_tables_list:
        with open(otu_tables_list) as f:
            for o in f:
                otus.add_otu_table_file(o.strip())
    if archive_otu_tables:
        for o in archive_otu_tables:
            otus.add_archive_otu_table_file(o.strip())
    if archive_otu_table_list:
        with open(archive_otu_table_list) as f:
            for o in f.readlines():
                otus.add_archive_otu_table_file(o)
    if gzip_archive_otu_table_list:
        with open(gzip_archive_otu_table_list) as f:
            for arc in f.readlines():
                otus.add_gzip_archive_otu_table_file(arc.strip())
    return otus

def get_min_orf_length(args, subparser='pipe'):
    if args.min_orf_length:
        return args.min_orf_length
    elif subparser == 'pipe' and args.genome_fasta_files:
        return SearchPipe.DEFAULT_GENOME_MIN_ORF_LENGTH
    elif subparser in ('pipe', 'renew'):
        return SearchPipe.DEFAULT_MIN_ORF_LENGTH
    else:
        raise Exception("Programming error")

def get_min_taxon_coverage(args, subparser='pipe'):
    if args.min_taxon_coverage:
        return args.min_taxon_coverage
    elif subparser == 'pipe' and args.genome_fasta_files:
        return CONDENSE_DEFAULT_GENOME_MIN_TAXON_COVERAGE
    else:
        return CONDENSE_DEFAULT_MIN_TAXON_COVERAGE

def parse_genome_fasta_files(args):
    genomes = []
    if getattr(args, 'genome_fasta_files', None):
        genomes.extend(args.genome_fasta_files)
    if getattr(args, 'genome_fasta_directory', None):
        extension = getattr(args, 'genome_fasta_extension', 'fna')
        for fn in sorted(os.listdir(args.genome_fasta_directory)):
            if fn.endswith('.' + extension):
                genomes.append(os.path.join(args.genome_fasta_directory, fn))
    if getattr(args, 'genome_fasta_list', None):
        with open(args.genome_fasta_list) as f:
            genomes.extend([line.strip() for line in f if line.strip()])
    args.genome_fasta_files = genomes if genomes else None
    if args.genome_fasta_files:
        args.forward = args.genome_fasta_files
    return args

def main():
    bird_argparser = BirdArgparser(
        program='SingleM',
        authors=[
            BEN_NAME_AND_CENTRE,
            "Samuel Aroney, "+CMR,
            "Raphael Eisenhofer, Centre for Evolutionary Hologenomics, University of Copenhagen, Denmark",
            "Rossen Zhao, "+CMR],
        version=singlem.__version__,
        raw_format=True,
        examples={'pipe': [
            Example(
                'Get a taxonomic profile from paired read input:',
                'singlem pipe -1 <fastq_or_fasta1> -2 <fastq_or_fasta2> -p <output.profile.tsv>'),
            Example(
                'Get a taxonomic profile Krona diagram from single read input (long or short read):',
                'singlem pipe -1 <fastq_or_fasta> --taxonomic-profile-krona <output.profile.html>'),
            Example(
                'Gather an OTU table (per marker sequence groupings) from paired reads:',
                'singlem pipe -1 <fastq_or_fasta1> -2 <fastq_or_fasta2> --otu-table <output.otu_table.tsv>'),
        ],
        'summarise': [
            Example(
                'Convert a taxonomic profile to a site-by-species table at the genus level:',
                'singlem summarise --input-taxonomic-profiles <profile1.tsv> <profile2.tsv> \\\n'
                    '    --output-species-by-site-relative-abundance <output.tsv> \\\n'
                    '    --output-species-by-site-level genus'),
            Example(
                'Create a Krona diagram from a taxonomic profile:',
                'singlem summarise --input-taxonomic-profiles <profile1.tsv> \\\n'
                    '    --output-taxonomic-profile-krona <output.html>'),
            Example(
                'Add extra coverage and relative abundance information to a taxonomic profile:',
                'singlem summarise --input-taxonomic-profiles <profile1.tsv> \\\n'
                    '    --output-taxonomic-profile-with-extras <output.tsv>'),
        ]}
    )

    data_description = 'Download reference metapackage data'
    data_parser = bird_argparser.new_subparser('data', data_description, parser_group='Tools')
    # TODO: Could make pipe invocation faster by moving DATA_ENVIRONMENT to a separate file
    data_parser.add_argument('--output-directory', help="Output directory [required unless {} is specified]".format(DATA_ENVIRONMENT_VARIABLE))
    data_parser.add_argument('--verify-only', help="Check that the data is up to date and each file has the correct checksum", action='store_true', default=False)

    pipe_description = 'Generate a taxonomic profile or OTU table from raw sequences'
    pipe_parser = bird_argparser.new_subparser('pipe', pipe_description, parser_group='Tools')

    common_pipe_arguments = pipe_parser.add_argument_group('Common options')
    add_common_pipe_arguments(common_pipe_arguments, extra_args=True)

    less_common_pipe_arguments = pipe_parser.add_argument_group('Less common options')
    add_less_common_pipe_arguments(less_common_pipe_arguments, extra_args=True)

    appraise_description = 'How much of the metagenome do the genomes or assembly represent?'
    appraise_parser = bird_argparser.new_subparser('appraise', appraise_description, parser_group='Tools')
    appraise_otu_table_options = appraise_parser.add_argument_group('Input OTU table options')
    appraise_otu_table_options.add_argument('--metagenome-otu-tables', nargs='+', help="output of 'pipe' run on metagenomes")
    appraise_otu_table_options.add_argument('--metagenome-archive-otu-tables', nargs='+', help="archive output of 'pipe' run on metagenomes")
    appraise_otu_table_options.add_argument('--genome-otu-tables', nargs='+', help="output of 'pipe' run on genomes")
    appraise_otu_table_options.add_argument('--genome-archive-otu-tables', nargs='+', help="archive output of 'pipe' run on genomes")
    appraise_otu_table_options.add_argument('--assembly-otu-tables', nargs='+', help="output of 'pipe' run on assembled sequence")
    appraise_otu_table_options.add_argument('--assembly-archive-otu-tables', nargs='+', help="archive output of 'pipe' run on assembled sequence")
    appraise_otu_table_options.add_argument('--metapackage', help='Metapackage used in the creation of the OTU tables')
    appraise_inexact_options = appraise_parser.add_argument_group('Inexact appraisal options')
    appraise_inexact_options.add_argument('--imperfect', action='store_true', help="use sequence searching to account for genomes that are similar to those found in the metagenome [default: False]", default=False)
    appraise_inexact_options.add_argument('--sequence-identity', type=float, help="sequence identity cutoff to use if --imperfect is specified [default: ~species level divergence i.e. %s]" % SPECIES_LEVEL_AVERAGE_IDENTITY, default=SPECIES_LEVEL_AVERAGE_IDENTITY)
    appraise_plot_group = appraise_parser.add_argument_group("Plotting-related options")
    appraise_plot_group.add_argument('--plot', help='Output plot SVG filename (marker chosen automatically unless --plot-marker is also specified)', default=None)
    appraise_plot_group.add_argument('--plot-marker', help='Marker gene to plot OTUs from', default=None)
    appraise_plot_group.add_argument('--plot-basename', help="Plot visualisation of appraisal results from all markers to this basename (one SVG per marker)", default=None)
    appraise_otu_table_group = appraise_parser.add_argument_group('Output summary OTU tables')
    appraise_otu_table_group.add_argument('--output-binned-otu-table', help="output OTU table of binned populations", default=None)
    appraise_otu_table_group.add_argument('--output-unbinned-otu-table', help="output OTU table of assembled but not binned populations", default=None)
    appraise_otu_table_group.add_argument('--output-assembled-otu-table', help="output OTU table of all assembled populations", default=None)
    appraise_otu_table_group.add_argument('--output-unaccounted-for-otu-table', help="Output OTU table of populations not accounted for", default=None)
    appraise_otu_table_group.add_argument('--output-found-in', action='store_true', help="Output sample name (genome or assembly) the hit was found in")
    appraise_otu_table_group.add_argument('--output-style', help="Style of output OTU tables", default=OTU_TABLE_OUTPUT_FORMAT,
                                          choices=[OTU_TABLE_OUTPUT_FORMAT, ARCHIVE_TABLE_OUTPUT_FORMAT])
    appraise_otu_table_group.add_argument('--stream-inputs', action='store_true', help="Stream input OTU tables, saving RAM. Only works with --output-otu-table and transformation options do not work [expert option].")
    default_appraise_threads = 1
    appraise_otu_table_group.add_argument('--threads', type=int, metavar='num_threads', help='Use this many threads when processing streaming inputs [default %i]' % default_appraise_threads, default=default_appraise_threads)

    seqs_description = 'Find the best window position for a SingleM package'
    seqs_parser = bird_argparser.new_subparser('seqs', seqs_description)

    seqs_parser.add_argument('--alignment', metavar='aligned_fasta', help="Protein sequences hmmaligned and converted to fasta format with seqmagick", required=True)
    seqs_parser.add_argument('--alignment-type', metavar='type', help="alignment is 'aa' or 'dna'", required=True)
    seqs_parser.add_argument('--window-size', metavar='INT',
        help='Number of nucleotides to use in continuous window [default: {}]'.format(DEFAULT_WINDOW_SIZE),
        default=DEFAULT_WINDOW_SIZE, type=int)
    seqs_parser.add_argument('--hmm', help="HMM file used to generate alignment, used here to rank windows according to their information content.")

    makedb_description = 'Create a searchable OTU sequence database from an OTU table'
    makedb_parser = bird_argparser.new_subparser('makedb', makedb_description)

    required_makedb_arguments = makedb_parser.add_argument_group('required arguments')

    required_makedb_arguments.add_argument('--otu-tables', '--otu-table', nargs='+', help="Make a db from these OTU tables")
    required_makedb_arguments.add_argument('--otu-tables-list', help="Make a db from the OTU table files newline separated in this file")
    required_makedb_arguments.add_argument('--archive-otu-tables', '--archive-otu-table', nargs='+', help="Make a db from these archive tables")
    required_makedb_arguments.add_argument('--archive-otu-table-list',
        help="Make a db from the archive tables newline separated in this file")
    required_makedb_arguments.add_argument('--gzip-archive-otu-table-list',
        help="Make a db from the gzip'd archive tables newline separated in this file")
    required_makedb_arguments.add_argument('--db', help="Name of database to create e.g. tundra.sdb", required=True)
    makedb_other_args = makedb_parser.add_argument_group('Other arguments')
    makedb_other_args.add_argument('--threads', help='Use this many threads where possible [default 1]')
    current_default = ['smafa-naive']
    makedb_other_args.add_argument('--sequence-database-methods',
        nargs='+',
        choices = ['smafa-naive','annoy','scann','nmslib','scann-naive','none'],
        default = current_default,
        help='Index sequences using these methods. Note that specifying "scann-naive" means "scann" databases will also be built [default {}]'.format(current_default))
    current_default = ['nucleotide']
    makedb_other_args.add_argument('--sequence-database-types', help='Index sequences using these types. [default: {}]'.format(current_default), nargs='+', default=['nucleotide'], choices=['nucleotide','protein'])
    makedb_other_args.add_argument('--pregenerated-otu-sqlite-db', help='[for internal usage] remake the indices using this input SQLite database')
    DEFAULT_ANNOY_NUCLEOTIDE_NTREES = 10
    makedb_other_args.add_argument('--num-annoy-nucleotide-trees', help='make annoy nucleotide sequence indices with this ntrees [default {}]'.format(DEFAULT_ANNOY_NUCLEOTIDE_NTREES),default=DEFAULT_ANNOY_NUCLEOTIDE_NTREES,type=int)
    DEFAULT_ANNOY_PROTEIN_NTREES = 10
    makedb_other_args.add_argument('--num-annoy-protein-trees', help='make annoy protein sequence indices with this ntrees [default {}]'.format(DEFAULT_ANNOY_PROTEIN_NTREES),default=DEFAULT_ANNOY_PROTEIN_NTREES,type=int)
    makedb_other_args.add_argument('--tmpdir', help='[for internal usage] use this directory internally for working')

    query_description = 'Find closely related sequences in a SingleM database.'
    query_parser = bird_argparser.new_subparser('query', query_description)
    query_db_args = query_parser.add_argument_group('Required arguments')
    query_db_args.add_argument('--db', help="Output from 'makedb' mode", required=True)
    query_otu_args = query_parser.add_argument_group('Database querying by OTU sequence')
    query_otu_args.add_argument('--query-otu-table', '--query-otu-tables', nargs='+', metavar='file', help="Query the database with all sequences in this OTU table")
    query_otu_args.add_argument('--query-otu-tables-list', help="Query the database with all sequences in OTU table files newline separated in this file")
    query_otu_args.add_argument('--query-archive-otu-tables', nargs='+', help="Query the database with all sequences in these archive tables")
    query_otu_args.add_argument('--query-archive-otu-table-list', help="Query the database with all sequences in archive tables newline separated in this file")
    query_otu_args.add_argument('--query-gzip-archive-otu-table-list', help="Query the database with all sequences in gzip'd archive tables newline separated in this file")
    current_default = 20
    query_otu_args.add_argument('--max-nearest-neighbours', help="How many nearest neighbours to report. Each neighbour is a distinct sequence from the DB. [default: {}]".format(current_default), type=int, default=current_default)
    query_otu_args.add_argument('--max-divergence', metavar='INT', help="Report sequences less than or equal to this divergence i.e. number of different bases/amino acids", type=int)
    current_default = 'smafa-naive'
    query_otu_args.add_argument('--search-method', help="Algorithm to perform search [default: {}]".format(current_default), default=current_default, choices=['smafa-naive','nmslib','annoy','scann','scann-naive'])
    current_default = 'nucleotide'
    query_otu_args.add_argument('--sequence-type', help="Which sequence types to compare (i.e. protein for blastp, nucleotide for blastn) [default: {}]".format(current_default), default=current_default,
        choices=['nucleotide','protein'])
    current_default = 100
    query_otu_args.add_argument('--max-search-nearest-neighbours', help="How many nearest neighbours to search for with approximate nearest neighbours. Of these hits, only --max-nearest-neighbours will actually be reported. Ignored for --search-method naive and scann-naive. [default: {}]".format(current_default), type=int, default=current_default)
    # query_otu_args.add_argument('--stream-output','--stream_output', help='Stream output. Results may not be sorted by divergence [default: do not]', action='store_true')
    current_default = 1
    query_otu_args.add_argument('--threads', help='Use this many threads where possible [default %i]' % current_default, default=current_default)
    query_otu_args.add_argument('--limit-per-sequence',type=int, help='How many entries (samples/genomes from DB with identical sequences) to report for each distinct, matched sequence (arbitrarily chosen) [default: No limit]')
    query_otu_args.add_argument('--preload-db', action='store_true', help='Cache all DB data in python-land instead of querying for it by SQL each time. This is faster particularly for querying many sequences, but uses more memory and has a larger start-up time for each marker gene.')
    query_other_args = query_parser.add_argument_group('Other database extraction methods')
    query_other_args.add_argument('--sample-names', metavar='name', help='Print all OTUs from these samples', nargs='+')
    query_other_args.add_argument('--sample-list', metavar='path', help='Print all OTUs from the samples listed in the file (newline-separated)')
    query_other_args.add_argument('--taxonomy', metavar='name', help='Print all OTUs assigned a taxonomy including this string e.g. \'Archaea\'')
    query_other_args.add_argument('--dump', action='store_true', help='Print all OTUs in the DB')
    # continue_on_missing_genes
    query_other_args.add_argument('--continue-on-missing-genes', action='store_true', help='Continue if a gene is missing from the DB. Only works with smafa/nuclotide search method.', default=False)
    current_default = None

    summarise_description = 'Summarise and transform taxonomic profiles and OTU tables.'
    summarise_parser = bird_argparser.new_subparser('summarise', summarise_description, parser_group='Tools')

    summarise_taxonomic_profile_input_args = summarise_parser.add_argument_group('Taxonomic profile input')
    summarise_taxonomic_profile_input_args.add_argument('--input-taxonomic-profiles', nargs='+', help='Input taxonomic profiles to be e.g. converted to krona HTML, or concatenated')

    summarise_taxonomic_profile_output_args = summarise_parser.add_argument_group('Taxonomic profile output')
    summarise_taxonomic_profile_output_args.add_argument('--output-taxonomic-profile', metavar='FILE', help="Output a single output file containing taxonomic profiles of all input taxonomic profile files. Requires --input-taxonomic-profiles")
    summarise_taxonomic_profile_output_args.add_argument('--output-taxonomic-profile-krona', metavar='FILE', help="Output taxonomic profile to this file in Krona format.")
    summarise_taxonomic_profile_output_args.add_argument('--output-species-by-site-relative-abundance', metavar='FILE', help="Output site by species relative abundance to this file")
    summarise_taxonomic_profile_output_args.add_argument('--output-species-by-site-level', help="Output site by species level to this file. Requires --output-species-by-site-relative-abundance.", choices=['species','genus','family','order','class','phylum','domain'], default='species')
    summarise_taxonomic_profile_output_args.add_argument('--output-species-by-site-relative-abundance-prefix', metavar='PATH_PREFIX', help="Output site by species relative abundance to this file prefix. One file will be written for each taxonomic level.")
    summarise_taxonomic_profile_output_args.add_argument('--output-filled-taxonomic-profile', metavar='FILE', help="Output a taxonomic profile where the coverage of each taxon includes the coverage of each of its descendent taxons e.g. the d__Bacteria entry includes the p__Patescibacteria entry.")
    summarise_taxonomic_profile_output_args.add_argument('--output-taxonomic-profile-with-extras', metavar='FILE', help="Output a taxonomic profile with extra information (coverage, 'filled' coverage, relative abundance, taxonomy level).")
    summarise_taxonomic_profile_output_args.add_argument('--num-decimal-places', metavar='INT', type=int, help="Number of decimal places to report in the coverage column of the --output-taxonomic-profile-with-extras [default: 2].")
    summarise_taxonomic_profile_output_args.add_argument('--output-taxonomic-level-coverage', metavar='FILE', help="Output summary of how much coverage has been assigned to each taxonomic level in a taxonomic profile to a TSV file.")
    
    summarise_otu_table_input_args = summarise_parser.add_argument_group('OTU table input')
    summarise_otu_table_input_args.add_argument('--input-otu-tables', '--input-otu-table', nargs='+', help="Summarise these tables")
    summarise_otu_table_input_args.add_argument('--input-otu-tables-list', help="Summarise the OTU table files newline separated in this file")
    summarise_otu_table_input_args.add_argument('--input-archive-otu-tables', '--input-archive-otu-table', nargs='+', help="Summarise these tables", default=[])
    summarise_otu_table_input_args.add_argument('--input-archive-otu-table-list',
        help="Summarise the archive tables newline separated in this file")
    summarise_otu_table_input_args.add_argument('--input-gzip-archive-otu-table-list',
        help="Summarise the list of newline-separated gzip-compressed archive OTU tables specified in this file")
    summarise_otu_table_input_args.add_argument('--stream-inputs', help='Stream input OTU tables, saving RAM. Only works with --output-otu-table and transformation options do not work [expert option].', action='store_true')

    summarise_transformation_args = summarise_parser.add_argument_group('OTU table transformation')
    summarise_transformation_args.add_argument('--cluster', action='store_true', help="Apply sequence clustering to the OTU table. Any dashes in OTU sequences will be replaced by N.")
    summarise_transformation_args.add_argument('--cluster-id', type=float, help="Sequence clustering identity cutoff if --cluster is used [default: {} i.e. {}%]".format(SPECIES_LEVEL_AVERAGE_IDENTITY, SPECIES_LEVEL_AVERAGE_IDENTITY*100), default=SPECIES_LEVEL_AVERAGE_IDENTITY)
    summarise_transformation_args.add_argument('--taxonomy', help="Restrict analysis to OTUs that have this taxonomy (exact taxonomy or more fully resolved)")
    summarise_transformation_args.add_argument('--rarefied-output-otu-table', help="Output rarefied output OTU table, where each gene and sample combination is rarefied")
    summarise_transformation_args.add_argument('--number-to-choose', type=int, help="Rarefy using this many sequences. Sample/gene combinations with an insufficient number of sequences are ignored with a warning [default: maximal number such that all samples have sufficient counts]")
    summarise_transformation_args.add_argument('--collapse-to-sample-name', help="Merge all OTUs into a single OTU table, using the given sample name. Requires archive OTU table input and output.")
    summarise_transformation_args.add_argument('--collapse-coupled', action='store_true', help="Merge forward and reverse read OTU tables into a unified table. Sample names of coupled reads must end in '1' and '2' respectively. Read names are ignored, so that if the forward and reverse from a pair contain the same OTU sequence, they will each count separately.")
    summarise_transformation_args.add_argument('--collapse-paired-with-unpaired-archive-otu-table', help="For archive OTU tables that have both paired and unpaired components, merge these into a single output archive OTU table")

    summarise_output_args = summarise_parser.add_argument_group('OTU table output')
    summarise_output_args.add_argument('--output-otu-table', help="Output combined OTU table to this file")
    summarise_output_args.add_argument('--output-archive-otu-table', help="Output combined OTU table to this file")
    summarise_output_args.add_argument('--output-translated-otu-table', help="Output combined OTU table to this file, with seqeunces translated into amino acids")
    summarise_output_args.add_argument('--output-extras', action='store_true', help="Output extra information in the standard output OTU table", default=False)
    summarise_output_args.add_argument('--krona', help="Name of krona file to generate. Note that this generates a krona file from the OTU table, not the taxonomic profile")
    summarise_output_args.add_argument('--wide-format-otu-table', help="Name of output species by site CSV file")
    summarise_output_args.add_argument('--strain-overview-table', help="Name of output strains table to generate")
    summarise_output_args.add_argument('--unifrac-by-otu', help="Output UniFrac format file where entries are OTU sequences")
    summarise_output_args.add_argument('--unifrac-by-taxonomy', help="Output UniFrac format file where entries are taxonomies (generally used for phylogeny-driven beta diversity when pipe was run with '--assignment_method diamond_example')")
    summarise_output_args.add_argument('--clustered-output-otu-table', help="Output an OTU table with extra information about the clusters. To simply cluster an OTU table, use --cluster with --output-otu-table instead.")
    summarise_output_args.add_argument('--exclude-off-target-hits', action='store_true', help="Exclude hits that are not in the target domain of each SingleM package")
    summarise_output_args.add_argument('--singlem-packages', nargs='+', help="Packages used in the creation of the OTU tables")
    summarise_output_args.add_argument('--metapackage', help='Metapackage used in the creation of the OTU tables')
    summarise_output_args.add_argument('--unaligned-sequences-dump-file',
        help="Output unaligned sequences from in put archive OTU table to this file. After each read name '~N' is added which corresponds to the order of the read in the archive OTU table, so that no two sequences have the same read name. N>1 can happen e.g. when the input file contains paired reads. ~0 does not necessarily correspond to the first read in the original input sequence set, but instead to the order in the input archive OTU table.")

    read_fraction_description = 'Estimate the fraction of reads from a metagenome that are assigned to Bacteria and Archaea compared to e.g. eukaryote or phage. Also estimate average genome size.'

    def add_prokaryotic_fraction_parser(name, description, deprecated=False):
        parser_group = 'exclude' if name == "microbial_fraction" else 'Tools'
        parser = bird_argparser.new_subparser(name, description, parser_group=parser_group)
        read_fraction_io_args = parser.add_argument_group('input')
        read_fraction_io_args.add_argument('-p', '--input-profile', help="Input taxonomic profile file [required]", required=True)
        read_fraction_sequence_input_group1 = parser.add_argument_group('Read information [1+ args required]')
        read_fraction_sequence_input_group = read_fraction_sequence_input_group1.add_mutually_exclusive_group(required=True)
        # Keep parity of these arguments with the 'pipe' command
        read_fraction_sequence_input_group.add_argument('-1','--forward','--reads','--sequences',
                                nargs='+',
                                metavar='sequence_file',
                                help='nucleotide read sequence(s) (forward or unpaired) to be searched. Can be FASTA or FASTQ format, GZIP-compressed or not. These must be the same ones that were used to generate the input profile.')
        read_fraction_sequence_input_group1.add_argument('-2', '--reverse',
                                nargs='+',
                                metavar='sequence_file',
                                help='reverse reads to be searched. Can be FASTA or FASTQ format, GZIP-compressed or not. These must be the same reads that were used to generate the input profile.')
        read_fraction_sequence_input_group.add_argument('--input-metagenome-sizes', help="TSV file with 'sample' and 'num_bases' as a header, where sample matches the input profile name, and num_reads is the total number (forward+reverse) of bases in the metagenome that was analysed with 'pipe'. These must be the same reads that were used to generate the input profile.")
        read_fraction_database_args = parser.add_argument_group('database')
        read_fraction_database_args.add_argument('--taxon-genome-lengths-file', help="TSV file with 'rank' and 'genome_size' as headers [default: Use genome lengths from the default metapackage]")
        read_fraction_database_args.add_argument('--metapackage', help="Metapackage containing genome lengths [default: Use genome lengths from the default metapackage]")
        read_fraction_uncommon_args = parser.add_argument_group('other options')
        read_fraction_uncommon_args.add_argument('--accept-missing-samples', action='store_true', help="If a sample is missing from the input-metagenome-sizes file, skip analysis of it without croaking.")
        read_fraction_uncommon_args.add_argument('--output-tsv', help="Output file [default: stdout]")
        read_fraction_uncommon_args.add_argument('--output-per-taxon-read-fractions', help="Output a fraction for each taxon to this TSV [default: Do not output anything]")
        return parser

    add_prokaryotic_fraction_parser('prokaryotic_fraction', read_fraction_description)
    add_prokaryotic_fraction_parser('microbial_fraction', read_fraction_description + ' [deprecated; use prokaryotic_fraction]')

    renew_description = 'Reannotate an OTU table with an updated taxonomy'
    renew_parser = bird_argparser.new_subparser('renew', renew_description, parser_group='Tools')
    renew_input_args = renew_parser.add_argument_group('input')
    renew_input_args.add_argument('--input-archive-otu-table', help="Renew this table", required=True)
    renew_input_args.add_argument('--ignore-missing-singlem-packages', help="Ignore OTUs which have been assigned to packages not in the metapackage being used for renewal [default: croak]", action='store_true')
    renew_common = renew_parser.add_argument_group("Common arguments in shared with 'pipe'")
    add_common_pipe_arguments(renew_common)
    renew_less_common = renew_parser.add_argument_group("Less common arguments shared with 'pipe'")
    add_less_common_pipe_arguments(renew_less_common)

    create_description = 'Create a SingleM package.'
    create_parser = bird_argparser.new_subparser('create', create_description)

    create_parser.add_argument('--input-graftm-package', metavar="PATH", help="Input GraftM package underlying the new SingleM package. The GraftM package is usually made with 'graftM create --no_tree --hmm <your.hmm>' where <your.hmm> is the one provided to 'singlem seqs'.", required=True)
    create_parser.add_argument('--input-taxonomy', metavar="PATH", help="Input taxonomy file in GreenGenes format (2 column tab separated, ID then taxonomy with taxonomy separated by ';' or '; '.", required=True)
    create_parser.add_argument('--output-singlem-package', metavar="PATH", help="Output package path", required=True)
    create_parser.add_argument('--hmm-position', metavar="INTEGER", help="Position in the GraftM alignment HMM where the SingleM window starts. To choose the best position, use 'singlem seqs'. Note that this position (both the one output by 'seqs' and the one specified here) is a 1-based index, but this positions stored within the SingleM package as a 0-based index.", required=True, type=int)
    create_parser.add_argument('--window-size', metavar="INTEGER",
        help="Length of NUCLEOTIDE residues in the window, counting only those that match the HMM [default: {}]".format(DEFAULT_WINDOW_SIZE),
        default=DEFAULT_WINDOW_SIZE, type=int)
    create_parser.add_argument('--target-domains', nargs='+', help="Input domains targeted by this package e.g. 'Archaea', 'Bacteria', 'Eukaryota' or 'Viruses'. Input with multiple domains must be space separated.", required=True)
    create_parser.add_argument('--gene-description', metavar="STRING", help="Input free form text description of this marker package, for use with 'singlem metapackage --describe'.", required=True, type=str)
    create_parser.add_argument('--force', action='store_true', help='Overwrite output path if it already exists [default: false]')

    get_tree_description = 'Extract path to Newick tree file in a SingleM package.'
    get_tree_parser = bird_argparser.new_subparser('get_tree', get_tree_description)

    required_get_tree_arguments = get_tree_parser.add_argument_group('required arguments')

    required_get_tree_arguments.add_argument('--singlem-packages', nargs='+', help='SingleM packages to use [default: use the default set]')

    regenerate_description = 'Update a SingleM package with new sequences and taxonomy (expert mode).'
    regenerate_parser = bird_argparser.new_subparser('regenerate', regenerate_description)
    current_default = singlem.CREATE_MIN_ALIGNED_PERCENT
    regenerate_parser.add_argument('--min-aligned-percent', metavar='percent', help="remove sequences from the alignment which do not cover this percentage of the HMM [default: {}]".format(current_default), type=int, default=current_default)
    regenerate_parser.add_argument('--window-position', help="change window position of output package [default: do not change]", type=int, default=False)
    regenerate_parser.add_argument('--sequence-prefix', help="add a prefix to sequence names", type=str, default="")
    regenerate_parser.add_argument('--candidate-decoy-sequences', '--euk-sequences', help='candidate amino acid sequences fasta file to search for decoys')
    regenerate_parser.add_argument('--candidate-decoy-taxonomy', '--euk-taxonomy', help='tab-separated sequence ID to taxonomy file of candidate decoy sequences')
    regenerate_parser.add_argument('--no-candidate-decoy-sequences', '--no-further-euks', help='Do not include any euk sequences beyond what is already in the current SingleM package', action='store_true')

    required_regenerate_arguments = regenerate_parser.add_argument_group('required arguments')
    required_regenerate_arguments.add_argument('--input-singlem-package', metavar="PATH", help="input package path", required=True)
    required_regenerate_arguments.add_argument('--output-singlem-package', metavar="PATH", help="output package path", required=True)
    required_regenerate_arguments.add_argument('--input-sequences', required=True, help='all on-target amino acid sequences fasta file for new package')
    required_regenerate_arguments.add_argument('--input-taxonomy', required=True, help='tab-separated sequence ID to taxonomy file of on-target sequence taxonomy')

    metapackage_description = 'Create or describe a metapackage (i.e. set of SingleM packages)'
    metapackage_parser = bird_argparser.new_subparser('metapackage', metapackage_description)
    metapackage_parser.add_argument('--metapackage', help='Path to write generated metapackage to')
    metapackage_parser.add_argument('--singlem-packages', nargs='+', help="Input packages")
    metapackage_parser.add_argument('--nucleotide-sdb', help="Nucleotide SingleM database for initial assignment pass")
    metapackage_parser.add_argument('--no-nucleotide-sdb', action='store_true', help="Skip nucleotide SingleM database")
    metapackage_parser.add_argument('--taxon-genome-lengths', help="TSV file of genome lengths for each taxon")
    metapackage_parser.add_argument('--no-taxon-genome-lengths', action='store_true', help="Skip taxon genome lengths")
    current_default = CUSTOM_TAXONOMY_DATABASE_NAME
    metapackage_parser.add_argument('--taxonomy-database-name', help='Name of the taxonomy database to use [default: %s]' % current_default, default=current_default)
    metapackage_parser.add_argument('--taxonomy-database-version', help='Version of the taxonomy database to use [default: unspecified]')
    metapackage_parser.add_argument('--diamond-prefilter-performance-parameters', help='Performance-type arguments to use when calling \'diamond blastx\' during the prefiltering. [default: \'%s\']' % SearchPipe.DEFAULT_PREFILTER_PERFORMANCE_PARAMETERS, default=SearchPipe.DEFAULT_PREFILTER_PERFORMANCE_PARAMETERS)
    metapackage_parser.add_argument('--diamond-taxonomy-assignment-performance-parameters', help='Performance-type arguments to use when calling \'diamond blastx\' during the taxonomy assignment. [default: \'%s\']' % SearchPipe.DEFAULT_DIAMOND_ASSIGN_TAXONOMY_PERFORMANCE_PARAMETERS, default=SearchPipe.DEFAULT_DIAMOND_ASSIGN_TAXONOMY_PERFORMANCE_PARAMETERS)

    metapackage_parser.add_argument('--describe', action='store_true', help='Describe a metapackage rather than create it')
    current_default = 1
    metapackage_parser.add_argument('--threads', type=int, metavar='num_threads', help='number of CPUS to use [default: %i]' % current_default, default=current_default)
    current_default = 0.6
    metapackage_parser.add_argument('--prefilter-clustering-threshold', type=float, metavar='fraction', help='ID for dereplication of prefilter DB [default: %s]' % current_default, default=current_default)
    metapackage_parser.add_argument('--prefilter-diamond-db', metavar='DMND', help='Dereplicated DIAMOND db for prefilter to use [default: dereplicate from input SingleM packages]')
    metapackage_parser.add_argument('--makeidx-sensitivity-params', metavar='PARAMS', help='DIAMOND sensitivity parameters to use when indexing the prefilter DIAMOND db. [default: None]', default=None)
    metapackage_parser.add_argument('--calculate-average-num-genes-per-species', action='store_true', help='Calculate the average number of genes per species in the metapackage. [default: False]', default=False)

    chainsaw_description = 'Remove tree information and trim unaligned sequences from a SingleM package (expert mode)'
    chainsaw_parser = bird_argparser.new_subparser('chainsaw', chainsaw_description)
    chainsaw_parser.add_argument('--keep-tree', action='store_true', help="Stop tree info from being removed", default=False)
    required_chainsaw_arguments = chainsaw_parser.add_argument_group("required arguments")
    required_chainsaw_arguments.add_argument('--input-singlem-package', required=True, help="Remove tree info and trim unaligned sequences from this package")
    required_chainsaw_arguments.add_argument('--output-singlem-package', required=True, help="Package to be created")
    required_chainsaw_arguments.add_argument('--sequence-prefix', default="", help="Rename the sequences by adding this at the front [default: '']")

    condense_description = 'Combine OTU tables across different markers into a single taxonomic profile. Note that while this mode can be run independently, it is often more straightforward to invoke its methodology by specifying -p / --taxonomic-profile when running pipe mode.'
    condense_parser = bird_argparser.new_subparser('condense', condense_description)
    add_condense_arguments(condense_parser)

    trim_package_hmms_description = 'Trim the width of HMMs to increase speed (expert mode)'
    trim_package_hmms_parser = bird_argparser.new_subparser('trim_package_hmms', trim_package_hmms_description)
    trim_package_hmms_parser.add_argument('--keep-tree', action='store_true', help="Stop tree info from being removed", default=False)
    required_trim_package_hmmsarguments = trim_package_hmms_parser.add_argument_group("required arguments")
    required_trim_package_hmmsarguments.add_argument('--input-singlem-package', required=True, help="Input package to trim HMMs from")
    required_trim_package_hmmsarguments.add_argument('--output-singlem-package', required=True, help="Package to be created")

    supplement_description = 'Create a new metapackage from a vanilla one plus new genomes'
    supplement_parser = bird_argparser.new_subparser('supplement', supplement_description, parser_group='Tools')
    supplement_parser.add_argument('--new-genome-fasta-files',
                                   help='FASTA files of new genomes',
                                   nargs='+')
    supplement_parser.add_argument('--new-genome-fasta-files-list',
                                   help='File containing FASTA file paths of new genomes',
                                   nargs='+')
    supplement_parser.add_argument('--input-metapackage', help='metapackage to build upon [default: Use default package]')
    supplement_parser.add_argument('--output-metapackage', help='output metapackage', required=True)
    supplement_parser.add_argument('--threads', help='parallelisation', type=int, default=1)

    # Taxonomy
    taxonomy_group = supplement_parser.add_argument_group('Taxonomy')
    taxonomy_mut_group = taxonomy_group.add_mutually_exclusive_group()
    taxonomy_mut_group.add_argument(
        '--new-fully-defined-taxonomies',
        help='newline separated file containing taxonomies of new genomes (path<TAB>taxonomy). Must be fully specified to species level. If not specified, the taxonomy will be inferred from the new genomes using GTDB-tk or read from --taxonomy-file [default: not set, run GTDBtk].'
    )
    taxonomy_mut_group.add_argument(
        '--taxonomy-file',
        help='A 2 column tab-separated file containing each genome\'s taxonomy as output by GTDBtk [default: not set, run GTDBtk]')
    taxonomy_mut_group.add_argument(
         '--gtdbtk-output-directory',
         help='use this GTDBtk result. Not used if --new-taxonomies is used [default: not set, run GTDBtk]')
    taxonomy_mut_group.add_argument('--pplacer-threads', help='for GTDBtk classify_wf', type=int)
    taxonomy_group.add_argument(
        '--output-taxonomies',
        help='TSV output file of taxonomies of new genomes, whether they are novel species or not.')
    taxonomy_group.add_argument('--skip-taxonomy-check', help='skip check which ensures that GTDBtk assigned taxonomies are concordant with the old metapackage\'s [default: do the check]', action='store_true')

    # Quality filtering
    quality_group = supplement_parser.add_argument_group('Quality filtering of new genomes')
    quality_group.add_argument('--checkm2-quality-file', help='CheckM2 quality file of new genomes')
    quality_group.add_argument('--no-quality-filter', help='skip quality filtering', action='store_true')
    current_default = 50
    quality_group.add_argument('--checkm2-min-completeness', help='minimum completeness for CheckM2 [default: {}]'.format(current_default), type=float, default=current_default)
    current_default = 10
    quality_group.add_argument('--checkm2-max-contamination', help='maximum contamination for CheckM2 [default: {}]'.format(current_default), type=float, default=current_default)

    # Dereplication
    dereplication_group = supplement_parser.add_argument_group('Dereplication')
    dereplication_mut_group = dereplication_group.add_mutually_exclusive_group(required=True)
    dereplication_mut_group.add_argument('--no-dereplication', help='Assume genome inputs are already dereplicated', action='store_true')
    dereplication_mut_group.add_argument('--dereplicate-with-galah', help='Run galah to dereplicate genomes at species level', action='store_true')

    # Rarer options
    supplement_rare_group = supplement_parser.add_argument_group('Less common options')
    current_default = '1e-20'
    supplement_rare_group.add_argument('--hmmsearch-evalue', help='evalue for hmmsearch run on proteins to gather markers [default: {}]'.format(current_default), default=current_default)
    supplement_rare_group.add_argument('--gene-definitions', help='Tab-separated file of genome_fasta<TAB>transcript_fasta<TAB>protein_fasta [default: undefined, call genes using Prodigal]')
    supplement_rare_group.add_argument('--working-directory', help='working directory [default: use a temporary directory]')
    supplement_rare_group.add_argument('--no-taxon-genome-lengths', help='Do not include taxon genome lengths in updated metapackage', action='store_true')
    supplement_rare_group.add_argument('--ignore-taxonomy-database-incompatibility', help='Do not halt when the old metapackage is not the default metapackage.', action='store_true')
    supplement_rare_group.add_argument('--new-taxonomy-database-name', help='Name of the taxonomy database to record in the created metapackage [default: %s]' % CUSTOM_TAXONOMY_DATABASE_NAME, default=CUSTOM_TAXONOMY_DATABASE_NAME)
    supplement_rare_group.add_argument('--new-taxonomy-database-version', help='Version of the taxonomy database to use [default: None]')

    args = bird_argparser.parse_the_args()
    parse_genome_fasta_files(args)

    if args.debug:
        loglevel = logging.DEBUG
        logging.getLogger('numba').setLevel(logging.INFO)
    elif args.quiet:
        logging.getLogger('nmslib').setLevel(logging.ERROR)
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
        logging.getLogger('nmslib').setLevel(logging.WARN)
        os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2'
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    logging.info("SingleM v{}".format(singlem.__version__))

    if args.subparser_name=='pipe':
        validate_pipe_args(args)
        singlem.pipe.SearchPipe().run(
            sequences = args.forward,
            reverse_read_files = args.reverse,
            input_sra_files = args.sra_files,
            read_chunk_size = args.read_chunk_size,
            read_chunk_number = args.read_chunk_number,
            otu_table = args.otu_table,
            archive_otu_table = args.archive_otu_table,
            sleep_after_mkfifo = args.sleep_after_mkfifo,
            threads = args.threads,
            known_otu_tables = args.known_otu_tables,
            assignment_method = args.assignment_method,
            assignment_threads = args.assignment_threads,
            output_jplace = args.output_jplace,
            evalue = args.evalue,
            min_orf_length = get_min_orf_length(args),
            translation_table = args.translation_table,
            restrict_read_length = args.restrict_read_length,
            filter_minimum_protein = args.filter_minimum_protein,
            filter_minimum_nucleotide = args.filter_minimum_nucleotide,
            output_extras = args.output_extras,
            include_inserts = args.include_inserts,
            working_directory = args.working_directory,
            working_directory_dev_shm = args.working_directory_dev_shm,
            force = args.force,
            metapackage_path = args.metapackage,
            singlem_packages = args.singlem_packages,
            assign_taxonomy = not args.no_assign_taxonomy,
            known_sequence_taxonomy = args.known_sequence_taxonomy,
            diamond_prefilter = not args.no_diamond_prefilter,
            diamond_prefilter_performance_parameters = args.diamond_prefilter_performance_parameters,
            diamond_package_assignment = not args.hmmsearch_package_assignment,
            diamond_prefilter_db = args.diamond_prefilter_db,
            diamond_taxonomy_assignment_performance_parameters = args.diamond_taxonomy_assignment_performance_parameters,
            assignment_singlem_db = args.assignment_singlem_db,
            output_taxonomic_profile = args.taxonomic_profile,
            output_taxonomic_profile_krona = args.taxonomic_profile_krona,
            viral_profile_output = False,
            exclude_off_target_hits = args.exclude_off_target_hits,
            min_taxon_coverage = get_min_taxon_coverage(args),
            max_species_divergence = args.max_species_divergence
        )

    elif args.subparser_name=='renew':
        from singlem.renew import Renew
        validate_pipe_args(args, subparser='renew')
        Renew().renew(
            input_archive_otu_table=args.input_archive_otu_table,
            ignore_missing_singlem_packages=args.ignore_missing_singlem_packages,
            otu_table = args.otu_table,
            output_archive_otu_table = args.archive_otu_table,
            threads = args.threads,
            assignment_method = args.assignment_method,
            output_jplace = args.output_jplace,
            output_extras = args.output_extras,
            metapackage_path = args.metapackage,
            singlem_packages = args.singlem_packages,
            diamond_taxonomy_assignment_performance_parameters = args.diamond_taxonomy_assignment_performance_parameters,
            evalue = args.evalue,
            min_orf_length = get_min_orf_length(args, subparser='renew'),
            restrict_read_length = args.restrict_read_length,
            filter_minimum_protein = args.filter_minimum_protein,
            assignment_singlem_db = args.assignment_singlem_db,
            output_taxonomic_profile = args.taxonomic_profile,
            output_taxonomic_profile_krona = args.taxonomic_profile_krona,
            exclude_off_target_hits = args.exclude_off_target_hits,
            translation_table = args.translation_table,
            max_species_divergence = args.max_species_divergence
            )

    elif args.subparser_name == 'summarise':
        from singlem.summariser import Summariser
        from singlem.strain_summariser import StrainSummariser
        from singlem.otu_table_collection import OtuTableCollection
        from singlem.clusterer import Clusterer
        from singlem.metapackage import Metapackage
        from singlem.singlem_package import SingleMPackage

        num_output_types = 0
        if args.strain_overview_table: num_output_types += 1
        if args.krona: num_output_types += 1
        if args.unifrac_by_otu: num_output_types += 1
        if args.unifrac_by_taxonomy: num_output_types += 1
        if args.output_otu_table: num_output_types += 1
        if args.output_archive_otu_table: num_output_types += 1
        if args.output_translated_otu_table: num_output_types += 1
        if args.clustered_output_otu_table: num_output_types += 1
        if args.rarefied_output_otu_table: num_output_types += 1
        if args.wide_format_otu_table: num_output_types += 1
        if args.collapse_paired_with_unpaired_archive_otu_table: num_output_types += 1
        if args.unaligned_sequences_dump_file: num_output_types += 1
        if args.output_taxonomic_profile: num_output_types += 1
        if args.output_taxonomic_profile_krona: num_output_types += 1
        if args.output_species_by_site_relative_abundance: num_output_types += 1
        if args.output_species_by_site_relative_abundance_prefix: num_output_types += 1
        if args.output_taxonomic_level_coverage: num_output_types += 1
        if args.output_filled_taxonomic_profile: num_output_types += 1
        if args.output_taxonomic_profile_with_extras: num_output_types += 1
        if num_output_types != 1:
            raise Exception("Exactly 1 output type must be specified, sorry, %i were provided" % num_output_types)
        if not args.input_otu_tables and \
            not args.input_otu_tables_list and \
            not args.input_archive_otu_tables and \
            not args.input_archive_otu_table_list and \
            not args.input_gzip_archive_otu_table_list and \
            not args.input_taxonomic_profiles:
            raise Exception("Summary requires input OTU tables, archive OTU tables, or taxonomic profiles")
        if args.exclude_off_target_hits:
            if args.singlem_packages and args.metapackage:
                raise Exception("Only one of --metapackage or --singlem-packages can be defined")
        if args.collapse_to_sample_name:
            if args.input_otu_tables:
                raise Exception("--collapse-to-sample-name currently only works with archive tables")
            elif not len(args.input_archive_otu_tables) >= 1 and not args.input_archive_otu_table_list and not args.input_gzip_archive_otu_table_list:
                raise Exception("--collapse-to-sample-name currently only works with archive tables as input")
        if args.collapse_paired_with_unpaired_archive_otu_table:
            if args.input_otu_tables:
                raise Exception("--collapse-paired-with-unpaired-archive-otu-table currently only works with archive tables")
            elif not len(args.input_archive_otu_tables) == 2:
                raise Exception("--collapse-paired-with-unpaired-archive-otu-table requires exactly two archive tables")
        if args.krona and args.input_taxonomic_profiles:
            raise Exception("--krona generates a krona file from an OTU table. Use --output-taxonomic-profile-krona when providing taxonomic profiles")
        if args.output_taxonomic_profile_krona and not args.input_taxonomic_profiles:
            raise Exception("--output-taxonomic-profile-krona requires --input-taxonomic-profiles to be defined")
        if args.output_archive_otu_table and not args.collapse_to_sample_name:
            raise Exception("--output-archive-otu-table requires --collapse-to-sample-name to be defined. Some other transfomations have output archive tables, but they are specified separately.")
        if args.output_species_by_site_relative_abundance and not args.input_taxonomic_profiles:
            raise Exception("--output-species-by-site-relative-abundance requires --input-taxonomic-profiles to be defined")
        if args.output_species_by_site_relative_abundance_prefix and not args.input_taxonomic_profiles:
            raise Exception("--output-species-by-site-relative-abundance-prefix requires --input-taxonomic-profiles to be defined")
        if args.output_taxonomic_level_coverage and not args.input_taxonomic_profiles:
            raise Exception("--output-taxonomic-level-coverage requires --input-taxonomic-profiles to be defined")
        if args.output_filled_taxonomic_profile and not args.input_taxonomic_profiles:
            raise Exception("--output-filled-taxonomic-profile requires --input-taxonomic-profiles to be defined")
        if args.output_taxonomic_profile_with_extras and not args.input_taxonomic_profiles:
            raise Exception("--output-taxonomic-profile-with-extras requires --input-taxonomic-profiles to be defined")
        if args.num_decimal_places and not args.output_taxonomic_profile_with_extras:
            raise Exception("--num-decimal-places currently requires --output-taxonomic-profile-with-extras to be defined")

        if args.stream_inputs or args.unaligned_sequences_dump_file:
            from singlem.otu_table_collection import StreamingOtuTableCollection
            if not args.output_otu_table and not args.unaligned_sequences_dump_file:
                raise Exception("--stream-inputs requires --output-otu-table or --unaligned-sequences-dump-file to be defined")
            if args.taxonomy:
                raise Exception("--stream-inputs does not currently support --taxonomy")
            require_archive_input = args.unaligned_sequences_dump_file is not None
            min_archive_otu_table_version = 2 if args.unaligned_sequences_dump_file else None
            otus = generate_streaming_otu_table_from_args(
                args,
                input_prefix=True,
                archive_only=require_archive_input,
                min_archive_otu_table_version=min_archive_otu_table_version)
        else:
            otus = OtuTableCollection()
            if args.input_otu_tables:
                for o in args.input_otu_tables:
                    otus.add_otu_table(open(o))
            if args.input_archive_otu_tables and not args.collapse_paired_with_unpaired_archive_otu_table and not args.collapse_to_sample_name:
                for o in args.input_archive_otu_tables:
                    otus.add_archive_otu_table(open(o))
            if args.input_otu_tables_list:
                with open(args.input_otu_tables_list) as f:
                    for o in f:
                        otus.add_otu_table(open(o.strip()))
            if args.input_archive_otu_table_list:
                with open(args.input_archive_otu_table_list) as f:
                    for o in f:
                        otus.add_archive_otu_table(open(o.strip()))
            if args.input_gzip_archive_otu_table_list:
                with open(args.input_gzip_archive_otu_table_list) as f:
                    for arc in f.readlines():
                        try:
                            otus.add_archive_otu_table(gzip.open(arc.strip()))
                        except json.decoder.JSONDecodeError:
                            logging.warning("Failed to parse JSON from archive OTU table {}, skipping".format(arc))
            otus.set_target_taxonomy_by_string(args.taxonomy)

        if args.collapse_coupled:
            if not args.output_otu_table:
                raise Exception("Collapsing is currently only implemented for regular OTU table outputs.")
            if args.stream_inputs:
                raise Exception("Streaming inputs is not currently known to work with collapse_coupled.")
            logging.info("Collapsing forward and reverse read OTU tables")
            o2 = OtuTableCollection()
            o2.otu_table_objects.append(otus.collapse_coupled())
            otus = o2

        if args.exclude_off_target_hits:
            if args.stream_inputs:
                raise Exception("Streaming inputs is not currently known to work with exclude_off_target_hits.")
            logging.info("Removing hits that are not assigned to their target domains")
            o2 = OtuTableCollection()
            if args.metapackage:
                pkgs = Metapackage.acquire(args.metapackage).singlem_packages
            elif args.singlem_packages:
                pkgs = []
                for p in args.singlem_packages:
                    pkg = SingleMPackage.acquire(p)
                    if pkg.version < 3:
                        raise Exception("Version 3 SingleM packages are required for --exclude_off_target_hits")
                    pkgs.append(pkg)
            else:
                mpkg = Metapackage.acquire_default()
                pkgs = mpkg.singlem_packages
            # Do not lose the extra columns if they are provided as input by
            # using return_archive_table=True. Currently an issue here that it
            # doesn't error when extra info is not provided as input.
            o2.otu_table_objects.append(otus.exclude_off_target_hits(pkgs, return_archive_table=True))
            otus = o2

        if args.cluster:
            if args.stream_inputs:
                raise Exception("Streaming inputs is not currently known to work with cluster.")
            logging.info("Clustering OTUs with clustering identity %f.." % args.cluster_id)
            o2 = OtuTableCollection()
            o2.otu_table_objects = [list(Clusterer().each_cluster(otus, args.cluster_id))]
            otus = o2
            logging.info("Finished clustering")

        if args.krona:
            Summariser.write_otu_table_krona(
                table_collection = otus,
                krona_output = args.krona)
        elif args.unifrac_by_otu:
            Summariser.write_unifrac_by_otu_format_file(
                table_collection = otus,
                unifrac_output_prefix = args.unifrac_by_otu)
        elif args.unifrac_by_taxonomy:
            Summariser.write_unifrac_by_taxonomy_format_file(
                table_collection = otus,
                unifrac_output_prefix = args.unifrac_by_taxonomy)
        elif args.output_otu_table:
            with open(args.output_otu_table, 'w') as f:
                Summariser.write_otu_table(
                    table_collection = otus,
                    output_table_io = f,
                    output_extras = args.output_extras)
        elif args.wide_format_otu_table:
            with open(args.wide_format_otu_table, 'w') as f:
                Summariser.write_wide_format_otu_table(
                    table_collection = otus,
                    output_table_io = f)
        elif args.strain_overview_table:
            with open(args.strain_overview_table, 'w') as f:
                StrainSummariser().summarise_strains(
                    table_collection = otus,
                    output_table_io = f)
        elif args.clustered_output_otu_table:
            if not args.cluster:
                raise Exception("If --clustered-output-otu-table is set, then clustering (--cluster) must be applied")
            with open(args.clustered_output_otu_table, 'w') as f:
                Summariser.write_clustered_otu_table(
                    table_collection = otus,
                    output_table_io = f)
        elif args.rarefied_output_otu_table:
            with open(args.rarefied_output_otu_table, 'w') as f:
                Summariser.write_rarefied_otu_table(
                    table_collection = otus,
                    output_table_io = f,
                    number_to_choose = args.number_to_choose)
        elif args.output_translated_otu_table:
            with open(args.output_translated_otu_table, 'w') as f:
                Summariser.write_translated_otu_table(
                    table_collection = otus,
                    output_table_io = f)
        elif args.collapse_to_sample_name:
            if not args.output_archive_otu_table:
                raise Exception("If --collapse-to-sample-name is set, then --output-archive-otu-table must be set")
            with open(args.output_archive_otu_table, 'w') as f:
                Summariser.write_collapsed_paired_with_unpaired_otu_table(
                    archive_otu_tables = args.input_archive_otu_tables,
                    archive_otu_table_list = args.input_archive_otu_table_list,
                    gzip_archive_otu_table_list = args.input_gzip_archive_otu_table_list,
                    output_table_io = f,
                    set_sample_name = args.collapse_to_sample_name)
        elif args.collapse_paired_with_unpaired_archive_otu_table:
            with open(args.collapse_paired_with_unpaired_archive_otu_table,'w') as output_io:
                Summariser.write_collapsed_paired_with_unpaired_otu_table(
                    archive_otu_tables = args.input_archive_otu_tables,
                    archive_otu_table_list = args.input_archive_otu_table_list,
                    gzip_archive_otu_table_list = args.input_gzip_archive_otu_table_list,
                    output_table_io = output_io)
        elif args.unaligned_sequences_dump_file:
            with open(args.unaligned_sequences_dump_file, 'w') as f:
                Summariser.dump_raw_sequences_from_archive_otu_table(
                    table_collection = otus,
                    output_table_io = f)
        elif args.input_taxonomic_profiles:
            from singlem.condense import CondensedCommunityProfile, CondensedCommunityProfileKronaWriter
            if args.output_taxonomic_profile_krona:
                profiles = []
                logging.info("Reading taxonomic profiles ..")
                for p in args.input_taxonomic_profiles:
                    with open(p) as f:
                        for sample_profile in CondensedCommunityProfile.each_sample_wise(f):
                            profiles.append(sample_profile)
                logging.info("Writing krona to %s" % args.output_taxonomic_profile_krona)
                CondensedCommunityProfileKronaWriter.write_krona(
                    profiles, args.output_taxonomic_profile_krona
                )
            elif args.output_species_by_site_relative_abundance:
                Summariser.write_species_by_site_table(
                    input_taxonomic_profiles = args.input_taxonomic_profiles,
                    output_species_by_site_relative_abundance_table = args.output_species_by_site_relative_abundance,
                    level = args.output_species_by_site_level)
            elif args.output_species_by_site_relative_abundance_prefix:
                Summariser.write_species_by_site_table(
                    input_taxonomic_profiles = args.input_taxonomic_profiles,
                    output_species_by_site_relative_abundance_table = args.output_species_by_site_relative_abundance_prefix,
                    level=None)
            elif args.output_taxonomic_level_coverage:
                Summariser.write_taxonomic_level_coverage_table(
                    input_taxonomic_profiles = args.input_taxonomic_profiles,
                    output_taxonomic_level_coverage_table = args.output_taxonomic_level_coverage)
            elif args.output_taxonomic_profile:
                with open(args.output_taxonomic_profile, 'w') as f:
                    Summariser.write_taxonomic_profile(
                        input_taxonomic_profiles = args.input_taxonomic_profiles,
                        output_taxonomic_profile_io = f)
            elif args.output_filled_taxonomic_profile:
                with open(args.output_filled_taxonomic_profile, 'w') as f:
                    Summariser.write_filled_taxonomic_profile(
                        input_taxonomic_profile_files = args.input_taxonomic_profiles,
                        output_filled_taxonomic_profile_io = f)
            elif args.output_taxonomic_profile_with_extras:
                with open(args.output_taxonomic_profile_with_extras, 'w') as f:
                    Summariser.write_taxonomic_profile_with_extras(
                        input_taxonomic_profile_files = args.input_taxonomic_profiles,
                        output_taxonomic_profile_extras_io = f,
                        num_decimal_places = args.num_decimal_places)
            else:
                raise Exception("Expected --output-taxonomic-profile-krona or --output-site-by-species-relative-abundance or --output-taxonomic-level-coverage to be defined, since --input-taxonomic-profiles was defined")

        else: raise Exception("Programming error")
        logging.info("Finished")

    elif args.subparser_name == 'create':
        from singlem.package_creator import PackageCreator
        PackageCreator().create(input_graftm_package = args.input_graftm_package,
                                input_taxonomy = args.input_taxonomy,
                                output_singlem_package = args.output_singlem_package,
                                hmm_position = args.hmm_position,
                                window_size = args.window_size,
                                target_domains = args.target_domains,
                                gene_description = args.gene_description,
                                force = args.force)
    elif args.subparser_name == 'appraise':
        from singlem.otu_table_collection import OtuTableCollection, StreamingOtuTableCollection
        from singlem.appraiser import Appraiser
        from singlem.metapackage import Metapackage

        if not args.metagenome_otu_tables and not args.metagenome_archive_otu_tables:
            raise Exception("Appraise requires metagenome OTU tables or archive tables")

        if args.output_style == ARCHIVE_TABLE_OUTPUT_FORMAT and not args.metagenome_archive_otu_tables:
            raise Exception("Appraise archive output format requires metagenome archive input")

        appraiser = Appraiser()

        if args.stream_inputs:
            logging.info("Preparing file IO for streaming inputs")
            metagenomes = StreamingOtuTableCollection()
            file_io = []
            if args.metagenome_otu_tables:
                for table in args.metagenome_otu_tables:
                    table_io = open(table)
                    metagenomes.add_otu_table(table_io)
                    file_io.append(open(table))
            if args.metagenome_archive_otu_tables:
                for table in args.metagenome_archive_otu_tables:
                    table_io = open(table)
                    metagenomes.add_archive_otu_table(table_io)
                    file_io.append(open(table))
        else:
            metagenomes = OtuTableCollection()
            if args.metagenome_otu_tables:
                for table in args.metagenome_otu_tables:
                    with open(table) as f:
                        metagenomes.add_otu_table(f)
            if args.metagenome_archive_otu_tables:
                for table in args.metagenome_archive_otu_tables:
                    with open(table) as f:
                        metagenomes.add_archive_otu_table(f)

        if args.metapackage:
            pkgs = Metapackage.acquire(args.metapackage).singlem_packages
        else:
            mpkg = Metapackage.acquire_default()
            pkgs = mpkg.singlem_packages

        if args.stream_inputs:
            logging.info("Removing non-target hits is not compatible with streaming inputs")
        else:
            logging.info("Removing hits that are not assigned to their target domains")
            o2 = OtuTableCollection()
            if args.metagenome_archive_otu_tables:
                o2.archive_table_objects.append(metagenomes.exclude_off_target_hits(pkgs, return_archive_table=True))
            else:
                o2.otu_table_objects.append(metagenomes.exclude_off_target_hits(pkgs))
            metagenomes = o2

        if args.genome_otu_tables or args.genome_archive_otu_tables:
            genomes = OtuTableCollection()
            if args.genome_otu_tables:
                for table in args.genome_otu_tables:
                    with open(table) as f:
                        genomes.add_otu_table(f)
            if args.genome_archive_otu_tables:
                for table in args.genome_archive_otu_tables:
                    with open(table) as f:
                        genomes.add_archive_otu_table(f)
            o2 = OtuTableCollection()
            if args.genome_archive_otu_tables:
                o2.archive_table_objects.append(genomes.exclude_off_target_hits(pkgs, return_archive_table=True))
            else:
                o2.otu_table_objects.append(genomes.exclude_off_target_hits(pkgs))
            genomes = o2
        else:
            genomes = None
        if args.assembly_otu_tables or args.assembly_archive_otu_tables:
            assemblies = OtuTableCollection()
            if args.assembly_otu_tables:
                for table in args.assembly_otu_tables:
                    with open(table) as f:
                        assemblies.add_otu_table(f)
            if args.assembly_archive_otu_tables:
                for table in args.assembly_archive_otu_tables:
                    with open(table) as f:
                        assemblies.add_archive_otu_table(f)
            o2 = OtuTableCollection()
            if args.assembly_archive_otu_tables:
                o2.archive_table_objects.append(assemblies.exclude_off_target_hits(pkgs, return_archive_table=True))
            else:
                o2.otu_table_objects.append(assemblies.exclude_off_target_hits(pkgs))
            assemblies = o2
        else:
            assemblies = None

        if genomes is None and assemblies is None:
            raise Exception("Appraise must be run with genomes and/or assemblies.")

        if args.output_binned_otu_table:
            output_binned_otu_table_io = open(args.output_binned_otu_table,'w')
        if args.output_unbinned_otu_table:
            output_unbinned_otu_table_io = open(args.output_unbinned_otu_table,'w')
        if args.output_assembled_otu_table:
            output_assembled_otu_table_io = open(args.output_assembled_otu_table,'w')
        if args.output_unaccounted_for_otu_table:
            output_unaccounted_for_otu_table_io = open(args.output_unaccounted_for_otu_table,'w')

        if args.stream_inputs:
            if args.assembly_otu_tables or args.assembly_archive_otu_tables:
                raise Exception("Streaming inputs is not currently known to work with assembly OTU tables")
            if args.output_assembled_otu_table:
                raise Exception("Streaming inputs is not currently known to work with --output-assembled-otu-table")
            if not (args.output_binned_otu_table and (args.output_unbinned_otu_table or args.output_unaccounted_for_otu_table)):
                raise Exception("--stream-inputs requires --output-binned-otu-table and --output-unbinned-otu-table to be defined")
            if args.output_unaccounted_for_otu_table and args.output_unbinned_otu_table:
                raise Exception("Cannot specify both --output-unaccounted-for-otu-table and --output-unbinned-otu-table")
            if args.output_unaccounted_for_otu_table:
                output_unbinned_otu_table_io = output_unaccounted_for_otu_table_io

            appraiser.streaming_appraise(
                genome_otu_table_collection=genomes,
                metagenome_otu_table_collection=metagenomes,
                output_found_in=args.output_found_in,
                sequence_identity=(args.sequence_identity if args.imperfect else None),
                window_size=DEFAULT_WINDOW_SIZE,
                binned_otu_table_io=output_binned_otu_table_io if args.output_binned_otu_table else None,
                unbinned_otu_table_io=output_unbinned_otu_table_io if (args.output_unbinned_otu_table or args.output_unaccounted_for_otu_table) else None,
                threads=args.threads,
            )
            for f in file_io:
                f.close()
        else:
            app = appraiser.appraise(genome_otu_table_collection=genomes,
                                    metagenome_otu_table_collection=metagenomes,
                                    assembly_otu_table_collection=assemblies,
                                    output_found_in = args.output_found_in,
                                    sequence_identity=(args.sequence_identity if args.imperfect else None),
                                    packages=pkgs,
                                    window_size=DEFAULT_WINDOW_SIZE)

            if args.plot_basename or args.plot:
                if args.plot and args.plot_basename:
                    raise Exception("Cannot specify both --plot and --plot-basename")
                if args.plot:
                    app.plot(
                        cluster_identity=args.sequence_identity,
                        doing_assembly=assemblies is not None,
                        doing_binning=genomes is not None,
                        gene_to_plot=args.plot_marker,
                        output_svg=args.plot)
                else:
                    app.plot(
                        output_svg_base=args.plot_basename,
                        cluster_identity=args.sequence_identity,
                        doing_assembly=assemblies is not None,
                        doing_binning=genomes is not None)

            appraiser.print_appraisal(
                app,
                packages=pkgs,
                doing_binning = genomes is not None,
                doing_assembly = assemblies is not None,
                output_found_in = args.output_found_in,
                output_style = args.output_style,
                binned_otu_table_io=output_binned_otu_table_io if args.output_binned_otu_table else None,
                unbinned_otu_table_io=output_unbinned_otu_table_io if args.output_unbinned_otu_table else None,
                assembled_otu_table_io=output_assembled_otu_table_io if args.output_assembled_otu_table else None,
                unaccounted_for_otu_table_io=output_unaccounted_for_otu_table_io \
                if args.output_unaccounted_for_otu_table else None)

        if args.output_binned_otu_table: output_binned_otu_table_io.close()
        if args.output_unbinned_otu_table: output_unbinned_otu_table_io.close()
        if args.output_assembled_otu_table: output_assembled_otu_table_io.close()
        if args.output_unaccounted_for_otu_table: output_unaccounted_for_otu_table_io.close()

    elif args.subparser_name == 'regenerate':
        from singlem.regenerator import Regenerator
        if not args.no_candidate_decoy_sequences and (not args.candidate_decoy_sequences or not args.candidate_decoy_taxonomy):
            raise Exception("Either --no-further-euks or euk taxonomy and sequences must be specified")
        Regenerator().regenerate(
            input_singlem_package = args.input_singlem_package,
            output_singlem_package = args.output_singlem_package,
            input_sequences = args.input_sequences,
            input_taxonomy = args.input_taxonomy,
            euk_sequences = args.candidate_decoy_sequences,
            euk_taxonomy = args.candidate_decoy_taxonomy,
            window_position = args.window_position,
            sequence_prefix = args.sequence_prefix,
            min_aligned_percent = args.min_aligned_percent,
            no_further_euks = args.no_candidate_decoy_sequences)

    elif args.subparser_name == 'get_tree':
        from singlem.metapackage import Metapackage

        hmm_db = Metapackage(args.singlem_packages)
        print("marker\ttree_file")
        for pkg in hmm_db.singlem_packages:
            print("{}\t{}".format(
                pkg.graftm_package_basename(),
                os.path.abspath(
                    pkg.graftm_package().reference_package_tree_path())))

    elif args.subparser_name == 'chainsaw':
        from singlem.chainsaw import Chainsaw
        Chainsaw().chainsaw(
            input_singlem_package_path = args.input_singlem_package,
            output_singlem_package_path = args.output_singlem_package,
            sequence_prefix = args.sequence_prefix,
            keep_tree = args.keep_tree)

    elif args.subparser_name == 'metapackage':
        from singlem.metapackage import Metapackage
        if args.describe:
            if args.metapackage:
                mpkg = Metapackage.acquire(args.metapackage)
            else:
                mpkg = Metapackage.acquire_default()
            mpkg.describe()
        else:
            if not args.metapackage:
                raise Exception("Must specify --metapackage as a path of the metapackage to create")
            if not args.singlem_packages:
                raise Exception("Must specify at least one singlem package to create a metapackage")
            if not args.nucleotide_sdb and not args.no_nucleotide_sdb:
                raise Exception("Must specify --nucleotide-sdb or --no-nucleotide-sdb")
            if not args.taxon_genome_lengths and not args.no_taxon_genome_lengths:
                raise Exception("Must specify --taxon-genome-lengths or --no-taxon-genome-lengths")
            if args.taxonomy_database_version and not args.taxonomy_database_name:
                raise Exception("--taxonomy-name must be specified if --taxonomy-version is specified")
            Metapackage.generate(
                singlem_packages = args.singlem_packages,
                nucleotide_sdb = args.nucleotide_sdb,
                taxon_genome_lengths = args.taxon_genome_lengths,
                output_path = args.metapackage,
                threads = args.threads,
                prefilter_clustering_threshold = args.prefilter_clustering_threshold,
                prefilter_diamond_db = args.prefilter_diamond_db,
                taxonomy_database_name = args.taxonomy_database_name,
                taxonomy_database_version = args.taxonomy_database_version,
                diamond_prefilter_performance_parameters = args.diamond_prefilter_performance_parameters,
                diamond_taxonomy_assignment_performance_parameters = args.diamond_taxonomy_assignment_performance_parameters,
                makeidx_sensitivity_params = args.makeidx_sensitivity_params,
                calculate_average_num_genes_per_species = args.calculate_average_num_genes_per_species,
            )

    elif args.subparser_name == 'condense':
        # from singlem.otu_table_collection import StreamingOtuTableCollection
        from singlem.condense import Condenser

        otus = generate_streaming_otu_table_from_args(args, input_prefix=True, archive_only=True, min_archive_otu_table_version=4)
        if not args.taxonomic_profile and not args.taxonomic_profile_krona:
            raise Exception("Either a krona or OTU table output must be specified for condense.")
        Condenser().condense(
            input_streaming_otu_table = otus,
            viral_mode = False,
            metapackage_path = args.metapackage,
            trim_percent = args.trim_percent,
            output_otu_table = args.taxonomic_profile,
            krona = args.taxonomic_profile_krona,
            min_taxon_coverage = args.min_taxon_coverage,
            output_after_em_otu_table = args.output_after_em_otu_table)

    elif args.subparser_name == 'trim_package_hmms':
        from singlem.trim_package_hmms import PackageHmmTrimmer

        PackageHmmTrimmer().trim(
            args.input_singlem_package,
            args.output_singlem_package,
        )

    elif args.subparser_name == 'seqs':
        seqs(args)

    elif args.subparser_name=='makedb':
        # Import here to avoid the slow tensorflow import in other subcommands
        from singlem.sequence_database import SequenceDatabase
        # from singlem.otu_table_collection import StreamingOtuTableCollection

        if args.pregenerated_otu_sqlite_db:
            otus = None
        else:
            otus = generate_streaming_otu_table_from_args(args)

        sequence_database_methods = args.sequence_database_methods
        if 'none' in sequence_database_methods:
            if len(sequence_database_methods) != 1:
                raise Exception("--sequence-database-methods does not except both 'none' and any another method at the same time")
            else:
                sequence_database_methods = []

        singlem.sequence_database.SequenceDatabase.create_from_otu_table(
            args.db,
            otus,
            num_threads=args.threads,
            pregenerated_sqlite3_db=args.pregenerated_otu_sqlite_db,
            num_annoy_nucleotide_trees=args.num_annoy_nucleotide_trees,
            num_annoy_protein_trees=args.num_annoy_protein_trees,
            tmpdir = args.tmpdir,
            sequence_database_methods = sequence_database_methods,
            sequence_database_types = args.sequence_database_types)
    elif args.subparser_name=='query':
        # Import here to avoid the slow tensorflow import in other subcommands
        from singlem.querier import Querier
        from singlem.sequence_database import SequenceDatabase
        # from singlem.otu_table_collection import StreamingOtuTableCollection
        querier = Querier()

        if args.dump:
            logging.debug("Dumping SingleM database {}".format(args.db))
            SequenceDatabase.dump(args.db)
        elif args.sample_names or args.taxonomy or args.sample_list:
            definition_type_count = 0
            if args.sample_names: definition_type_count += 1
            if args.taxonomy: definition_type_count += 1
            if args.sample_list: definition_type_count += 1
            if definition_type_count > 1:
                raise Exception("Only one of --sample-names, --sample-list and --taxonomy can be specified")
            sample_names = args.sample_names
            if args.sample_list:
                with open(args.sample_list) as f:
                    sample_names = list([line.strip() for line in f])
                    logging.info("Read in {} sample names from {}".format(len(sample_names), args.sample_list))
            querier.print_samples(
                db=args.db,
                sample_names = sample_names,
                taxonomy = args.taxonomy,
                output_io=sys.stdout)
        else:
            if args.db:
                querier.query(
                    db = args.db,
                    max_divergence = args.max_divergence,
                    output_style = 'sparse',#args.otu_table_type, # There is only sparse type atm.
                    query_otu_table = generate_streaming_otu_table_from_args(args, query_prefix=True),
                    num_threads = args.threads,
                    search_method = args.search_method,
                    sequence_type = args.sequence_type,
                    # stream_output = args.stream_output,
                    max_nearest_neighbours = args.max_nearest_neighbours,
                    max_search_nearest_neighbours = args.max_search_nearest_neighbours,
                    preload_db = args.preload_db,
                    limit_per_sequence = args.limit_per_sequence,
                    continue_on_missing_genes = args.continue_on_missing_genes)

    elif args.subparser_name=='data':
        from singlem.metapackage import Metapackage
        if args.verify_only:
            Metapackage.verify(output_directory = args.output_directory)
        else:
            Metapackage.download(
                output_directory = args.output_directory,
            )

    elif args.subparser_name in ('prokaryotic_fraction', 'microbial_fraction'):
        from singlem.read_fraction import ReadFractionEstimator
        ReadFractionEstimator().calculate_and_report_read_fraction(
            input_profile = args.input_profile,
            metagenome_sizes = args.input_metagenome_sizes,
            forward_read_files = args.forward,
            reverse_read_files = args.reverse,
            taxonomic_genome_lengths_file = args.taxon_genome_lengths_file,
            metapackage = args.metapackage,
            accept_missing_samples = args.accept_missing_samples,
            output_tsv = args.output_tsv,
            output_per_taxon_read_fractions = args.output_per_taxon_read_fractions,)

    elif args.subparser_name == 'supplement':
        from singlem.supplement import Supplementor
        if not args.no_taxon_genome_lengths and not args.checkm2_quality_file:
            raise Exception("Either --no-taxon-genome-lengths or --checkm2-quality-file must be specified")
        if not args.checkm2_quality_file and not args.no_quality_filter:
            raise Exception("Either --no-quality-filter or --checkm2-quality-file must be specified")
        if args.gtdbtk_output_directory and args.taxonomy_file:
            raise Exception("Options --gtdbtk-output-directory and --taxonomy-file cannot both be specified")

        Supplementor().supplement(
            new_genome_fasta_files=args.new_genome_fasta_files,
            new_genome_fasta_files_list=args.new_genome_fasta_files_list,
            new_taxonomies=args.new_fully_defined_taxonomies,
            input_metapackage=args.input_metapackage,
            output_metapackage=args.output_metapackage,
            threads=args.threads,
            pplacer_threads=args.pplacer_threads,
            working_directory=args.working_directory,
            gtdbtk_output_directory=args.gtdbtk_output_directory,
            taxonomy_file=args.taxonomy_file,
            output_taxonomies=args.output_taxonomies,
            checkm2_quality_file=args.checkm2_quality_file,
            no_quality_filter=args.no_quality_filter,
            checkm2_min_completeness=args.checkm2_min_completeness,
            checkm2_max_contamination=args.checkm2_max_contamination,
            hmmsearch_evalue=args.hmmsearch_evalue,
            no_dereplication=args.no_dereplication,
            dereplicate_with_galah=args.dereplicate_with_galah,
            skip_taxonomy_check=args.skip_taxonomy_check,
            no_taxon_genome_lengths=args.no_taxon_genome_lengths,
            gene_definitions=args.gene_definitions,
            ignore_taxonomy_database_incompatibility=args.ignore_taxonomy_database_incompatibility,
            new_taxonomy_database_name=args.new_taxonomy_database_name,
            new_taxonomy_database_version=args.new_taxonomy_database_version,
        )

    else:
        raise Exception("Programming error")

if __name__ == '__main__':
    main()
