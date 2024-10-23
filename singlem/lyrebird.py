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

def main():
    bird_argparser = BirdArgparser(
        program='Lyrebird',
        authors=[
            BEN_NAME_AND_CENTRE,
            "Samuel Aroney, "+CMR,
            "Raphael Eisenhofer, Centre for Evolutionary Hologenomics, University of Copenhagen, Denmark",
            "Rossen Zhao, "+CMR],
        version=singlem.__version__, # TODO: change to lyrebird version 0.1.0
        raw_format=True,
        examples={'pipe': [
            Example(
                'Get a taxonomic profile for dsDNA phages from paired read input:',
                'lyrebird pipe -1 <fastq_or_fasta1> -2 <fastq_or_fasta2> -p <output.profile.tsv>'),
            Example(
                'Get a taxonomic profile Krona diagram for dsDNA phages from single read input:',
                'lyrebird pipe -i <fastq_or_fasta> --taxonomic-profile-krona <output.profile.html>'),
            Example(
                'Gather an OTU table (per marker sequence groupings) from paired reads:',
                'lyrebird pipe -1 <fastq_or_fasta1> -2 <fastq_or_fasta2> --otu-table <output.otu_table.tsv>'),
        ]},
    )

    data_description = 'Download reference metapackage data'
    data_parser = bird_argparser.new_subparser('data', data_description, parser_group='Tools')
    # TODO: Could make pipe invocation faster by moving DATA_ENVIRONMENT to a separate file
    data_parser.add_argument('--output-directory', help="Output directory [required unless {} is specified]".format(DATA_ENVIRONMENT_VARIABLE))
    data_parser.add_argument('--verify-only', help="Check that the data is up to date and each file has the correct checksum", action='store_true', default=False)

    pipe_description = 'Generate a taxonomic profile or OTU table for dsDNA phages from raw sequences'
    pipe_parser = bird_argparser.new_subparser('pipe', pipe_description, parser_group='Tools')

    # Make a function here so the code can be re-used between pipe and renew
    def add_common_pipe_arguments(argument_group):
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
    common_pipe_arguments = pipe_parser.add_argument_group('Common options')
    sequence_input_group = common_pipe_arguments.add_mutually_exclusive_group(required=True)
    # Keep parity of these arguments with the 'read_fraction' command
    sequence_input_group.add_argument('-1','--forward','--reads','--sequences',
                                nargs='+',
                                metavar='sequence_file',
                                help='nucleotide read sequence(s) (forward or unpaired) to be searched. Can be FASTA or FASTQ format, GZIP-compressed or not.')
    common_pipe_arguments.add_argument('-2', '--reverse',
                                nargs='+',
                                metavar='sequence_file',
                                help='reverse reads to be searched. Can be FASTA or FASTQ format, GZIP-compressed or not.')
    sequence_input_group.add_argument('--genome-fasta-files',
                                nargs='+',
                                metavar='sequence_file',
                                help='nucleotide genome sequence(s) to be searched')
    sequence_input_group.add_argument('--sra-files',
                                nargs='+',
                                metavar='sra_file',
                                help='"sra" format files (usually from NCBI SRA) to be searched')
    add_common_pipe_arguments(common_pipe_arguments)

    def add_less_common_pipe_arguments(argument_group):
        argument_group.add_argument('--archive-otu-table', metavar='filename', help='output OTU table in archive format for making DBs etc. [default: unused]')
        argument_group.add_argument('--output-jplace', metavar='filename', help='Output a jplace format file for each singlem package to a file starting with this string, each with one entry per OTU. Requires \'%s\' as the --assignment_method [default: unused]' % pipe.PPLACER_ASSIGNMENT_METHOD)
        argument_group.add_argument('--metapackage', help='Set of SingleM packages to use [default: use the default set]')
        argument_group.add_argument('--singlem-packages', nargs='+', help='SingleM packages to use [default: use the set from the default metapackage]')
        argument_group.add_argument('--assignment-singlem-db', '--assignment_singlem_db', help='Use this SingleM DB when assigning taxonomy [default: not set, use the default]')
        argument_group.add_argument('--diamond-taxonomy-assignment-performance-parameters',
                                    help='Performance-type arguments to use when calling \'diamond blastx\' during the taxonomy assignment step. [default: use setting defined in metapackage when set, otherwise use \'%s\']' % SearchPipe.DEFAULT_DIAMOND_ASSIGN_TAXONOMY_PERFORMANCE_PARAMETERS,
                                    default=None)
        argument_group.add_argument('--evalue',
                                    help='HMMSEARCH e-value cutoff to use for sequence gathering [default: %s]' % SearchPipe.DEFAULT_HMMSEARCH_EVALUE, default=SearchPipe.DEFAULT_HMMSEARCH_EVALUE)
        argument_group.add_argument('--min-orf-length',
                                    metavar='length',
                                    help='When predicting ORFs require this many base pairs uninterrupted by a stop codon [default: %i when input is reads, %i when input is genomes]' % (SearchPipe.DEFAULT_MIN_ORF_LENGTH, SearchPipe.DEFAULT_GENOME_MIN_ORF_LENGTH),
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

    less_common_pipe_arguments = pipe_parser.add_argument_group('Less common options')
    add_less_common_pipe_arguments(less_common_pipe_arguments)

    less_common_pipe_arguments.add_argument('--working-directory', metavar='directory', help='use intermediate working directory at a specified location, and do not delete it upon completion [default: not set, use a temporary directory]')
    less_common_pipe_arguments.add_argument('--working-directory-dev-shm', default=False, action='store_true', help='use an intermediate results temporary working directory in /dev/shm rather than the default [default: the usual temporary working directory, currently {}]'.format(
        tempfile.gettempdir()
    ))
    less_common_pipe_arguments.add_argument('--force', action='store_true', help='overwrite working directory if required [default: not set]')
    less_common_pipe_arguments.add_argument('--filter-minimum-nucleotide',
                                metavar='length',
                                help='Ignore reads aligning in less than this many positions to each nucleotide HMM [default: %i]' % SearchPipe.DEFAULT_FILTER_MINIMUM_NUCLEOTIDE,
                                type=int, default=SearchPipe.DEFAULT_FILTER_MINIMUM_NUCLEOTIDE)
    less_common_pipe_arguments.add_argument('--include-inserts', action='store_true',
                                help='print the entirety of the sequences in the OTU table, not just the aligned nucleotides [default: not set]', default=False)
    less_common_pipe_arguments.add_argument('--known-otu-tables', nargs='+',
                                help='OTU tables previously generated with trusted taxonomies for each sequence [default: unused]')
    less_common_pipe_arguments.add_argument('--no-assign-taxonomy', action='store_true',
                                help='Do not assign any taxonomy except for those already known [default: not set]',
                                default=False)
    less_common_pipe_arguments.add_argument('--known-sequence-taxonomy', metavar='FILE',
                                help='A 2-column "sequence<tab>taxonomy" file specifying some sequences that have known taxonomy [default: unused]')
    less_common_pipe_arguments.add_argument('--no-diamond-prefilter', action='store_true',
                                help='Do not parse sequence data through DIAMOND blastx using a database constructed from the set of singlem packages. Should be used with --hmmsearch-package-assignment. NOTE: ignored for nucleotide packages [default: protein packages: use the prefilter, nucleotide packages: do not use the prefilter]',
                                default=False)
    less_common_pipe_arguments.add_argument('--diamond-prefilter-performance-parameters',
                                help='Performance-type arguments to use when calling \'diamond blastx\' during the prefiltering. By default, SingleM should run in <4GB of RAM except in very large (>100Gbp) metagenomes. [default: use setting defined in metapackage when set, otherwise use \'%s\']' % SearchPipe.DEFAULT_PREFILTER_PERFORMANCE_PARAMETERS,
                                default=None)
    less_common_pipe_arguments.add_argument('--hmmsearch-package-assignment', '--hmmsearch_package_assignment', action='store_true',
                                help='Assign each sequence to a SingleM package using HMMSEARCH, and a sequence may then be assigned to multiple packages. [default: not set]',
                                default=False)
    less_common_pipe_arguments.add_argument('--diamond-prefilter-db',
                                help='Use this DB when running DIAMOND prefilter [default: use the one in the metapackage, or generate one from the SingleM packages]')
    less_common_pipe_arguments.add_argument('--assignment-threads',type=int,
                                help='Use this many processes in parallel while assigning taxonomy [default: %i]' % SearchPipe.DEFAULT_ASSIGNMENT_THREADS,
                                default=SearchPipe.DEFAULT_ASSIGNMENT_THREADS)
    less_common_pipe_arguments.add_argument('--sleep-after-mkfifo', type=int,
                                help='Sleep for this many seconds after running os.mkfifo [default: None]')

    renew_description = 'Reannotate an OTU table with an updated taxonomy'
    renew_parser = bird_argparser.new_subparser('renew', renew_description, parser_group='Tools')
    renew_input_args = renew_parser.add_argument_group('input')
    renew_input_args.add_argument('--input-archive-otu-table', help="Renew this table", required=True)
    renew_common = renew_parser.add_argument_group("Common arguments in shared with 'pipe'")
    add_common_pipe_arguments(renew_common)
    renew_less_common = renew_parser.add_argument_group("Less common arguments shared with 'pipe'")
    add_less_common_pipe_arguments(renew_less_common)

    condense_description = 'Combine OTU tables across different markers into a single taxonomic profile. Modified for non-universal markers and requires a Lyrebird metapackage. Note that while this mode can be run independently, it is often more straightforward to invoke its methodology by specifying -p / --taxonomic-profile when running pipe mode.'
    condense_parser = bird_argparser.new_subparser('condense', condense_description)

    input_condense_arguments = condense_parser.add_argument_group("Input arguments (1+ required)")

    input_condense_arguments.add_argument('--input-archive-otu-tables', '--input-archive-otu-table', nargs='+', help="Condense from these archive tables")
    input_condense_arguments.add_argument('--input-archive-otu-table-list',
        help="Condense from the archive tables newline separated in this file")
    input_condense_arguments.add_argument('--input-gzip-archive-otu-table-list',
        help="Condense from the gzip'd archive tables newline separated in this file")

    output_condense_arguments = condense_parser.add_argument_group("Output arguments (1+ required)")
    output_condense_arguments.add_argument('-p', '--taxonomic-profile', metavar='filename', help="output OTU table")
    output_condense_arguments.add_argument('--taxonomic-profile-krona', metavar='filename', help='name of krona file to generate.')
    output_condense_arguments.add_argument('--output-after-em-otu-table', metavar='filename', help="output OTU table after expectation maximisation has been applied. Note that this table usually contains multiple rows with the same window sequence.")

    optional_condense_arguments = condense_parser.add_argument_group("Other options")
    optional_condense_arguments.add_argument('--metapackage', help='Set of SingleM packages to use [default: use the default set]')
    current_default = CONDENSE_DEFAULT_MIN_TAXON_COVERAGE
    optional_condense_arguments.add_argument('--min-taxon-coverage',metavar='FRACTION',
        help='Set taxons with less coverage to coverage=0. [default: {}]'.format(current_default), default=current_default, type=float)
    current_default = CONDENSE_DEFAULT_TRIM_PERCENT
    optional_condense_arguments.add_argument('--trim-percent', type=float, default=current_default, help="percentage of markers to be trimmed for each taxonomy [default: {}]".format(current_default))

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

    args = bird_argparser.parse_the_args()

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

    def get_min_orf_length(args, subparser='pipe'):
        if args.min_orf_length:
            return args.min_orf_length
        elif subparser=='pipe' and (args.forward or args.sra_files):
            return SearchPipe.DEFAULT_MIN_ORF_LENGTH
        elif subparser=='pipe' and args.genome_fasta_files:
            return SearchPipe.DEFAULT_GENOME_MIN_ORF_LENGTH
        elif subparser=='renew':
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

    if args.subparser_name=='pipe':
        validate_pipe_args(args)
        singlem.pipe.SearchPipe().run(
            sequences = args.forward,
            reverse_read_files = args.reverse,
            genomes = args.genome_fasta_files,
            input_sra_files = args.sra_files,
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
            viral_profile_output = True,
            exclude_off_target_hits = args.exclude_off_target_hits,
            min_taxon_coverage = get_min_taxon_coverage(args),
            max_species_divergence = args.max_species_divergence,
        )

    elif args.subparser_name=='renew':
        from singlem.renew import Renew
        validate_pipe_args(args, subparser='renew')
        Renew().renew(
            input_archive_otu_table=args.input_archive_otu_table,
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
            max_species_divergence = args.max_species_divergence,
            )

    elif args.subparser_name == 'condense':
        from singlem.otu_table_collection import StreamingOtuTableCollection
        from singlem.condense import Condenser

        otus = generate_streaming_otu_table_from_args(args, input_prefix=True, archive_only=True, min_archive_otu_table_version=4)
        if not args.taxonomic_profile and not args.taxonomic_profile_krona:
            raise Exception("Either a krona or OTU table output must be specified for condense.")
        Condenser().condense(
            input_streaming_otu_table = otus,
            viral_mode = True,
            metapackage_path = args.metapackage,
            trim_percent = args.trim_percent,
            output_otu_table = args.taxonomic_profile,
            krona = args.taxonomic_profile_krona,
            min_taxon_coverage = args.min_taxon_coverage,
            output_after_em_otu_table = args.output_after_em_otu_table)

    elif args.subparser_name=='data':
        from singlem.metapackage import Metapackage
        if args.verify_only:
            Metapackage.verify(output_directory = args.output_directory)
        else:
            Metapackage.download(
                output_directory = args.output_directory,
            )

    else:
        raise Exception("Programming error")
