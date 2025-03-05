import logging
import tempfile

from .pipe import SearchPipe
from .metapackage import Metapackage
from .archive_otu_table import ArchiveOtuTable
from .pipe_sequence_extractor import ExtractedReads, ExtractedReadSet
from .sequence_classes import Sequence, UnalignedAlignedNucleotideSequence
from .otu_table_collection import StreamingOtuTableCollection

class Renew:
    @staticmethod
    def renew(**kwargs):
        '''Renew an OTU table, annotating old OTU tables with new taxonomy info'''
        input_archive_otu_table = kwargs.pop('input_archive_otu_table')
        output_archive_otu_table = kwargs.pop('output_archive_otu_table')
        output_otu_table = kwargs.pop('otu_table')
        output_extras = kwargs.pop('output_extras')
        threads = kwargs.pop('threads')
        singlem_assignment_method = kwargs.pop('assignment_method')
        output_jplace = kwargs.pop('output_jplace')
        diamond_taxonomy_assignment_performance_parameters = kwargs.pop('diamond_taxonomy_assignment_performance_parameters')
        evalue = kwargs.pop('evalue')
        min_orf_length = kwargs.pop('min_orf_length')
        restrict_read_length = kwargs.pop('restrict_read_length')
        translation_table = kwargs.pop('translation_table')
        filter_minimum_protein = kwargs.pop('filter_minimum_protein')
        assignment_singlem_db = kwargs.pop('assignment_singlem_db')
        output_taxonomic_profile = kwargs.pop('output_taxonomic_profile')
        output_taxonomic_profile_krona = kwargs.pop('output_taxonomic_profile_krona')
        exclude_off_target_hits = kwargs.pop('exclude_off_target_hits')
        max_species_divergence = kwargs.pop('max_species_divergence')
        ignore_missing_singlem_packages = kwargs.pop('ignore_missing_singlem_packages')

        logging.info("Acquiring singlem packages ..")
        metapackage = SearchPipe()._parse_packages_or_metapackage(**kwargs)
        kwargs.pop('singlem_packages', None)
        kwargs.pop('metapackage_path', None)

        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        kwargs['metapackage_object'] = metapackage

        if assignment_singlem_db is None:
            assignment_singlem_db = metapackage.nucleotide_sdb_path()
        
        # Extract diamond_taxonomy_assignment_performance_parameters from metapackage (v5 metapackages only)
        if diamond_taxonomy_assignment_performance_parameters == None:
            diamond_taxonomy_assignment_performance_parameters = metapackage.diamond_taxonomy_assignment_performance_parameters()
            if diamond_taxonomy_assignment_performance_parameters == None:
                diamond_taxonomy_assignment_performance_parameters = SearchPipe.DEFAULT_DIAMOND_ASSIGN_TAXONOMY_PERFORMANCE_PARAMETERS

        # Create a dict of singlem package name (i.e. the string in the OTU
        # table) to singlem package object
        marker_name_to_spkg = {}
        for spkg in metapackage:
            marker_name = spkg.graftm_package_basename().replace('.gpkg','')
            if marker_name in marker_name_to_spkg:
                raise Exception("Duplicate packages for marker name {}".format(marker_name))
            marker_name_to_spkg[marker_name] = spkg
        logging.info("Acquired {} singlem packages with name e.g. {}".format(
            len(marker_name_to_spkg),
            list(marker_name_to_spkg.keys())[0]))

        # Read in archive OTU table, requiring a minimum version
        #
        # TODO: Probably this can converted to a streaming input, following
        # https://stackoverflow.com/questions/54560154/streaming-json-parser
        logging.info("Reading in archive OTU table ..")
        with open(input_archive_otu_table) as f:
            input_otus = ArchiveOtuTable.read(f)
        if input_otus.version < 2:
            raise Exception("Currently only version 2+ archive otu tables are supported")
        logging.info("Read in {} OTUs".format(len(input_otus.data)))

        # Sort the archive OTU table, because sometimes the OTUs are not in
        # marker-wise order, which causes repeated calls to DIAMOND with too few
        # seqs. Only really need marker and sample name to be together.
        input_otus.sort()
        
        # Generate ExtractedReads. Never analysing pairs because no 2 reads with
        # the same name should be in the same OTU.
        extracted_reads = ExtractedReads(False)

        class State:
            def __init__(self, marker_name_to_spkg):
                self.marker_name_to_spkg = marker_name_to_spkg
                self.reset()

            def reset(self):
                self.current_singlem_package = None
                self.current_sample_name = None
                self.current_sequences = []
                self.current_unaligned_aligned_nuc_seqs = []

        last_marker_and_sample = None
        state = State(marker_name_to_spkg)

        def process_otu_batch(state):
            extracted_reads.add(ExtractedReadSet(
                state.current_sample_name,
                state.marker_name_to_spkg[state.current_marker_name],
                state.current_sequences,
                [],
                state.current_unaligned_aligned_nuc_seqs))

        read_unaligned_sequences_field = ArchiveOtuTable.FIELDS_VERSION2.index('read_unaligned_sequences')
        nucleotides_aligned_field = ArchiveOtuTable.FIELDS_VERSION2.index('nucleotides_aligned')

        num_ignored_otus = 0
        for otu in input_otus:
            # Ensure that the packages in the archive exist
            if not otu.marker in marker_name_to_spkg:
                if ignore_missing_singlem_packages:
                    logging.debug("Ignoring marker '{}' as it was not found in the specified singlem packages".format(otu.marker))
                    num_ignored_otus += 1
                    continue
                else:
                    raise Exception("Found marker '{}' that was not one of the specified singlem packages".format(otu.marker))

            marker_and_sample = [otu.marker, otu.sample_name]
            if marker_and_sample != last_marker_and_sample:
                if last_marker_and_sample is not None:
                    process_otu_batch(state)
                    state.reset()
                last_marker_and_sample = marker_and_sample

            state.current_marker_name = otu.marker
            state.current_sample_name = otu.sample_name

            read_names = otu.read_names()
            seqs = otu.data[read_unaligned_sequences_field]
            nucleotides_aligned = otu.data[nucleotides_aligned_field]
            if len(read_names) != len(nucleotides_aligned) or len(seqs) != len(nucleotides_aligned):
                raise Exception("Unexpected format of otu found: {}".format(otu))
            for (name, seq, num_aligned) in zip(read_names, seqs, nucleotides_aligned):
                state.current_sequences.append(Sequence(name, seq))
                state.current_unaligned_aligned_nuc_seqs.append(
                    UnalignedAlignedNucleotideSequence(
                        name, None, otu.sequence, seq, num_aligned))
        if ignore_missing_singlem_packages:
            logging.info("Ignored {} OTUs that were not found in the specified singlem packages".format(num_ignored_otus))
        
        # process last batch
        if last_marker_and_sample is not None:
            process_otu_batch(state)

        # To save RAM (untested)
        del input_otus

        # Assign taxonomy and process
        pipe = SearchPipe()
        with tempfile.TemporaryDirectory(prefix='singlem-renew') as td:
            pipe._working_directory = td #FIXME, shouldn't be referencing underscore variables.
            pipe._graftm_verbosity = '5' if logging.getLevelName(logging.getLogger().level) == 'DEBUG' else '2'
            pipe._evalue = evalue
            pipe._restrict_read_length = restrict_read_length
            pipe._min_orf_length = min_orf_length
            pipe._num_threads = threads
            pipe._filter_minimum_protein = filter_minimum_protein
            pipe._translation_table = translation_table
            pipe._max_species_divergence = max_species_divergence

            logging.info("Running taxonomy assignment and post-processing ..")
            otu_table_object = pipe.assign_taxonomy_and_process(
                extracted_reads=extracted_reads,
                analysing_pairs=False,
                assign_taxonomy=True,
                singlem_assignment_method=singlem_assignment_method,
                threads=threads,
                diamond_taxonomy_assignment_performance_parameters=diamond_taxonomy_assignment_performance_parameters,
                known_sequence_taxonomy=None,
                known_taxes=None,
                output_jplace=output_jplace,
                assignment_singlem_db=assignment_singlem_db,
            )

        # Write outputs
        if output_otu_table or output_archive_otu_table:
            logging.info("Writing output files ..")
            pipe.write_otu_tables(
                otu_table_object,
                output_otu_table,
                output_archive_otu_table,
                output_extras,
                metapackage,
                exclude_off_target_hits)
        
        if output_taxonomic_profile or output_taxonomic_profile_krona:
            from .condense import Condenser
            otu_table_collection = StreamingOtuTableCollection()
            otu_table_collection.add_archive_otu_table_object(otu_table_object)
            Condenser().condense(
                input_streaming_otu_table = otu_table_collection,
                output_otu_table = output_taxonomic_profile,
                krona = output_taxonomic_profile_krona,
                metapackage = metapackage)

        logging.info("Renew is finished")
