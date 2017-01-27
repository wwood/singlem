import tempdir
import logging
import os.path
import shutil
import extern
import itertools
import tempfile
import subprocess
import json
import re
from string import split

from singlem import HmmDatabase, TaxonomyFile, OrfMUtils
from otu_table import OtuTable
from known_otu_table import KnownOtuTable
from metagenome_otu_finder import MetagenomeOtuFinder
from sequence_classes import SeqReader
from diamond_parser import DiamondResultParser
from graftm_result import GraftMResult

from graftm.sequence_extractor import SequenceExtractor
from graftm.greengenes_taxonomy import GreenGenesTaxonomy
from graftm.sequence_search_results import HMMSearchResult, SequenceSearchResult

PPLACER_ASSIGNMENT_METHOD = 'pplacer'
DIAMOND_ASSIGNMENT_METHOD = 'diamond'
DIAMOND_EXAMPLE_BEST_HIT_ASSIGNMENT_METHOD = 'diamond_example'

class SearchPipe:
    def run(self, **kwargs):
        forward_read_files = kwargs.pop('sequences')
        output_otu_table = kwargs.pop('otu_table', None)
        archive_otu_table = kwargs.pop('archive_otu_table', None)
        num_threads = kwargs.pop('threads')
        known_otu_tables = kwargs.pop('known_otu_tables')
        singlem_assignment_method = kwargs.pop('assignment_method')
        output_jplace = kwargs.pop('output_jplace')
        output_extras = kwargs.pop('output_extras')
        evalue = kwargs.pop('evalue')
        min_orf_length = kwargs.pop('min_orf_length')
        restrict_read_length = kwargs.pop('restrict_read_length')
        filter_minimum_protein = kwargs.pop('filter_minimum_protein')
        filter_minimum_nucleotide = kwargs.pop('filter_minimum_nucleotide')
        include_inserts = kwargs.pop('include_inserts')
        singlem_packages = kwargs.pop('singlem_packages')
        window_size = kwargs.pop('window_size')
        assign_taxonomy = kwargs.pop('assign_taxonomy')
        known_sequence_taxonomy = kwargs.pop('known_sequence_taxonomy')

        working_directory = kwargs.pop('working_directory')
        force = kwargs.pop('force')
        if len(kwargs) > 0:
            raise Exception("Unexpected arguments detected: %s" % kwargs)
        
        self._num_threads = num_threads
        self._evalue = evalue
        self._min_orf_length = min_orf_length
        self._restrict_read_length = restrict_read_length
        self._filter_minimum_protein = filter_minimum_protein
        self._filter_minimum_nucleotide = filter_minimum_nucleotide

        hmms = HmmDatabase(singlem_packages)
        if singlem_assignment_method == DIAMOND_EXAMPLE_BEST_HIT_ASSIGNMENT_METHOD:
            graftm_assignment_method = DIAMOND_ASSIGNMENT_METHOD
        else:
            graftm_assignment_method = singlem_assignment_method
            
        if logging.getLevelName(logging.getLogger().level) == 'DEBUG':
            self._graftm_verbosity = '5'
        else:
            self._graftm_verbosity = '2'

        using_temporary_working_directory = working_directory is None
        if using_temporary_working_directory:
            shared_mem_directory = '/dev/shm'
            if os.path.exists(shared_mem_directory):
                logging.debug("Using shared memory as a base directory")
                tmp = tempdir.TempDir(basedir=shared_mem_directory)
                tempfiles_path = os.path.join(tmp.name, 'tempfiles')
                os.mkdir(tempfiles_path)
                os.environ['TEMP'] = tempfiles_path
            else:
                logging.debug("Shared memory directory not detected, using default temporary directory instead")
                tmp = tempdir.TempDir()
            working_directory = tmp.name
        else:
            working_directory = working_directory
            if os.path.exists(working_directory):
                if force:
                    logging.info("Overwriting directory %s" % working_directory)
                    shutil.rmtree(working_directory)
                    os.mkdir(working_directory)
                else:
                    raise Exception("Working directory '%s' already exists, not continuing" % working_directory)
            else:
                os.mkdir(working_directory)
        logging.debug("Using working directory %s" % working_directory)
        self._working_directory = working_directory

        extracted_reads = None
        def return_cleanly():
            if extracted_reads: extracted_reads.cleanup()
            if using_temporary_working_directory: tmp.dissolve()
            logging.info("Finished")

        #### Search
        self._singlem_package_database = hmms
        search_result = self._search(hmms, forward_read_files)
        sample_names = search_result.samples_with_hits()
        if len(sample_names) == 0:
            logging.info("No reads identified in any samples, stopping")
            return_cleanly()
            return
        logging.debug("Recovered %i samples with at least one hit e.g. '%s'" \
                     % (len(sample_names), sample_names[0]))

        #### Alignment
        align_result = self._align(search_result)

        ### Extract reads that have already known taxonomy
        if known_otu_tables:
            logging.info("Parsing known taxonomy OTU tables")
            known_taxes = KnownOtuTable()
            known_taxes.parse_otu_tables(known_otu_tables)
            logging.debug("Read in %i sequences with known taxonomy" % len(known_taxes))
        else:
            known_taxes = []
        if known_sequence_taxonomy:
            logging.debug("Parsing sequence-wise taxonomy..")
            tax1 = GreenGenesTaxonomy.read(open(known_sequence_taxonomy)).taxonomy
            known_sequence_tax = {}
            for seq_id, tax in tax1.items():
                known_sequence_tax[seq_id] = '; '.join(tax)
            logging.info("Read in %i taxonomies from the GreenGenes format taxonomy file" % len(known_sequence_tax))

        ### Extract other reads which do not have known taxonomy
        extracted_reads = self._extract_relevant_reads(
            align_result, include_inserts, known_taxes)
        logging.info("Finished extracting aligned sequences")

        #### Taxonomic assignment
        if assign_taxonomy:
            logging.info("Running taxonomic assignment with graftm..")
            assignment_result = self._assign_taxonomy(
                extracted_reads, graftm_assignment_method)

        #### Process taxonomically assigned reads
        # get the sequences out for each of them
        otu_table_object = OtuTable()
        regular_output_fields = split('gene sample sequence num_hits coverage taxonomy')
        otu_table_object.fields = regular_output_fields + \
                                  split('read_names nucleotides_aligned taxonomy_by_known?')

        for sample_name, singlem_package, tmp_graft, known_sequences, unknown_sequences in extracted_reads:
            def add_info(infos, otu_table_object, known_tax):
                for info in infos:
                    to_print = [
                        singlem_package.graftm_package_basename(),
                        sample_name,
                        info.seq,
                        info.count,
                        info.coverage,
                        info.taxonomy,
                        info.names,
                        info.aligned_lengths,
                        known_tax]
                    otu_table_object.data.append(to_print)
            known_infos = self._seqs_to_counts_and_taxonomy(
                known_sequences,
                known_taxes,
                False,
                True)
            add_info(known_infos, otu_table_object, True)
            
            if tmp_graft: # if any sequences were aligned (not just already known)
                tmpbase = os.path.basename(tmp_graft.name[:-6])#remove .fasta
                
                if assign_taxonomy:
                    is_known_taxonomy = False
                    aligned_seqs = self._get_windowed_sequences(
                        assignment_result.prealigned_sequence_file(
                            sample_name, singlem_package, tmpbase),
                        assignment_result.nucleotide_hits_file(
                            sample_name, singlem_package, tmpbase),
                        singlem_package,
                        include_inserts)
                    if singlem_assignment_method == DIAMOND_EXAMPLE_BEST_HIT_ASSIGNMENT_METHOD:
                        tax_file = assignment_result.diamond_assignment_file(
                            sample_name, singlem_package, tmpbase)
                    else:
                        tax_file = assignment_result.read_tax_file(
                            sample_name, singlem_package, tmpbase)
                    logging.debug("Reading taxonomy from %s" % tax_file)

                    if singlem_assignment_method == DIAMOND_EXAMPLE_BEST_HIT_ASSIGNMENT_METHOD:
                        taxonomies = DiamondResultParser(tax_file)
                        use_first = True
                    else:
                        if not os.path.isfile(tax_file):
                            logging.warn("Unable to find tax file for gene %s from sample %s "
                                         "(likely do to min length filtering), skipping" % (
                                             os.path.basename(singlem_package.base_directory()),
                                             sample_name))
                            taxonomies = {}
                        else:
                            taxonomies = TaxonomyFile(tax_file)
                        use_first = False
                        
                else: # Taxonomy has not been assigned.
                    aligned_seqs = unknown_sequences
                    if known_sequence_taxonomy:
                        taxonomies = known_sequence_tax
                    else:
                        taxonomies = {}
                    use_first = False # irrelevant
                    is_known_taxonomy = True
                    
                new_infos = list(self._seqs_to_counts_and_taxonomy(
                    aligned_seqs, taxonomies, use_first, False))
                add_info(new_infos, otu_table_object, is_known_taxonomy)
                
                if output_jplace:
                    base_dir = assignment_result._base_dir(
                        sample_name, singlem_package, tmpbase)
                    input_jplace_file = os.path.join(base_dir, "placements.jplace")
                    output_jplace_file = os.path.join(base_dir, "%s_%s_%s.jplace" % (
                        output_jplace, sample_name, singlem_package.graftm_package_basename()))
                    logging.debug("Converting jplace file %s to singlem jplace file %s" % (
                        input_jplace_file, output_jplace_file))
                    with open(output_jplace_file, 'w') as output_jplace_io:
                        self._write_jplace_from_infos(
                            open(input_jplace_file), new_infos, output_jplace_io)

                            
        if output_otu_table:
            with open(output_otu_table, 'w') as f:
                if output_extras:
                    otu_table_object.write_to(f, otu_table_object.fields)
                else:
                    otu_table_object.write_to(f, regular_output_fields)
        if archive_otu_table:
            with open(archive_otu_table, 'w') as f:
                otu_table_object.archive(hmms.singlem_packages).write_to(f)
        return_cleanly()

    def _get_windowed_sequences(self, protein_sequences_file, nucleotide_sequence_file,
                                singlem_package, include_inserts):
        if not os.path.exists(nucleotide_sequence_file) or \
            os.stat(nucleotide_sequence_file).st_size == 0: return []
        nucleotide_sequences = SeqReader().read_nucleotide_sequences(nucleotide_sequence_file)
        protein_alignment = self._align_proteins_to_hmm(
            protein_sequences_file,
            singlem_package.graftm_package().alignment_hmm_path())
        return MetagenomeOtuFinder().find_windowed_sequences(
            protein_alignment,
            nucleotide_sequences,
            60,
            include_inserts,
            singlem_package.is_protein_package(),
            best_position=singlem_package.singlem_position())

    def _extract_relevant_reads(self, alignment_result, include_inserts, known_taxonomy):
        '''Given a SingleMPipeAlignSearchResult, extract reads that will be used as
        part of the singlem choppage process as tempfiles in a hash'''
        
        extracted_reads = ExtractedReads(self._singlem_package_database)
        
        for sample_name in alignment_result.sample_names():
            for singlem_package in self._singlem_package_database:
                for prealigned_file in alignment_result.prealigned_sequence_files(
                        sample_name, singlem_package):
                    if singlem_package.is_protein_package():
                        nucleotide_sequence_fasta_file = alignment_result.nucleotide_sequence_file(
                            sample_name, singlem_package)
                    else:
                        nucleotide_sequence_fasta_file = prealigned_file

                    logging.debug("Aligning %s to the HMM" % prealigned_file)
                    aligned_seqs = self._get_windowed_sequences(
                        prealigned_file,
                        nucleotide_sequence_fasta_file,
                        singlem_package,
                        include_inserts)

                    known_sequences = []
                    unknown_sequences = []
                    for s in aligned_seqs:
                        if s.aligned_sequence in known_taxonomy:
                            known_sequences.append(s)
                        else:
                            unknown_sequences.append(s)
                    logging.debug("For sample %s, spkg %s, found %i known and %i unknown OTU sequences" %(
                        sample_name,
                        singlem_package.base_directory(),
                        len(known_sequences),
                        len(unknown_sequences)))
                            
                    if len(unknown_sequences) > 0:
                        tmp = tempfile.NamedTemporaryFile(
                            prefix='singlem.%s.' % sample_name,
                            suffix='.fasta')
                        SequenceExtractor().extract([s.name for s in unknown_sequences],
                                                    nucleotide_sequence_fasta_file,
                                                    tmp.name)
                    else:
                        tmp = None
                    extracted_reads.add(sample_name, singlem_package, tmp, known_sequences, unknown_sequences)

        return extracted_reads

    def _align_proteins_to_hmm(self, proteins_file, hmm_file):
        '''hmmalign proteins to hmm, and return an alignment object'''

        with tempfile.NamedTemporaryFile(prefix="singlem", suffix=".fasta") as f:
            cmd = "hmmalign %s %s |seqmagick convert --input-format stockholm - %s" % (hmm_file,
                                                              proteins_file,
                                                              f.name)
            extern.run(cmd)
            return SeqReader().alignment_from_alignment_file(f.name)

    def _seqs_to_counts_and_taxonomy(self, sequences, taxonomies, use_first_taxonomy=False,
                                     taxonomy_is_known_taxonomy=False):
        '''given an array of UnalignedAlignedNucleotideSequence objects, and hash of
        taxonomy file, yield over 'Info' objects that contain e.g. the counts of
        the aggregated sequences and corresponding median taxonomies.
        
        Parameters
        ----------
        use_first_taxonomy: boolean
            False: get a median taxonomy. True: use the taxonomy of the first encountered sequence

        '''
        class CollectedInfo:
            def __init__(self):
                self.count = 0
                self.taxonomies = []
                self.names = []
                self.coverage = 0.0
                self.aligned_lengths = []

        seq_to_collected_info = {}
        for s in sequences:
            try:
                if taxonomy_is_known_taxonomy:
                    tax = taxonomies[s.aligned_sequence].taxonomy
                else:
                    tax = taxonomies[s.name]
            except KeyError:
                # happens sometimes when HMMER picks up something where
                # diamond does not, or when --no_assign_taxonomy is specified.
                logging.debug("Did not find any taxonomy information for %s" % s.name)
                tax = ''
            try:
                collected_info = seq_to_collected_info[s.aligned_sequence]
            except KeyError:
                collected_info = CollectedInfo()
                seq_to_collected_info[s.aligned_sequence] = collected_info

            collected_info.count += 1
            if taxonomies: collected_info.taxonomies.append(tax)
            collected_info.names.append(s.name)
            collected_info.coverage += s.coverage_increment()
            collected_info.aligned_lengths.append(s.aligned_length)

        class Info:
            def __init__(self, seq, count, taxonomy, names, coverage, aligned_lengths):
                self.seq = seq
                self.count = count
                self.taxonomy = taxonomy
                self.names = names
                self.coverage = coverage
                self.aligned_lengths = aligned_lengths

        for seq, collected_info in seq_to_collected_info.iteritems():
            if use_first_taxonomy:
                tax = collected_info.taxonomies[0]
                if tax is None: tax = ''
            else:
                tax = self._median_taxonomy(collected_info.taxonomies)
            yield Info(seq,
                       collected_info.count,
                       tax,
                       collected_info.names,
                       collected_info.coverage,
                       collected_info.aligned_lengths)



    def _median_taxonomy(self, taxonomies):
        levels_to_counts = []
        for tax_string in taxonomies:
            for i, tax in enumerate(tax_string.split(';')):
                tax = tax.strip()
                if i >= len(levels_to_counts):
                    levels_to_counts.append({})
                try:
                    levels_to_counts[i][tax] += 1
                except KeyError:
                    levels_to_counts[i][tax] = 1


        median_tax = []
        for level_counts in levels_to_counts:
            max_count = 0
            max_tax = None
            for tax, count in level_counts.iteritems():
                if count > max_count:
                    max_count = count
                    max_tax = tax
            if float(max_count) / len(taxonomies) >= 0.5:
                median_tax.append(max_tax)
            else:
                break
        return '; '.join(median_tax)
    
    def _write_jplace_from_infos(self, input_jplace_io, infos, output_jplace_io):
        
        jplace = json.load(input_jplace_io)
        if jplace['version'] != 3:
            raise Exception("SingleM currently only works with jplace version 3 files, sorry")
        
        name_to_info = {}
        for info in infos:
            for name in info.names:
                name_to_info[name] = info
                
        # rewrite placements to be OTU-wise instead of sequence-wise
        orfm_utils = OrfMUtils()
        another_regex = re.compile(u'_\d+$')
        sequence_to_count = {}
        sequence_to_example_p = {}
        for placement in jplace['placements']:
            if 'nm' not in placement:
                raise Exception("Unexpected jplace format detected in placement %s" % placement)
            for name_and_count in placement['nm']:
                if len(name_and_count) != 2:
                    raise Exception("Unexpected jplace format detected in nm %s" % name_and_count)
                name, count = name_and_count
                real_name = another_regex.sub('', orfm_utils.un_orfm_name(name))
                info = name_to_info[real_name]
                sequence = info.seq
                
                try:
                    sequence_to_count[sequence] += count
                except KeyError:
                    sequence_to_count[sequence] = count
                    
                if real_name == info.names[0]:
                    sequence_to_example_p[sequence] = placement['p']
            
        new_placements = {}
        for sequence, example_p in sequence_to_example_p.items():
            new_placements[sequence] = {}
            new_placements[sequence]['nm'] = [[sequence, sequence_to_count[sequence]]]
            new_placements[sequence]['p'] = example_p
            
        jplace['placements'] = new_placements
        json.dump(jplace, output_jplace_io)

    def _graftm_command_prefix(self, is_protein):
        # --min_orf_length is unused for nucleotide HMMs but does no harm.
        cmd = "graftM graft "\
              "--verbosity %s "\
              "--input_sequence_type nucleotide " % self._graftm_verbosity
        if self._evalue: cmd += ' --evalue %s' % self._evalue
        if self._restrict_read_length: cmd += ' --restrict_read_length %i' % self._restrict_read_length
        
        if is_protein:
            cmd += " --min_orf_length %s " % self._min_orf_length
            if self._filter_minimum_protein:
                cmd += "  --filter_minimum %i" % self._filter_minimum_protein
        elif self._filter_minimum_nucleotide:
            cmd += "  --filter_minimum %i" % self._filter_minimum_nucleotide
            
        return cmd+' '

    def _search(self, singlem_package_database, forward_read_files):
        '''Find all reads that match one or more of the search HMMs in the
        singlem_package_database.

        Parameters
        ----------
        singlem_package_database: HmmDatabase
            packages to search the reads for
        forward_read_files: list of str
            paths to the sequences to be searched

        Returns
        -------
        SingleMPipeSearchResult
        '''
        logging.info("Using as input %i different sequence files e.g. %s" % (
            len(forward_read_files), forward_read_files[0]))
        graftm_protein_search_directory = os.path.join(
            self._working_directory, 'graftm_protein_search')
        graftm_nucleotide_search_directory = os.path.join(
            self._working_directory,'graftm_nucleotide_search')

        def run(hmm_paths, output_directory, is_protein):
            cmd = self._graftm_command_prefix(is_protein) + \
                  "--threads %i "\
                  "--forward %s "\
                  "--search_only "\
                  "--search_hmm_files %s "\
                  "--output_directory %s "\
                  "--aln_hmm_file %s " % (
                      self._num_threads,
                      ' '.join(forward_read_files),
                      ' '.join(hmm_paths),
                      output_directory,
                      hmm_paths[0])
            extern.run(cmd)
        
        # Run searches for proteins
        hmms = singlem_package_database.protein_search_hmm_paths()
        doing_proteins = False
        if len(hmms) > 0:
            doing_proteins = True
            logging.info("Searching for reads matching %i different protein HMMs" % len(hmms))
            run(hmms, graftm_protein_search_directory, True)

        # Run searches for nucleotides
        hmms = singlem_package_database.nucleotide_search_hmm_paths()
        doing_nucs = False
        if len(hmms) > 0:
            doing_nucs = True
            logging.info("Searching for reads matching %i different nucleotide HMMs" % len(hmms))
            run(hmms, graftm_nucleotide_search_directory, False)
            
        logging.info("Finished search phase")
        protein_graftm = GraftMResult(graftm_protein_search_directory, hmms) if \
                         doing_proteins else None
        nuc_graftm = GraftMResult(graftm_nucleotide_search_directory, hmms) if \
                     doing_nucs else None
        return SingleMPipeSearchResult(protein_graftm, nuc_graftm)

    def _align(self, search_result):
        graftm_separate_directory_base = os.path.join(self._working_directory, 'graftm_separates')
        os.mkdir(graftm_separate_directory_base)
        logging.info("Running separate alignments in GraftM..")
        commands = []

        def command(singlem_package, hit_files, is_protein):
            return self._graftm_command_prefix(is_protein) + \
                "--threads %i "\
                "--forward %s "\
                "--graftm_package %s --output_directory %s/%s "\
                "--search_only" % (
                    1, #use 1 thread since most likely better to parallelise processes with extern
                    ' '.join(hit_files),
                    singlem_package.graftm_package_path(),
                    graftm_separate_directory_base,
                    os.path.basename(singlem_package.graftm_package_path()))

        # Gather commands for aligning protein packages
        for singlem_package in self._singlem_package_database.protein_packages():
            commands.append(command(singlem_package, search_result.protein_hit_paths().values(), True))
        # Gather commands for aligning nucleotide packages.
        for singlem_package in self._singlem_package_database.nucleotide_packages():
            temporary_hit_files = [tf for _, tf in \
                search_result.direction_corrected_nucleotide_read_files()]
            commands.append(command(singlem_package, temporary_hit_files, False))

        extern.run_many(commands, num_threads=self._num_threads)
        return SingleMPipeAlignSearchResult(
            graftm_separate_directory_base, search_result.samples_with_hits())

    def _assign_taxonomy(self, extracted_reads, assignment_method):
        graftm_align_directory_base = os.path.join(self._working_directory, 'graftm_aligns')
        os.mkdir(graftm_align_directory_base)
        commands = []
        for singlem_package, sample_names, tmp_grafts in extracted_reads.each_package_wise():
            tmpnames = list([tg.name for tg in tmp_grafts if tg])
            if len(tmpnames) > 0:
                cmd = "%s "\
                      "--threads %i "\
                      "--forward %s "\
                      "--graftm_package %s "\
                      "--output_directory %s/%s "\
                      "--max_samples_for_krona 0 "\
                      "--assignment_method %s" % (
                          self._graftm_command_prefix(singlem_package.is_protein_package()),
                          self._num_threads,
                          ' '.join(tmpnames),
                          singlem_package.graftm_package_path(),
                          graftm_align_directory_base,
                          singlem_package.graftm_package_basename(),
                          assignment_method)
                commands.append(cmd)
        extern.run_many(commands, num_threads=1)
        logging.info("Finished running taxonomic assignment with graftm")
        return SingleMPipeTaxonomicAssignmentResult(graftm_align_directory_base)

class SingleMPipeSearchResult:
    def __init__(self, graftm_protein_result, graftm_nucleotide_result):
        self._protein_result = graftm_protein_result
        self._nucleotide_result = graftm_nucleotide_result

    def protein_hit_paths(self):
        '''Return a dict of sample name to corresponding '_hits.fa' files generated in
        the search step. Do not return those samples where there were no hits.

        '''
        if self._protein_result is None:
            return {}
        else:
            return self._protein_result.unaligned_sequence_paths(require_hits=True)

    def direction_corrected_nucleotide_read_files(self):
        '''For nucleotide HMMs: Iterate over the sample names plus a fasta filename per
        sample, fasta files that are 'direction-corrected' i.e. contain
        sequences in the direction that they were aligned. These tempfiles must
        be closed by code using this function. Do not use this method for
        protein HMMs.

        '''
        if self._nucleotide_result is None:
            # No nucleotide singlem packages
            return
        for sample_name in self._nucleotide_result.sample_names(require_hits=True):
            forward_reads = set()
            reverse_reads = set()
            for hmmout in self._nucleotide_result.hmmout_paths_from_sample_name(sample_name):
                hmmout_result = HMMSearchResult.import_from_nhmmer_table(hmmout)
                for hit in hmmout_result.each(
                        [SequenceSearchResult.QUERY_ID_FIELD,
                         SequenceSearchResult.ALIGNMENT_DIRECTION]):
                    if hit[1]:
                        forward_reads.add(hit[0])
                    else:
                        reverse_reads.add(hit[0])
            nucs = self._nucleotide_result.unaligned_sequences_path_from_sample_name(sample_name)
            
            yieldme = os.path.join(self._nucleotide_result.output_directory,
                                   "%s_hits.fa" % sample_name)
            SequenceExtractor().extract_forward_and_reverse_complement(
                forward_reads, reverse_reads, nucs, yieldme)
            if os.stat(yieldme).st_size > 0:
                yield sample_name, yieldme

    def samples_with_hits(self):
        '''Return a list of sample names that had at least one hit'''
        return list(set(itertools.chain(
            self.protein_hit_paths().keys(),
            self._nucleotide_result.sample_names(require_hits=True) if \
                self._nucleotide_result else [])))

class SingleMPipeAlignSearchResult:
    def __init__(self, graftm_separate_directory_base, sample_names):
        self._graftm_separate_directory_base = graftm_separate_directory_base
        self._sample_names = sample_names

    def _base_dir(self, sample_name, singlem_package):
        return os.path.join(
            self._graftm_separate_directory_base,
            os.path.basename(singlem_package.graftm_package_path()),
            '%s_hits' % sample_name)

    def prealigned_sequence_files(self, sample_name, singlem_package):
        '''Yield a path to the sequences that were aligned to the HMM. In the nucleotide
        case, this is a temporary file(name) which is a concatenation of forward
        and reverse.

        '''
        if singlem_package.is_protein_package():
            yield os.path.join(
                self._base_dir(sample_name, singlem_package),
                "%s_hits_orf.fa" % sample_name)
        else:
            yield self.nucleotide_sequence_file(sample_name, singlem_package)

    def nucleotide_sequence_file(self, sample_name, singlem_package):
        return os.path.join(
            self._base_dir(sample_name, singlem_package),
            "%s_hits_hits.fa" % sample_name)

    def sample_names(self):
        return self._sample_names

class SingleMPipeTaxonomicAssignmentResult:
    def __init__(self, graftm_output_directory):
        self._graftm_output_directory = graftm_output_directory

    def _base_dir(self, sample_name, singlem_package, tmpbase):
        return os.path.join(self._graftm_output_directory,
                            '%s/%s' % (
                                singlem_package.graftm_package_basename(),
                                re.sub('\.fasta$','',tmpbase)))

    def protein_orf_file(self, sample_name, singlem_package, tmpbase):
        return os.path.join(self._base_dir(sample_name, singlem_package, tmpbase),
                            "%s_orf.fa" % tmpbase)
    
    def prealigned_sequence_file(self, sample_name, singlem_package, tmpbase):
        '''path to the sequences that were aligned (ORF for proteins, regular seqs for
        nucleotide).

        '''
        if singlem_package.is_protein_package():
            return self.protein_orf_file(sample_name, singlem_package, tmpbase)
        else:
            return self.nucleotide_hits_file(sample_name, singlem_package, tmpbase)
    
    def nucleotide_hits_file(self, sample_name, singlem_package, tmpbase):
        return os.path.join(self._base_dir(sample_name, singlem_package, tmpbase),
                            "%s_hits.fa" % tmpbase)

    def diamond_assignment_file(self, sample_name, singlem_package, tmpbase):
        return os.path.join(self._base_dir(sample_name, singlem_package, tmpbase),
                            '%s_diamond_assignment.daa' % tmpbase)

    def read_tax_file(self, sample_name, singlem_package, tmpbase):
        return os.path.join(self._base_dir(sample_name, singlem_package, tmpbase),
                            '%s_read_tax.tsv' % tmpbase)

class ExtractedReads:
    def __init__(self, singlem_packages):
        self._singlem_packages = singlem_packages
        self._sample_to_array = {}

    def add(self, sample_name, singlem_package, tmp_file, known_sequences, unknown_sequences):
        arr = [singlem_package, tmp_file, known_sequences, unknown_sequences]
        if sample_name in self._sample_to_array:
            self._sample_to_array[sample_name].append(arr)
        else:
            self._sample_to_array[sample_name] = [arr]

    def __iter__(self):
        '''yield sample, single_package, tmp_file, known_sequences, unknown_sequences'''
        for sample, arrs in self._sample_to_array.items():
            for arr in arrs:
                yield sample, arr[0], arr[1], arr[2], arr[3]

    def each_package_wise(self):
        '''yield once per pkg: [singlem_package, list of samples, list of tmp_files]'''
        pkg_basename_to_pkg_samples_tmp_files = {}
        for sample, arrs in self._sample_to_array.items():
            for arr in arrs:
                pkg_base = arr[0].base_directory()
                if pkg_base in pkg_basename_to_pkg_samples_tmp_files:
                    pkg_basename_to_pkg_samples_tmp_files[pkg_base][1].append(sample)
                    pkg_basename_to_pkg_samples_tmp_files[pkg_base][2].append(arr[1])
                else:
                    pkg_basename_to_pkg_samples_tmp_files[pkg_base] = \
                        [arr[0], [sample], [arr[1]]]
        for arr in pkg_basename_to_pkg_samples_tmp_files.values():
            yield arr

    def each_sample(self):
        for sample in self._sample_to_array.keys():
            yield sample

    def cleanup(self):
        # remove these tempfiles because otherwise errors are spewed
        # when they are cleaned up after the tempdir is gone
        for sample, arrs in self._sample_to_array.items():
            for arr in arrs:
                if arr[1] is not None: arr[1].close()
