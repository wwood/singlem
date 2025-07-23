import os
import logging
import extern
from Bio import SeqIO
from io import StringIO
import multiprocessing
import tempfile

from graftm.hmmsearcher import HmmSearcher

from .sequence_classes import SeqReader, AlignedProteinSequence, Sequence
from .metagenome_otu_finder import MetagenomeOtuFinder
from . import sequence_extractor as singlem_sequence_extractor
from .streaming_hmm_search_result import StreamingHMMSearchResult
from .utils import OrfMUtils


# Must be defined outside a class so that it is pickle-able, so multiprocessing can work
def _run_individual_extraction(sample_name, singlem_package, sequence_files_for_alignment, separate_search_result, include_inserts, known_taxonomy):
    if singlem_package.is_protein_package():
        nucleotide_sequence_fasta_file = \
            separate_search_result.nucleotide_sequence_file(
                sample_name, singlem_package)
    else:
        nucleotide_sequence_fasta_file = sequence_files_for_alignment

    if separate_search_result.analysing_pairs:
        logging.debug("Extracting forward reads")
        if os.path.exists(sequence_files_for_alignment[0]):
            prots = SeqReader().readfq(open(sequence_files_for_alignment[0]))
        else:
            prots = []
        readset1 = _extract_reads(
            singlem_package,
            sample_name,
            prots,
            nucleotide_sequence_fasta_file[0],
            include_inserts,
            known_taxonomy,
            'forward')
        logging.debug("Extracting reverse reads")
        if os.path.exists(sequence_files_for_alignment[1]):
            prots = SeqReader().readfq(open(sequence_files_for_alignment[1]))
        else:
            prots = []
        readset2 = _extract_reads(
            singlem_package,
            sample_name,
            prots,
            nucleotide_sequence_fasta_file[1],
            include_inserts,
            known_taxonomy,
            'reverse')
        return [readset1, readset2]
    else:
        if os.path.exists(sequence_files_for_alignment):
            prots = SeqReader().readfq(open(sequence_files_for_alignment))
        else:
            prots = []
        readset = _extract_reads(
            singlem_package,
            sample_name,
            prots,
            nucleotide_sequence_fasta_file,
            include_inserts,
            known_taxonomy,
            'forward')
        return readset

def _extract_reads(
        singlem_package,
        sample_name,
        prealigned_protein_sequences,
        nucleotide_sequence_fasta_file,
        include_inserts,
        known_taxonomy,
        read_direction):

    aligned_seqs = _get_windowed_sequences(
        prealigned_protein_sequences,
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
    logging.debug("For sample {} ({}), spkg {}, found {} known and {} unknown OTU sequences".format(
        sample_name,
        read_direction,
        singlem_package.base_directory(),
        len(known_sequences),
        len(unknown_sequences)))

    if len(unknown_sequences) > 0:
        extractor = singlem_sequence_extractor.SequenceExtractor()
        logging.debug("Extracting reads for {} from {} ..".format(
            singlem_package.base_directory(), nucleotide_sequence_fasta_file
        ))
        seqs = extractor.extract_and_read(
            [s.name for s in unknown_sequences],
            nucleotide_sequence_fasta_file)
    else:
        seqs = []
    readset = ExtractedReadSet(
        sample_name, singlem_package,
        seqs, known_sequences, unknown_sequences)
    return readset

def _get_windowed_sequences(protein_sequences, nucleotide_sequence_file,
                            singlem_package, include_inserts):
    if not os.path.exists(nucleotide_sequence_file) or \
        os.stat(nucleotide_sequence_file).st_size == 0: return []
    nucleotide_sequences = SeqReader().read_nucleotide_sequences(nucleotide_sequence_file)
    protein_alignment = _align_proteins_to_hmm(
        protein_sequences,
        singlem_package.graftm_package().alignment_hmm_path())
    return MetagenomeOtuFinder().find_windowed_sequences(
        protein_alignment,
        nucleotide_sequences,
        singlem_package.window_size(),
        include_inserts,
        singlem_package.is_protein_package(),
        best_position=singlem_package.singlem_position())

def _align_proteins_to_hmm(protein_sequences, hmm_file):
    '''hmmalign proteins to hmm, and return an alignment object

    Parameters
    ----------
    protein_sequences: generator / list of tuple(name,sequence) objects
    from SeqReader().

    '''
    cmd = "hmmalign '{}' /dev/stdin".format(hmm_file)
    logging.debug("Running command: {}".format(cmd))
    output = extern.run(cmd, stdin=''.join([
        ">{}\n{}\n".format(s[0], s[1]) for s in protein_sequences]))
    logging.debug("Finished command: {}".format(cmd))
    protein_alignment = []
    for record in SeqIO.parse(StringIO(output), 'stockholm'):
        protein_alignment.append(AlignedProteinSequence(record.name, str(record.seq)))
    if len(protein_alignment) > 0:
        logging.debug("Read in %i aligned sequences e.g. %s %s" % (
            len(protein_alignment),
            protein_alignment[0].name,
            protein_alignment[0].seq))
    else:
        logging.debug("No aligned sequences found for this HMM")
    return protein_alignment


def _filter_sequences_through_hmmsearch(
    singlem_package,
    prefilter_result,
    min_orf_length,
    translation_table,
    evalue):
    """Return a generator of Sequence objects that match one or more of the
    search HMMs in the graftm_package. To save RAM, the full set of sequences is
    never read in.
    """
    
    graftm_package = singlem_package.graftm_package()
    target_sequence_ids = set(singlem_package.get_sequence_ids())

    # Can re-use HmmSearcher, useful for when there is >1 search HMM
    searcher = HmmSearcher(len(graftm_package.search_hmm_paths()), '--domE {}'.format(evalue))
    
    output_tempfiles = list([
        tempfile.NamedTemporaryFile(prefix='singlem_hmmsearch') for _ in graftm_package.search_hmm_paths()
    ])
    with tempfile.NamedTemporaryFile(prefix='singlem_hmmsearch_input') as input_tf:
        count, example = _generate_package_specific_fasta_input(target_sequence_ids, prefilter_result, input_tf)
        logging.debug("Running {} sequences through HMMSEARCH e.g. {}".format(count, example))
        input_tf.flush()
        if count == 0:
            return

        # With some hoop jumping it should be possible to stream this, but eh
        # for now.
        searcher.hmmsearch(
            'orfm -c {} -m {} {}'.format(translation_table, min_orf_length, input_tf.name),
            graftm_package.search_hmm_paths(),
            list([tf.name for tf in output_tempfiles]))

    # Stream reading of the hmmout file. Only need the query ID?
    seqs_to_extract = set()
    for output_tempfile in output_tempfiles:
        for orfm_seq_id in StreamingHMMSearchResult.yield_from_hmmsearch_table(output_tempfile.name):
            seqs_to_extract.add(OrfMUtils().un_orfm_name(orfm_seq_id))

    for (qseqid, seq) in _yield_target_sequences(target_sequence_ids, prefilter_result):
        if qseqid in seqs_to_extract:
            yield Sequence(qseqid, seq)

def _yield_target_sequences(target_sequence_ids, prefilter_result):
    with open(prefilter_result.query_sequences_file) as f:
        for (qseqid, seq, _) in SeqReader().readfq(f):
            if prefilter_result.best_hits[qseqid] in target_sequence_ids:
                yield (qseqid, seq)

def _generate_package_specific_fasta_input(
    target_sequence_ids, prefilter_result, output_io):
    """Write sequences that are in the prefilter result that match the singlem
    package to output_io.
    """
    
    count = 0
    example = None


    for (qseqid, seq) in _yield_target_sequences(target_sequence_ids, prefilter_result):
        count += 1
        if example is not None:
            example = qseqid
        output_io.write(">{}\n{}\n".format(qseqid, seq).encode())

    return count, example


# Function to chunk sequences by max bp
# there is an edge case where a sequence is longer than the max bp
# in this case, the sequence will be in a chunk by itself
def chunk_sequences_by_bp(sequences, soft_max_bp):
    chunks = []
    current_chunk = []
    current_bp = 0

    for seq in sequences:
        seq_length = len(seq.seq)
        
        if current_bp + seq_length > soft_max_bp:
            chunks.append(current_chunk)
            current_chunk = []
            current_bp = 0

        current_chunk.append(seq)
        current_bp += seq_length

    if current_chunk:
        chunks.append(current_chunk)

    return chunks


def _extract_reads_by_diamond_for_package_and_sample(prefilter_result, spkg,
    sample_name, min_orf_length, include_inserts, translation_table, hmmsearch_evalue):

    # Happens when there is no hits
    if len(prefilter_result.best_hits) == 0:
        # Add something so that when analysing pairs and one side has no hits, indexing errors don't happen.
        return ExtractedReadSet(
            sample_name, spkg,
            [], [], []
        )

    # Sequences have to be run through hmmsearch, because sometimes DIAMOND
    # returns some less than good quality hits, and taking those hits directly
    # to hmmalign causes odd sequences to be counted in.

    sequences = list(_filter_sequences_through_hmmsearch(
        spkg,
        prefilter_result,
        min_orf_length,
        translation_table,
        hmmsearch_evalue))

    # On very rare occasions, the same sequence can be returned twice, causing
    # hmmalign to give unexpected format output. So dedup. For instance using
    # S3.2.1 and
    # /home/woodcrob/m/abisko/data/flat20230929_all/201607_SubstrateInc_9to19.G_F_S_T50.1.fq.gz
    # &
    # /home/woodcrob/m/abisko/data/flat20230929_all/201607_SubstrateInc_9to19.G_F_S_T50.2.fq.gz
    # causes the eventual bug.

    sequences = list(set(sequences))

    # limit the number of bp per chunk to avoid memory issues
    soft_max_bp_per_chunk = 1000000
    window_seqs = []
    
    # to help debug when multi-processing
    pid = os.getpid()
    logging.debug("PID: {}".format(pid))

    # Chunk sequences by max bp
    sequence_chunks = chunk_sequences_by_bp(sequences, soft_max_bp_per_chunk)

    logging.debug("[PID: {}] Processing {} chunks".format(pid, len(sequence_chunks)))
    for chunk_sequences in sequence_chunks:
        # create temporary files to store input seqs and output alignment
        with tempfile.NamedTemporaryFile(prefix='singlem_hmmalign') as output_tempfile, \
             tempfile.NamedTemporaryFile(prefix='singlem_input') as input_tempfile:
            
            # write sequences to input file
            input_tempfile.write(">dummy\n{}\n".format('A'*min_orf_length).encode())
            for s in chunk_sequences:
                input_tempfile.write(">{}\n{}\n".format(s.name, s.seq).encode())
            input_tempfile.flush()      
            
            # define command
            cmd = "orfm -c {} -m {} {} | hmmalign -o '{}' '{}' /dev/stdin".format( 
                translation_table, 
                min_orf_length, 
                input_tempfile.name,
                output_tempfile.name,
                spkg.graftm_package().alignment_hmm_path())
            
            # run command
            extern.run(cmd)

            # Convert to AlignedProteinSequence
            protein_alignment = []
            for record in SeqIO.parse(output_tempfile.name, 'stockholm'):
                protein_alignment.append(AlignedProteinSequence(record.name, str(record.seq)))

        if len(protein_alignment) > 0:
            logging.debug("[PID: %s] Read in %i aligned sequences from this chunk e.g. %s %s" % (
                pid,
                len(protein_alignment),
                protein_alignment[0].name,
                protein_alignment[0].seq))
        else:
            logging.debug("[PID: {}] No aligned sequences found for this HMM".format(pid))

        # Extract OTU sequences
        nucleotide_sequence_hash = {s.name: s.seq for s in sequences}
        logging.debug("[PID: {}] First sequence: {} / {}".format(pid, protein_alignment[0].name, protein_alignment[0].seq))
        # Window sequences must be found for each chunk, otherwise the
        # alignments won't line up re insert characters, between chunks.
        window_seqs.extend(MetagenomeOtuFinder().find_windowed_sequences(
            protein_alignment,
            nucleotide_sequence_hash,
            spkg.window_size(),
            include_inserts,
            spkg.is_protein_package(), # Always true
            best_position=spkg.singlem_position()))


    logging.debug("[PID: {}] Found {} window sequences for spkg {}".format(pid, len(window_seqs),spkg.base_directory()))
    return ExtractedReadSet(
        sample_name, spkg,
        sequences, [], window_seqs
    )


class PipeSequenceExtractor:
    '''Part of the singlem pipe, abstracted out here so that it can be run in
    parallel
    '''

    def extract_relevant_reads_from_separate_search_result(self, singlem_package_database, num_threads, separate_search_result, include_inserts, known_taxonomy):
        '''Given a SingleMPipeSeparateSearchResult, extract reads that will be used as
        part of the singlem choppage process.

        Returns
        -------
        ExtractedReads object
        '''

        extracted_reads = ExtractedReads(separate_search_result.analysing_pairs)

        # Collect all possibilities
        to_iterate = []
        for sample_name in separate_search_result.sample_names():
            for singlem_package in singlem_package_database:
                for sequence_files_for_alignment in separate_search_result.sequence_files_for_alignment(
                        sample_name, singlem_package):
                    to_iterate.append((sample_name, singlem_package, sequence_files_for_alignment, separate_search_result, include_inserts, known_taxonomy))

        # Multiprocess across all instances if num_threads > 1
        if num_threads > 1:
            pool = multiprocessing.Pool(num_threads)
            extraction_processes = [pool.apply_async(_run_individual_extraction, args=myargs) for myargs in to_iterate]
            for readset_possibly_paired_process in extraction_processes:
                extracted_reads.add(readset_possibly_paired_process.get())
            pool.close()
            pool.join()
        else:
            for myargs in to_iterate:
                extracted_reads.add(_run_individual_extraction(*myargs))

        return extracted_reads


    def extract_relevant_reads_from_diamond_prefilter(self,
        num_threads, singlem_package_database,
        diamond_forward_search_results, diamond_reverse_search_results,
        analysing_pairs, include_inserts, min_orf_length, translation_table,
        hmmsearch_evalue):
        '''Return an ExtractedReads object built by running hmmalign on
        sequences, aligning to each HMM only sequences that
        have a best hit to sequences from that singlem package

        Returns
        -------
        ExtractedReads object
        '''
        extracted_reads = ExtractedReads(analysing_pairs)

        if num_threads > 1:
            # Multiprocessing incurs a RAM overhead, so only do it if required.
            pool = multiprocessing.Pool(num_threads)
        else:
            pool = None

        logging.debug("Aligning and extracting forward reads ..")
        logging.debug("Extracting reads from {} forward search results".format(len(diamond_forward_search_results)))
        forward_extraction_process_lists_per_sample = []
        for diamond_search_result in diamond_forward_search_results:
            extraction_processes = self._extract_relevant_reads_from_diamond_prefilter_from_one_search_result(
                pool, singlem_package_database, diamond_search_result, include_inserts, min_orf_length, translation_table, hmmsearch_evalue
            )
            forward_extraction_process_lists_per_sample.append(extraction_processes)

        logging.debug("Aligning and extracting reverse reads ..")
        reverse_extraction_process_lists_per_sample = []
        if analysing_pairs:
            for diamond_search_result in diamond_reverse_search_results:
                extraction_processes = self._extract_relevant_reads_from_diamond_prefilter_from_one_search_result(
                    pool, singlem_package_database, diamond_search_result, include_inserts, min_orf_length, translation_table, hmmsearch_evalue
                )
                reverse_extraction_process_lists_per_sample.append(extraction_processes)

        # Could pickle the extracted reads an instead return a list of the files 
        # to read from, or write to a fasta.
        # This could improve memory even more, but might slow it down.
        if analysing_pairs:
            for (fwds, revs) in zip(forward_extraction_process_lists_per_sample,reverse_extraction_process_lists_per_sample):
                for (fwd, rev) in zip(fwds, revs):
                    if pool is None:
                        extracted_reads.add((fwd, rev))
                    else:
                        extracted_reads.add((fwd.get(), rev.get()))
        else:
            for fwds in forward_extraction_process_lists_per_sample:
                for fwd in fwds:
                    if pool is None:
                        extracted_reads.add(fwd)
                    else:
                        extracted_reads.add(fwd.get())
        if pool is not None:
            pool.close()
            pool.join()
        logging.debug("Finished aligning and extracting reads")

        return extracted_reads

    def _extract_relevant_reads_from_diamond_prefilter_from_one_search_result(
        self, 
        pool, 
        singlem_package_database, 
        diamond_search_result, 
        include_inserts, 
        min_orf_length, 
        translation_table, 
        hmmsearch_evalue):
        '''Return results for one search result (sample). If pool is a
        multiprocessing pool, run individual extractions with that, otherwise it
        should be None, which means just run each extraction serially without
        invoking multiprocessing at all.

        Returns array of ExtractedReadSet objects or multiprocessing processes
        that yield those objects.
        '''

        # Determine sample name. In order to have compatible sample names with
        # the hmmsearch mode, remove filename suffixes.
        sample_name = diamond_search_result.sample_name()

        # TODO: I think there are some efficiency gains to be made here as many 
        # of the processes started here have no reads and just exit.
        # This requires starting up processes for no reason.
        # Would just need to check if there are any reads for the singlem 
        # package before starting the process.


        # Align each read via hmmsearch and pick windowed sequences
        extraction_of_read_set_processes = []
        for spkg in singlem_package_database:
            if pool is None:
                extraction_of_read_set_processes.append(
                    _extract_reads_by_diamond_for_package_and_sample(
                        diamond_search_result, spkg, sample_name, min_orf_length, include_inserts, translation_table, hmmsearch_evalue))
            else:
                extraction_of_read_set_processes.append(
                    pool.apply_async(
                        _extract_reads_by_diamond_for_package_and_sample, args=(
                    diamond_search_result, spkg, sample_name, min_orf_length, include_inserts, translation_table, hmmsearch_evalue)))

        return extraction_of_read_set_processes


class ExtractedReads:
    '''Collection class for ExtractedReadSet objects'''
    def __init__(self, analysing_pairs):
        self._sample_to_extracted_read_objects = {}
        self.analysing_pairs = analysing_pairs

    def add(self, extracted_read_set):
        '''Add an ExtractedReadSet, or if analysing_pairs, a pair of them'''
        if self.analysing_pairs:
            sample_name = extracted_read_set[0].sample_name
        else:
            sample_name = extracted_read_set.sample_name

        if sample_name in self._sample_to_extracted_read_objects:
            self._sample_to_extracted_read_objects[sample_name].append(extracted_read_set)
        else:
            self._sample_to_extracted_read_objects[sample_name] = [extracted_read_set]

    def __iter__(self):
        '''yield sample, extracted_read_sets'''
        for readsets in self._sample_to_extracted_read_objects.values():
            for readset in readsets:
                yield readset

    def empty(self):
        '''True if all readsets have no sequences, else False'''
        for readsets in self._sample_to_extracted_read_objects.values():
            for readset in readsets:
                if self.analysing_pairs:
                    if len(readset[0].unknown_sequences) > 0  or len(readset[0].known_sequences) > 0 or \
                       len(readset[1].unknown_sequences) > 0  or len(readset[1].known_sequences) > 0:
                        return False
                else:
                    if len(readset.unknown_sequences) > 0  or len(readset.known_sequences) > 0:
                        return False
        return True

    def each_package_wise(self):
        '''yield once per pkg: [singlem_package, ExtractedReadSet objects with all
        samples / sequences that have extracted sequences from it]

        '''
        pkg_basename_to_readsets = {}
        for sample, readsets in self._sample_to_extracted_read_objects.items():
            for readset in readsets:
                if self.analysing_pairs:
                    pkg_base = readset[0].singlem_package.base_directory()
                else:
                    pkg_base = readset.singlem_package.base_directory()
                if pkg_base in pkg_basename_to_readsets:
                    pkg_basename_to_readsets[pkg_base].append(readset)
                else:
                    pkg_basename_to_readsets[pkg_base] = [readset]
        for readsets in pkg_basename_to_readsets.values():
            if self.analysing_pairs:
                yield readsets[0][0].singlem_package, readsets
            else:
                yield readsets[0].singlem_package, readsets

    def each_sample(self):
        return self._sample_to_array.keys()


class ExtractedReadSet:
    def __init__(self, sample_name, singlem_package, sequences,
                 known_sequences, unknown_sequences):
        '''
        * sequences: list of Sequence objects
        * unknown_sequences: list of UnalignedAlignedNucleotideSequence objects
        '''

        self.sample_name = sample_name
        self.singlem_package = singlem_package
        self.sequences = sequences
        self.known_sequences = known_sequences
        self.unknown_sequences = unknown_sequences
        self.tmpfile_basename = None # Used as part of pipe, making this object
                                     # not suitable for use outside that
                                     # setting.

    def remove_duplicate_sequences(self):
        '''When OrfM is run on a genome, then sometimes the same stretch of the
        genome will be in 2 transcripts, so is inappropriately included twice.
        Remove these cases.'''

        # Do not import at the top level because it isn't compatible with newest Python.
        import pyranges as pr
        
        logging.debug("Before duplicate removal for {}, have {} sequences".format(
            self.singlem_package.graftm_package_basename(), len(self.unknown_sequences)))
        
        # Order the unknown_sequences in order of decreasing length, as longer
        # ones are more likely to be correct.
        sorted_unknowns = sorted(self.unknown_sequences, key=lambda x: len(x.unaligned_sequence), reverse=True)

        unique_aligned_to_genome_ranges = {}
        unknown_to_return = []

        orfm_utils = OrfMUtils()

        # Add each
        for seq in sorted_unknowns:
            start, _, _ = orfm_utils.un_orfm_start_frame_number(seq.orf_name)
            grange = pr.PyRanges(
                chromosomes=[orfm_utils.un_orfm_name(seq.name)],
                starts=[start],
                ends=[start + len(seq.unaligned_sequence) - 1], # minus one since edges are inclusive
            )
            if seq.aligned_sequence not in unique_aligned_to_genome_ranges:
                unique_aligned_to_genome_ranges[seq.aligned_sequence] = grange
                unknown_to_return.append(seq)
            else:
                if len(unique_aligned_to_genome_ranges[seq.aligned_sequence].intersect(grange)) == 0:
                    unique_aligned_to_genome_ranges[seq.aligned_sequence] = \
                        unique_aligned_to_genome_ranges[seq.aligned_sequence].join(grange)
                    unknown_to_return.append(seq)
        
        # Reconstruct
        logging.debug("After duplicate removal, now have {} sequences".format(len(unknown_to_return)))
        self.unknown_sequences = unknown_to_return

