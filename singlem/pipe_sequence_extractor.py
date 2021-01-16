import os
import logging
import extern
from Bio import SeqIO
from io import StringIO
import multiprocessing

from .sequence_classes import SeqReader, AlignedProteinSequence
from .metagenome_otu_finder import MetagenomeOtuFinder
from . import sequence_extractor as singlem_sequence_extractor

# Must be defined outside a class so that it is pickle-able, so multiprocessing can work
def _run_individual_extraction(sample_name, singlem_package, prealigned_file, alignment_result, include_inserts, known_taxonomy):
    if singlem_package.is_protein_package():
        nucleotide_sequence_fasta_file = \
            alignment_result.nucleotide_sequence_file(
                sample_name, singlem_package)
    else:
        nucleotide_sequence_fasta_file = prealigned_file

    if alignment_result.analysing_pairs:
        logging.debug("Extracting forward reads")
        if os.path.exists(prealigned_file[0]):
            prots = SeqReader().readfq(open(prealigned_file[0]))
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
        if os.path.exists(prealigned_file[1]):
            prots = SeqReader().readfq(open(prealigned_file[1]))
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
        if os.path.exists(prealigned_file):
            prots = SeqReader().readfq(open(prealigned_file))
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


class PipeSequenceExtractor:
    '''Part of the singlem pipe, abstracted out here so that it can be run in parallel

    '''

    def extract_relevant_reads(self, singlem_package_database, num_threads, alignment_result, include_inserts, known_taxonomy):
        '''Given a SingleMPipeAlignSearchResult, extract reads that will be used as
        part of the singlem choppage process.

        Returns
        -------
        ExtractedReads object
        '''

        extracted_reads = ExtractedReads(alignment_result.analysing_pairs)

        # Collect all possibilities
        to_iterate = []
        for sample_name in alignment_result.sample_names():
            for singlem_package in singlem_package_database:
                for prealigned_file in alignment_result.prealigned_sequence_files(
                        sample_name, singlem_package):
                    to_iterate.append((sample_name, singlem_package, prealigned_file, alignment_result, include_inserts, known_taxonomy))

        # Multiprocess across all instances
        pool = multiprocessing.Pool(num_threads)
        extracted_readsets = [pool.apply(_run_individual_extraction, args=myargs) for myargs in to_iterate]
        for readset_possibly_paired in extracted_readsets:
            extracted_reads.add(readset_possibly_paired)

        return extracted_reads

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
        self.sample_name = sample_name
        self.singlem_package = singlem_package
        self.sequences = sequences
        self.known_sequences = known_sequences
        self.unknown_sequences = unknown_sequences
        self.tmpfile_basename = None # Used as part of pipe, making this object
                                     # not suitable for use outside that
                                     # setting.
