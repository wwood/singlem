from Bio.Seq import Seq
import logging
import re
from singlem import OrfMUtils

class Sequence:
    '''Simple name+sequence object'''
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq

class AlignedProteinSequence(Sequence):
    def un_orfm_name(self):
        return OrfMUtils().un_orfm_name(self.name)

    def orfm_nucleotides(self, nucleotide_sequence):
        m = re.search('_(\d+)_(\d+)_\d+$', self.name)
        start = int(m.groups(0)[0])-1
        translated_seq = nucleotide_sequence[start:(start+3*self.unaligned_length())]
        logging.debug("Returning orfm nucleotides %s" % translated_seq)
        if int(m.groups(0)[1]) > 3:
            # revcomp type frame
            return(str(Seq(translated_seq).reverse_complement()))
        else:
            return(translated_seq)

    def unaligned_length(self):
        return len(re.sub('-','',self.seq))

class UnalignedAlignedNucleotideSequence:
    '''Represent a nucleotide sequence (aligned in protein space or nucleotide
    space), together with the nucleotide sequence that it came from.

    '''

    def __init__(self, name, orf_name, aligned_sequence, unaligned_sequence, aligned_length):
        '''
        Parameters
        ---------
        name: str
            name of the sequence
        orf_name: str
            name of the ORF
        aligned_sequence: str
            aligned nucleotide sequence
        unaligned_sequence: str
            unaligned nucleotide sequence
        aligned_length:
            the number of nucleotides used in the alignment, including columns
            that were removed as not aligned
        '''
        self.name = name
        self.orf_name = orf_name
        self.aligned_sequence = aligned_sequence
        self.unaligned_sequence = unaligned_sequence
        self.aligned_length = aligned_length

    def coverage_increment(self):
        '''Given the alignment came from a read of length
        original_nucleotide_sequence_length, how much coverage does the
        observation of this aligned sequence indicate?'''
        return float(len(self.unaligned_sequence))/\
            (len(self.unaligned_sequence)-self.aligned_length+1)


class SeqReader:
    # Stolen from https://github.com/lh3/readfq/blob/master/readfq.py
    def readfq(self, fp): # this is a generator function
        last = None # this is a buffer keeping the last unprocessed line
        while True: # mimic closure; is it a bad idea?
            if not last: # the first record or a record following a fastq
                for l in fp: # search for the start of the next record
                    if l[0] in '>@': # fasta/q header line
                        last = l[:-1] # save this line
                        break
            if not last: break
            name, seqs, last = last[1:].partition(" ")[0], [], None
            for l in fp: # read the sequence
                if l[0] in '@+>':
                    last = l[:-1]
                    break
                seqs.append(l[:-1])
            if not last or last[0] != '+': # this is a fasta record
                yield name, ''.join(seqs), None # yield a fasta record
                if not last: break
            else: # this is a fastq record
                seq, leng, seqs = ''.join(seqs), 0, []
                for l in fp: # read the quality
                    seqs.append(l[:-1])
                    leng += len(l) - 1
                    if leng >= len(seq): # have read enough quality
                        last = None
                        yield name, seq, ''.join(seqs); # yield a fastq record
                        break
                if last: # reach EOF before reading enough quality
                    yield name, seq, None # yield a fasta record instead
                    break

    def read_nucleotide_sequences(self, nucleotide_file):
        nucleotide_sequences = {}
        for name, seq, _ in self.readfq(open(nucleotide_file)):
            nucleotide_sequences[name] = seq
        return nucleotide_sequences

    def alignment_from_alignment_file(self, alignment_file):
        protein_alignment = []
        for name, seq, _ in self.readfq(open(alignment_file)):
            protein_alignment.append(AlignedProteinSequence(name, seq))
        if len(protein_alignment) > 0:
            logging.debug("Read in %i aligned sequences e.g. %s %s" % (
                len(protein_alignment),
                protein_alignment[0].name,
                protein_alignment[0].seq))
        else:
            logging.debug("No aligned sequences found for this HMM")
        return protein_alignment
