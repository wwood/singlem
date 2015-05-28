#!/usr/bin/env python

from Bio.Seq import Seq
import itertools
import logging
import re
import os
            
            
class Sequence:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        
    def unaligned_length(self):
        return len(re.sub('-','',self.seq))
        
    def un_orfm_name(self):
        return re.sub('_\d+_\d+_\d+$', '', self.name)
    
    def orfm_nucleotides(self, original_nucleotide_sequence):
        m = re.search('_(\d+)_(\d+)_\d+$', self.name)
        start = int(m.groups(0)[0])-1
        translated_seq = original_nucleotide_sequence[start:(start+3*self.unaligned_length())]
        logging.debug("Returning orfm nucleotides %s" % translated_seq)
        if int(m.groups(0)[1]) > 3:
            # revcomp type frame
            return(str(Seq(translated_seq).reverse_complement()))
        else:
            return(translated_seq)
        
        
class HmmDatabase:
    def __init__(self):  
        self.hmms_and_positions = {}
        self.hmm_directory = os.path.join(
                                     os.path.dirname(os.path.realpath(__file__)),
                                     '..',
                                     'hmms')
                                     
        for array in [
         ['DNGNGWU00001', 20],
         ['DNGNGWU00024', 15],
         ['DNGNGWU00036', 51],
          ]:
            self.hmms_and_positions[array[0]] = HmmAndPostion(
                               os.path.join(self.hmm_directory, array[0])+".hmm",
                               array[1]
                               )
            
            
    def hmm_paths(self):
        'return an array of absolute paths to the hmms in this database'
        return [hp.hmm_filename for hp in self.hmms_and_positions.values()]
    
    def hmm_basenames(self):
        return self.hmms_and_positions.keys()
    
    def base_directory(self):
        return self.hmm_directory
    
    def best_position(self, hmm_basename):
        return self.hmms_and_positions[hmm_basename].best_position
            
        
class HmmAndPostion:
    def __init__(self, hmm_filename, best_position):
        self.hmm_filename = hmm_filename
        self.best_position = best_position
        
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
    
class MetagenomeOtuFinder:
    def __init__(self):
        pass
    
    def find_windowed_sequences(self,
                                aligned_sequences,
                                nucleotide_sequences,
                                stretch_length,
                                best_position=None,
                                ):
        ignored_columns = self._find_lower_case_columns(aligned_sequences)
        logging.debug("Ignoring columns %s", str(ignored_columns))
        
        # Find the stretch of the protein that has the most number of aligned bases in a 20 position stretch,
        # excluding sequences that do not aligned to the first and last bases
        if best_position:
            start_position = self._upper_case_position_to_alignment_position(best_position, ignored_columns)
            logging.info("Using pre-defined best section of the alignment starting from %i" % (start_position+1))
        else:
            start_position = self._find_best_window(aligned_sequences, stretch_length, ignored_columns)
            logging.info("Found best section of the alignment starting from %i" % (start_position+1))
        
        chosen_positions = self._best_position_to_chosen_positions(start_position, stretch_length, ignored_columns)
        logging.debug("Found chosen positions %s", chosen_positions)
        
        # For each read aligned to that region i.e. has the first and last bases,
        # print out the corresponding nucleotide sequence
        windowed_sequences = []
        for s in aligned_sequences:
            if s.seq[chosen_positions[0]] != '-' and s.seq[chosen_positions[-1]] != '-':
                nuc = nucleotide_sequences[s.un_orfm_name()]
                windowed_sequences.append(
                                          self._nucleotide_alignment(s, s.orfm_nucleotides(nuc), chosen_positions)
                                          )
        return windowed_sequences
    
    def _find_lower_case_columns(self, protein_alignment):
        lower_cases = [False]*len(protein_alignment[0].seq)
        lower_case_chars = re.compile(r'[a-z]')
        for pro in protein_alignment:
            for i, aa in enumerate(pro.seq):
                if lower_case_chars.match(aa):
                    lower_cases[i] = True
        return [i for i, is_lower in enumerate(lower_cases) if is_lower]
    
    def _find_best_window(self, protein_alignment, stretch_length, ignored_columns):
        '''Return the position in the alignment that has the most bases aligned
        only counting sequences that overlap the entirety of the stretch.'''
        
        # Convert the alignment into a True/False matrix for ease,
        # True meaning that there is something aligned, else False
        binary_alignment = []
        for s in protein_alignment:
            aln = []
            for i, base in enumerate(s.seq):
                if base=='-':
                    aln.append(False)
                else:
                    aln.append(True)
            binary_alignment.append(aln)
        
        # Find the number of aligned bases at each position
        current_best_position = 0
        current_max_num_aligned_bases = 0
        for i in range(0, len(binary_alignment[0])-stretch_length+1):
            if i in ignored_columns: continue #don't start from ignored columns
            positions = self._best_position_to_chosen_positions(i, stretch_length, ignored_columns)
            logging.debug("Testing positions %s" % str(positions))
            num_bases_covered_here = 0
            end_index = positions[-1]
            if end_index >= len(binary_alignment[0]): continue #if ignored char is within stretch length of end of aln
            for s in binary_alignment:
                if not s[i] or not s[end_index]: continue #ignore reads that don't cover the entirety
                num_bases_covered_here += sum(s[i:end_index])
            logging.debug("Found %i aligned bases at position %i" % (num_bases_covered_here, i))
            if num_bases_covered_here > current_max_num_aligned_bases:
                current_best_position = i
                current_max_num_aligned_bases = num_bases_covered_here
        logging.info("Found a window starting at position %i with %i bases aligned" % (current_best_position,
                                                                                       current_max_num_aligned_bases
                                                                                       ))
        return current_best_position
    
    def _best_position_to_chosen_positions(self, best_position, stretch_length, ignored_columns):
        '''Given a position to start from, and the number of positions to index,
        return the consecutive indices that are not in the ignored_columns list'''
        chosens = []
        i = best_position
        while len(chosens) < stretch_length:
            if i not in ignored_columns:
                chosens.append(i)
            i += 1
        return chosens
    
    def _upper_case_position_to_alignment_position(self, position, ignored_columns):
        target = position
        for i in ignored_columns:
            if i <= target:
                target += 1
            else:
                return target
        return target
    
    def _nucleotide_alignment(self, protein_sequence, nucleotides, chosen_positions):
        '''Line up the nucleotides and the proteins, and return the alignment 
        across start_position for stretch length amino acids'''
        
        codons = []
        # For each position in the amino acid sequence
        # If non-dash character, take 3 nucleotides off the nucleotide sequence and
        # add that as the codon
        # else add None
        for aa in protein_sequence.seq:
            if aa=='-':
                codons.append('---')
            else:
                if len(nucleotides) < 3: raise Exception("Insufficient nucleotide length found")
                codons.append(nucleotides[:3])
                if len(nucleotides)>2: nucleotides = nucleotides[3:]
        if len(nucleotides) > 0: raise Exception("Insufficient protein length found")
        
        return ''.join(itertools.chain(codons[i] for i in chosen_positions))

