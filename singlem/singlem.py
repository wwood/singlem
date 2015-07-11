from Bio.Seq import Seq
import itertools
import logging
import re
import os
import csv
            
            
class Sequence:
    '''Simple name+sequence object'''
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
        
class AlignedProteinSequence:
    def __init__(self, name, seq):
        self.name = name
        self.seq = seq
    
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
        
class AlignedNucleotideSequence:
    '''Represent a nucleotide sequence aligned in protein space, together with the 
    nucleotide sequence that it came from'''
    
    def __init__(self, name, aligned_sequence, unaligned_sequence, aligned_length):
        '''
        Parameters
        ---------
        name: str
            name of the sequence
        aligned_sequence: str
            aligned nucleotide sequence
        unaligned_sequence: str
            unaligned nucleotide sequence
        aligned_length:
            the number of nucleotides used in the alignment, including columns
            that were removed as not aligned
        '''
        self.name = name
        self.aligned_sequence = aligned_sequence
        self.unaligned_sequence = unaligned_sequence
        self.aligned_length = aligned_length
        
    def coverage_increment(self):
        '''Given the alignment came from a read of length 
        original_nucleotide_sequence_length, how much coverage does the 
        observation of this aligned sequence indicate?'''
        return float(len(self.unaligned_sequence))/\
            (len(self.unaligned_sequence)-self.aligned_length+1)
       
class OrfMUtils:
    def un_orfm_name(self, name):
        return re.sub('_\d+_\d+_\d+$', '', name)
       
class TaxonomyFile:
    def __init__(self, taxonomy_file_path):
        self.sequence_to_taxonomy = {}
        utils = OrfMUtils()
        with open(taxonomy_file_path) as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                self.sequence_to_taxonomy[\
                      utils.un_orfm_name(row[0])] = row[1]
                      
    def __getitem__(self, item):
        return self.sequence_to_taxonomy[item]
        
class HmmDatabase:
    def __init__(self):
        # Array of gpkg names to HmmAndPostion
        self.hmms_and_positions = {}
                                     
        for array in [
         ['ribosomal_protein_S2_rpsB_gpkg','DNGNGWU00001', 20],
         ['ribosomal_protein_L11_rplK_gpkg','DNGNGWU00024', 15],
         ['ribosomal_protein_S17_gpkg','DNGNGWU00036', 51],
          ]:
            hmm_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                            '..', 'db', array[0])
            self.hmms_and_positions[os.path.basename(hmm_path)] = \
                HmmAndPostion(hmm_path,                                         
                               os.path.join(hmm_path, array[1])+".hmm",
                               array[2]
                               )
            
    def hmm_paths(self):
        'return an array of absolute paths to the hmms in this database'
        return [hp.hmm_filename for hp in self.hmms_and_positions.values()]
    
    def gpkg_basenames(self):
        return self.hmms_and_positions.keys()
    
    def gpkg_paths(self):
        return [h.gpkg_path for _, h in self.hmms_and_positions.iteritems()]
        
class HmmAndPostion:
    def __init__(self, gpkg_path, hmm_filename, best_position):
        self.gpkg_path = gpkg_path
        self.hmm_filename = hmm_filename
        self.best_position = best_position
        
    def hmm_path(self):
        return os.path.join(self.gpkg_path, self.hmm_filename)
        
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
                name = s.un_orfm_name()
                nuc = nucleotide_sequences[name]
                align, aligned_length = self._nucleotide_alignment(\
                    s, s.orfm_nucleotides(nuc), chosen_positions)
                
                windowed_sequences.append(AlignedNucleotideSequence(name, align, nuc, aligned_length))
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
        at the chosen_positions, and the length in nucleotides that the chosen
        positions stretch across.
        
        Assumes the chosen positions are ascending monotonically and has length
        >0'''
        
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
        
        aligned_length = 0
        for i in range(chosen_positions[0],chosen_positions[-1]+1):
            if protein_sequence.seq[i] != '-':
                aligned_length += 3
                
        return ''.join(itertools.chain(codons[i] for i in chosen_positions)), \
            aligned_length

