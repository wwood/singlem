import logging
import re
from sequence_classes import AlignedNucleotideSequence
import itertools

class MetagenomeOtuFinder:
    def __init__(self):
        pass

    def find_windowed_sequences(self,
                                aligned_sequences,
                                nucleotide_sequences,
                                stretch_length,
                                include_inserts,
                                is_protein_alignment,
                                best_position=None,
                                ):
        '''Return an array of AlignedNucleotideSequence objects containing sequences
        aligned at the given best_position, finding that best position if None
        is given.

        Parameters
        ----------

        '''
        ignored_columns = self._find_lower_case_columns(aligned_sequences)
        logging.debug("Ignoring columns %s", str(ignored_columns))

        # Find the stretch of the protein that has the most number of aligned bases in a 20 position stretch,
        # excluding sequences that do not aligned to the first and last bases
        if best_position:
            start_position = self._upper_case_position_to_alignment_position(best_position, ignored_columns)
            logging.debug("Using pre-defined best section of the alignment starting from %i" % (start_position+1))
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
                align, aligned_length = self._nucleotide_alignment(
                    s, s.orfm_nucleotides(nuc), chosen_positions, True, include_inserts=include_inserts)

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

    def _nucleotide_alignment(self,
                              protein_sequence,
                              nucleotides,
                              chosen_positions,
                              is_protein_alignment,
                              include_inserts=False):
        '''Line up the nucleotides and the proteins, and return the alignment
        at the chosen_positions, and the length in nucleotides that the chosen
        positions stretch across.

        Assumes the chosen positions are ascending monotonically and have length
        >0
        
        Parameters
        ----------
        protein_sequence: AlignedProteinSequence
            aligned amino acid sequence of the aligned ORF
        nucleotides: str
            nucleotide sequence of the unaligned ORF
        chosen_positions: list of int
            positions to return, in ascending order
        is_protein_alignment: boolean
            True for proteins, False for nucleotide alignments
        include_inserts: boolean
            if False, remove inserts in nucleotide sequence relative to the 
            alignment. If True, include them
            
        Returns
        -------
        list of 2: the nucleotides string, and the length of nucleotide sequence
            used to cover the alignment
        '''
        if is_protein_alignment:
            length_ratio = 3
            empty_codon = '---'
        else:
            length_ratio = 1
            empty_codon = '-'
            
        codons = []
        # For each position in the amino acid sequence
        # If non-dash character, take 3 nucleotides off the nucleotide sequence and
        # add that as the codon
        # else add None
        for aa in protein_sequence.seq:
            if aa=='-':
                codons.append(empty_codon)
            else:
                if len(nucleotides) < length_ratio: raise Exception("Insufficient nucleotide length found")
                codons.append(nucleotides[:length_ratio])
                if len(nucleotides)>=length_ratio: nucleotides = nucleotides[length_ratio:]
                if nucleotides[:length_ratio] == empty_codon: raise Exception("Input nucleotide sequence had gap characters, didn't expect this")
        if len(nucleotides) > 0: raise Exception("Insufficient protein length found")

        aligned_length = 0
        for i in range(chosen_positions[0],chosen_positions[-1]+1):
            if protein_sequence.seq[i] != '-':
                aligned_length += length_ratio

        if include_inserts:
            to_return = []
            for i in range(chosen_positions[0], chosen_positions[-1]+1):
                if i in chosen_positions:
                    to_return += codons[i]
                elif codons[i] == empty_codon:
                    pass
                else:
                    to_return += codons[i].lower()
            return ''.join(to_return), aligned_length
        else:
            return ''.join(itertools.chain(codons[i] for i in chosen_positions)), \
                aligned_length
