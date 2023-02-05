import re
import os
import logging

from .sequence_classes import SeqReader
from .pipe_sequence_extractor import ExtractedReads, ExtractedReadSet

class KingfisherSra:
    def _split_regex(self):
        return re.compile(r'^(.*)\.([012])$')

    def split_fasta(self, fasta_path, output_directory):
        '''fasta_path points to a fasta sequence that contains unordered
        sequences with names like "<seq_id>.X" where X is 0, 1 or 2, which
        signify (unpaired), (forward or unpaired) and (reverse) respectively.

        This function creates new files for (forward or unpaired) and (reverse)
        in the output directory where the .X is removed, and returns a tuple of
        (forward_path, reverse_path)
        '''
        forward_output = None
        reverse_output = None

        regex = self._split_regex()
        forward_count = 0
        reverse_count = 0

        with open(fasta_path) as input:
            for name, seq, _ in SeqReader().readfq(input):
                m = regex.match(name)
                if m is None:
                    raise Exception("Unexpected format for kingfisher SRA readname: {}".format(
                        name
                    ))
                else:
                    new_name = m[1]
                    if m[2] == '1' or m[1] == '0':
                        if forward_output is None:
                            forward_output = open(os.path.join(output_directory, 'forward.fna'),'w')
                        forward_output.write(">{}\n{}\n".format(
                            new_name, seq
                        ))
                        forward_count += 1
                    elif m[2] == '2':
                        if reverse_output is None:
                            reverse_output = open(os.path.join(output_directory, 'reverse.fna'),'w')
                        reverse_output.write(">{}\n{}\n".format(
                            new_name, seq
                        ))
                        reverse_count += 1
                    else:
                        raise Exception("Programming error")

        logging.debug("SRA split found {} forward/unpaired reads and {} reverse reads".format(
            forward_count, reverse_count
        ))

        to_return_forward = None
        to_return_reverse = None
        if forward_output is not None:
            forward_output.close()
            to_return_forward =  forward_output.name
        if reverse_output is not None:
            reverse_output.close()
            to_return_reverse =  reverse_output.name
        return (to_return_forward, to_return_reverse)

    def split_extracted_reads(self, extracted_reads):
        """Given an ExtractedReads object, return a copy of the data that has
        been split into forward and reverse (if necessary) and sequences have
        been renamed accordingly.
        """
        regex = self._split_regex()

        # Go through the sequences objects from all the samples. If any end in
        # .2 we are paired.
        analysing_pairs = False
        for readset in extracted_reads:
            for s in readset.sequences:
                m = regex.match(s.name)
                if m is None:
                    raise Exception("Unexpected format for kingfisher SRA readname: {}".format(
                        s.name
                    ))
                elif m[2] == '2':
                    analysing_pairs = True
                    break
            if analysing_pairs:
                break
        logging.debug("split_extracted_reads: Found analysing pairs {}".format(analysing_pairs))

        to_return = ExtractedReads(analysing_pairs)

        for readset in extracted_reads:
            # Go through the unaligned sequences, renaming them and putting them
            # into the correct thing.
            new_unknown_sequences_forward = []
            new_unknown_sequences_reverse = []
            for u in readset.unknown_sequences:
                m = regex.match(u.name)
                if m is None:
                    raise Exception("Unexpected format for kingfisher SRA readname: {}".format(
                        u.name
                    ))
                u.name = m[1]
                if analysing_pairs and m[2] == '2':
                    new_unknown_sequences_reverse.append(u)
                else:
                    new_unknown_sequences_forward.append(u)

            # Rename the sequences as well
            new_sequences_forward = []
            new_sequences_reverse = []
            for u in readset.sequences:
                m = regex.match(u.name)
                if m is None:
                    raise Exception("Unexpected format for kingfisher SRA readname: {}".format(
                        u.name
                    ))
                u.name = m[1]
                if analysing_pairs and m[2] == '2':
                    new_sequences_reverse.append(u)
                else:
                    new_sequences_forward.append(u)
            
            if analysing_pairs:
                to_return.add([
                    ExtractedReadSet(
                        readset.sample_name,
                        readset.singlem_package,
                        new_sequences_forward,
                        readset.known_sequences,
                        new_unknown_sequences_forward),
                    ExtractedReadSet(
                        readset.sample_name,
                        readset.singlem_package,
                        new_sequences_reverse,
                        readset.known_sequences,
                        new_unknown_sequences_reverse)
                ])
            else:
                to_return.add(ExtractedReadSet(
                    readset.sample_name,
                    readset.singlem_package,
                    new_sequences_forward,
                    readset.known_sequences,
                    new_unknown_sequences_forward
                ))
        
        return analysing_pairs, to_return


