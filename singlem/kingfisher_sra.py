import re
import os
import logging

from .sequence_classes import SeqReader

class KingfisherSra:
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

        regex = re.compile('^(.*)\.([012])$')
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


