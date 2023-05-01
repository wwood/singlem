#!/usr/bin/env python3

__prog_name__ = 'sim_reads.py'
__prog_desc__ = 'simulate reads across a fasta file.'

__author__ = 'David Wood'
__copyright__ = 'Copyright 2018'
__credits__ = ['Donovan Parks, David Wood, Ben Woodcroft']
__license__ = 'GPL3'
__version__ = '0.0.2'
__maintainer__ = 'David Wood'
__email__ = 'davidlawood@gmail.com'
__status__ = 'Development'

import os
import sys
import csv
import shutil
import tempfile
import argparse
import uuid
import logging

csv.field_size_limit(sys.maxsize)

class SequenceIO:
    # Stolen from https://github.com/lh3/readfq/blob/master/readfq.py
    def each(self, fp): # this is a generator function
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

class ReadSimulator(object):

    @staticmethod
    def simulate_reads(input_io, read_len, read_step_size, insert_size, fout1, fout2):
        """Simulate reads."""
        fragment_len = 2*read_len + insert_size
        pair_count = 0
        for seq_id, seq, _ in SequenceIO().each(input_io):
            for i in range(0, len(seq) - fragment_len, read_step_size):
                pair_count += 1
                fout1.write('>%d\n' % pair_count)
                fout1.write('%s\n' % seq[i:i+read_len])
                fout2.write('>%d\n' % pair_count)
                fout2.write('%s\n' % seq[i+read_len+insert_size:i+2*read_len+insert_size])
        print("Write %i read pairs" % pair_count)


if __name__ == '__main__':
    print(__prog_name__ + ' v' + __version__ + ': ' + __prog_desc__)

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', help='file to simulate reads across', required=True)
    parser.add_argument('-1', '--output_read1', help='output read1 file', required=True)
    parser.add_argument('-2', '--output_read2', help='output read2 file', required=True)

    parser.add_argument('--read_len', type=int, help='length of simulated reads', default=150)
    parser.add_argument('--step_size', type=int, help='simulate a read each X bp', default=10)
    parser.add_argument('--insert_size', type=int, help='insert length of simulated reads', default=200)
    args = parser.parse_args()

    with open(args.i) as in_io:
        with open(args.output_read1,'w') as o1:
            with open(args.output_read2,'w') as o2:
                ReadSimulator.simulate_reads(in_io, args.read_len, args.step_size, args.insert_size, o1, o2)
