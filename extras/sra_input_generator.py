#!/usr/bin/env python3

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2021"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3+"
__maintainer__ = "Ben Woodcroft"
__email__ = "b.woodcroft near uq.edu.au"
__status__ = "Development"

import argparse
import logging
import sys
import os
import math
import re

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

import singlem
from singlem.orf_length_checker import OrfLengthChecker

# Works out which files should be used in a singlem run, when data is input from fasterq-dump

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--version', help='output version information and quit',  action='version', version=singlem.__version__)
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument('--min-orf-length',default=72,type=int)
    parser.add_argument('--fastq-dump-outputs',required=True,nargs='+')
    args = parser.parse_args()

    if len(args.fastq_dump_outputs) > 3:
        raise Exception("Unexpectedly found >3 input files")

    # Go through each file and determine which have sufficient ORF length
    passable_paths = []
    for path in args.fastq_dump_outputs:
        if OrfLengthChecker.check_sequence_file_contains_an_orf(
            path, args.min_orf_length):

            logging.info("{} contains 1 or more ORFs of length {}".format(
                path, args.min_orf_length
            ))
            passable_paths.append(path)
        else:
            logging.warning("{} contained no ORFs of length {}".format(
                path, args.min_orf_length
            ))

    forwards = []
    reverses = []
    unpaireds = []

    unpaired_re = re.compile('.*\d+\.fastq')
    for path in passable_paths:
        if path.endswith('_1.fastq') or path.endswith('_1.fastq.gz'):
            forwards.append(path)
        elif path.endswith('_2.fastq') or path.endswith('_2.fastq.gz'):
            reverses.append(path)
        elif unpaired_re.match(path):
            unpaireds.append(path)
        else:
            raise Exception("Unexpected file name {}".format(path))

    if len(forwards) > 1: raise Exception(">1 forward")
    if len(reverses) > 1: raise Exception(">1 reverse")
    if len(unpaireds) > 1: raise Exception(">1 unpaired")

    forward = forwards[0] if len(forwards) == 1 else None
    reverse = reverses[0]  if len(reverses) == 1 else None
    unpaired = unpaireds[0] if len(unpaireds) == 1 else None

    if forward and reverse:
        print('--forward {} --reverse {}'.format(forward, reverse))
    elif forward and not reverse:
        print('--forward {}'.format(forward))
    elif reverse and not forward:
        print('--forward {}'.format(reverse))
    elif unpaired:
        print('--forward {}'.format(unpaired))
    else:
        logging.warn("No passable sequence files found")
        # Print nothing