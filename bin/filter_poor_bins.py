#!/usr/bin/env python

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2017"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3+"
__maintainer__ = "Ben Woodcroft"
__email__ = "b.woodcroft near uq.edu.au"
__status__ = "Development"

import argparse
import logging
import os
import sys
import csv
import tempfile
import extern
import re

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

import singlem.pipe as pipe
from singlem.querier import Querier, QueryInputSequence
from singlem.sequence_database import SequenceDatabase
from singlem.sequence_classes import SeqReader
from singlem.otu_table import OtuTable
from singlem.checkm import CheckM

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--debug', help='output debug information', action="store_true")
    #parser.add_argument('--version', help='output version information and quit',  action='version', version=singlem.__version__)
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument('--checkm_csv', required=True)
    parser.add_argument('--bin_otu_table', required=True)

    parser.add_argument('--min_completeness', type=int, default=80)
    parser.add_argument('--max_contamination', type=int, default=5)

    args = parser.parse_args()

    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # Read checkm CSV file
    with open(args.checkm_csv) as f:
        quality_bins = CheckM.read_checkm_stats(
            f, args.min_completeness, args.max_contamination)
    logging.info("Found %i bin(s) of sufficient quality according to CheckM statistics" % len(quality_bins))

    to_print = OtuTable()

    with open(args.bin_otu_table) as b:
        for otu in OtuTable.each(b):
            if otu.sample_name in quality_bins:
                to_print.add([otu])
    logging.info("Printing {} OTUs..".format(len(to_print.data)))
    to_print.write_to(sys.stdout)
