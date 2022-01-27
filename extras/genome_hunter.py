#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2020 Ben Woodcroft
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2022"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import argparse
from enum import unique
import logging
import sys
import os
import pandas as pd
from io import StringIO
import re

import extern

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    #parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument('--protein-query-against-reads',required=True)
    parser.add_argument('--sample-list',required=True)
    parser.add_argument('--assembly-files',nargs='+',required=True)
    parser.add_argument('--metapackage',required=True)

    args = parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # Read query OTU table
    # Remove any OTUs that are not in the sample list
    protein_queries = pd.read_csv(args.protein_query_against_reads,sep='\t')
    with open(args.sample_list) as f:
        samples_of_interest = set([s.strip() for s in f.readlines()])
    logging.info("Found {} samples of interest".format(len(samples_of_interest)))
    target_sequences = set(pd.merge(pd.DataFrame({'sample':list(samples_of_interest)}), protein_queries, on='sample')['hit_sequence'])
    logging.info("Found {} target sequences from those samples of interest".format(len(target_sequences)))
    if len(target_sequences) == 0:
        raise Exception("No target sequences found")
    

    # Run SingleM pipe on assembly files, putting results in STDOUT and reading, I guess
    assembly_file_otu_tables = {}
    for assembly_file in args.assembly_files:
        logging.info("Running SingleM on {}".format(assembly_file))
        cmd = f'singlem pipe --genome-fasta-files {assembly_file} --otu_table /dev/stdout --threads 1 --no-assign-taxonomy --output-extras --singlem-metapackage {args.metapackage}'
        logging.debug("Running {}".format(cmd))
        output = extern.run(cmd)
        assembly_file_otu_tables[assembly_file] = pd.read_csv(StringIO(output),sep='\t')
    logging.info("Finished running SingleM on {} assemblies".format(len(assembly_file_otu_tables)))

    # Find OTUs in assembly singlem that are in the trimmed query OTU table
    good_otus = None
    for assembly_file,otu_table in assembly_file_otu_tables.items():
        logging.debug("Inspecting OTU table for {}".format(assembly_file))
        found_otus = pd.merge(otu_table, pd.DataFrame({'sequence':list(target_sequences)}), on='sequence')
        logging.info("Found {} assembled matching OTUs in {}".format(len(found_otus), assembly_file))
        new_frame = pd.DataFrame(found_otus)
        new_frame['assembly_file'] = assembly_file

        if good_otus is None:
            good_otus = new_frame
        else:
            good_otus = pd.concat(good_otus, new_frame)
    logging.info("Found {} matching OTUs".format(len(good_otus)))

    # Print the length of the contigs that match OTU sequences in the query OTU table
    last_assembly_file = None
    last_assembly_file_contig_lengths = None
    r = re.compile('(.*)_\d+_\d+_\d+')
    for otu_tuple in good_otus.iterrows():
        otu = otu_tuple[1]

        if last_assembly_file is None or last_assembly_file != otu['assembly_file']:
            last_assembly_file = otu['assembly_file']
            logging.info("Finding contig lengths from {} ..".format(last_assembly_file))
            last_assembly_file_contig_lengths_str = extern.run("bbbin describe <{}".format(last_assembly_file))
            last_assembly_file_contig_lengths = pd.read_csv(StringIO(last_assembly_file_contig_lengths_str),sep='\t',header=None)

        # Read sequence length of each contig in the original assembly file - unaligned sequences in the OTU table are from the transcripts, not the contigs
        assembly_file = otu['assembly_file']

        sequence_names = otu['read_names'].split(' ')
        for name in sequence_names:
            contig_name = r.match(name).group(1)
            sequence_length = last_assembly_file_contig_lengths[last_assembly_file_contig_lengths[0] == contig_name][1].values[0]

            print("\t".join((otu['assembly_file'], str(sequence_length), '\t'.join([str(s) for s in otu.values]))))
