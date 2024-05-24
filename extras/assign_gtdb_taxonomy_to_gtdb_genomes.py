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

import os
import sys
import argparse
import logging

import pandas as pd

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path
from singlem.otu_table import OtuTable

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser()
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--otu-table-list', help='list of files to directly assign taxonomy to', required=True)
    parent_parser.add_argument('--gtdb-bacterial-metadata', help='metadata for bacterial genomes', required=True)
    parent_parser.add_argument('--gtdb-archaeal-metadata', help='metadata for archaeal genomes', required=True)
    
    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # Read GTDB metadata
    logging.info("Parsing GTDB metadata ..")
    gtdb = pd.concat([
        pd.read_csv(args.gtdb_bacterial_metadata, sep="\t"),
        pd.read_csv(args.gtdb_archaeal_metadata, sep="\t")
    ])
    # gtdb = gtdb[gtdb["gtdb_representative"] == "t"]
    # gtdb["phylum"] = gtdb["gtdb_taxonomy"].apply(lambda x: x.split(";")[1])
    # gtdb["class"] = gtdb["gtdb_taxonomy"].apply(lambda x: x.split(";")[2])
    # gtdb["order"] = gtdb["gtdb_taxonomy"].apply(lambda x: x.split(";")[3])
    # gtdb["family"] = gtdb["gtdb_taxonomy"].apply(lambda x: x.split(";")[4])
    # gtdb["genus"] = gtdb["gtdb_taxonomy"].apply(lambda x: x.split(";")[5])
    # gtdb["sp"] = gtdb["gtdb_taxonomy"].apply(lambda x: x.split(";")[6])
    gtdb_id_to_taxonomy = {}
    for i, row in gtdb.iterrows():
        gtdb_id_to_taxonomy[row["accession"]] = row["gtdb_taxonomy"]
    logging.info("Parsed in {} metadata entries from gtdb".format(len(gtdb_id_to_taxonomy)))

    # Iterate over OTU tables, assigning taxonomy to each line
    is_first = True
    
    taxonomy_field_index = OtuTable.DEFAULT_OUTPUT_FIELDS.index('taxonomy')
    sample_name_field_index = OtuTable.DEFAULT_OUTPUT_FIELDS.index('sample')
    coverage_field_index = OtuTable.DEFAULT_OUTPUT_FIELDS.index('coverage')

    with open(args.otu_table_list) as f: 
        otu_table_files = list(f.readlines())

    for otu_table_file in otu_table_files:
        otu_table_file = otu_table_file.strip()
        logging.debug("Opening {}".format(otu_table_file))
        with open(otu_table_file) as f:
            otu_table = OtuTable.read(f)

        for otu in otu_table.data:
            gtdb_id = otu[sample_name_field_index].replace('_protein','')
            taxonomy = gtdb_id_to_taxonomy[gtdb_id]

            otu[sample_name_field_index] = gtdb_id
            otu[taxonomy_field_index] = taxonomy
            otu[coverage_field_index] = 1.0

        otu_table.write_to(sys.stdout, fields_to_print=OtuTable.DEFAULT_OUTPUT_FIELDS, print_header=is_first)
        is_first = False

    logging.info("Finished processing {} OTU tables".format(len(otu_table_files)))
