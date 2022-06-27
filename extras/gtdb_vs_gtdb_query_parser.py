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

import csv

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path
from singlem.otu_table import OtuTable

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser()
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--gtdb-vs-gtdb-query', help='GTDB vs GTDB query result', required=True)
    parent_parser.add_argument('--max-lineage-matching', type=int, help='minimum difference in lineage between GTDB and GTDB query e.g. 5 for genus')
    
    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    # So now want to parse those lines, looking for, for each hit:

    # closest hit sequence distance
    # taxonomy difference

    # Then get a (marker, distance, taxon_level_difference, count) out at the end,
    # both for nearest hit and generally. Then graph that to work out the average
    # distance of the nearest same-genus, same-family, etc. for each of the markers.
    # Then should get a distribution of same-genus, and anything at say less than
    # the 95%+ difference on genus level can be designated as same species when
    # assigning taxonomy.

    headers = None
    state = "waiting_for_self_hit"
    self_hit_taxonomy = None
    results = {}
    
    with open(args.gtdb_vs_gtdb_query) as f:
        for row in csv.reader(f, delimiter='\t'):
            if headers is None:
                headers = row
                # query_name      query_sequence  divergence      num_hits        sample  marker  hit_sequence    taxonomy
                divergence_index = headers.index('divergence')
                marker_index = headers.index('marker')
                hit_taxonomy_index = headers.index('taxonomy')
                query_name_index = headers.index('query_name')
                hit_name_index = headers.index('sample')
                continue

            divergence = int(row[divergence_index])
            
            if row[query_name_index] == row[hit_name_index]:
                state = "finished_self_hit"
                self_hit_taxonomy = row[hit_taxonomy_index]
            elif state == "finished_self_hit":
                # best non-self hit
                query_taxonomy2 = self_hit_taxonomy.split(';')
                hit_taxonomy2 = row[hit_taxonomy_index].split(';')
                
                num_matching = 0
                for i in range(len(query_taxonomy2)):
                    if query_taxonomy2[i] == hit_taxonomy2[i]:
                        num_matching += 1
                    else:
                        break

                if args.max_lineage_matching is not None and num_matching > args.max_lineage_matching:
                    continue

                state = "finished_best_hit"
                # if row[marker_index] == 'S3.10.ribosomal_protein_S19_rpsS' and query_taxonomy2[0] == 'd__Bacteria':
                logging.debug('query {} hit {} divergence {}\t{}'.format(row[query_name_index], row[hit_name_index], divergence, num_matching))

                key = '{}:{}:{}:{}'.format(row[marker_index], query_taxonomy2[0], divergence, num_matching)
                if key in results:
                    results[key] += 1
                else:
                    results[key] = 1
            else:
                # Could technically be a equal-best non-self hit, but eh,
                # that's probably rare enough it doesn't affect results.
                pass

    print('\t'.join(['marker', 'query_domain', 'divergence', 'num_matching', 'count']))
    for (key, count) in results.items():
        print('{}\t{}'.format(key.replace(':','\t'), count))
                
                

