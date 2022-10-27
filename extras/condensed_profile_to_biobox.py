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
import logging
import sys
import os

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path
from singlem.condense import CondensedCommunityProfile

def write_biobox(input_condensed_table_path, output_biobox_path, template_biobox_path, no_fill):
    # # Taxonomic Profiling Output
    # @SampleID:SAMPLEID
    # @Version:0.9.1
    # @Ranks:superkingdom|phylum|class|order|family|genus|species
    # @TaxonomyID:ncbi-taxonomy_DATE
    # @@TAXID	RANK	TAXPATH	TAXPATHSN	PERCENTAGE
    # 2	superkingdom	2	Bacteria	98.81211
    # 2157	superkingdom	2157	Archaea	1.18789
    # 1239	phylum	2|1239	Bacteria|Firmicutes	59.75801
    # 1224	phylum	2|1224	Bacteria|Proteobacteria	18.94674
    # 28890	phylum	2157|28890	Archaea|Euryarchaeotes	1.18789

    levels = ['kingdom','phylum','class','order','family','genus','species']
    total_coverages = [0.0]*len(levels)

    with open(output_biobox_path,'w') as out:
        for s in [
            "# Taxonomic Profiling Output",
            "@SampleID:SAMPLEID",
            "@Version:0.9.1",
            "@Ranks:kingdom|phylum|class|order|family|genus|species",
            "@TaxonomyID:ncbi-taxonomy_DATE",
            "@@TAXID	RANK	TAXPATH	TAXPATHSN	PERCENTAGE"
        ]:
            print(s, file=out)
            
        with open(input_condensed_table_path) as f:
            for condensed_table in CondensedCommunityProfile.each_sample_wise(f):
                # Use template biobox if required
                tax_numbers = {}
                if template_biobox_path:
                    with open(template_biobox_path) as template:
                        for line in template:
                            if not line.startswith('@') and not line.startswith('#'):
                                splits = line.strip().split('\t')
                                if len(splits) != 5:
                                    raise Exception("Unexpected biobox line in template: {}".format(line))
                                taxpathsn = splits[3]
                                tax_numbers[taxpathsn] = len(tax_numbers)+1
                    logging.info("Read in {} taxonomic ranks from template biobox".format(len(tax_numbers)))
                    
                # Pass 1 - collect total coverage for top level
                for wordnode in condensed_table.breadth_first_iter():
                    level = wordnode.calculate_level()-1
                    if level == -1: continue
                    if level > 0: break # Only need total coverage in entire sample
                    if no_fill:
                        total_coverages[level] += wordnode.coverage
                    else:
                        total_coverages[level] += wordnode.get_full_coverage()

                total_coverage = total_coverages[0]
                # Pass 2 - calculate percents
                for wordnode in condensed_table.breadth_first_iter():
                    # Percent includes unassigned
                    # logging.debug("calculating for taxonomy: {}".format(wordnode.get_taxonomy()))
                    level = wordnode.calculate_level()-1
                    if level == -1: continue
                    if no_fill:
                        cov = wordnode.coverage
                    else:
                        cov = wordnode.get_full_coverage()
                    percent = cov/total_coverages[0]
                    taxonomy = wordnode.get_taxonomy()[1:]
                    taxpathsn = '|'.join(taxonomy)
                    if taxpathsn not in tax_numbers: # If not cached
                        tax_numbers[taxpathsn] = len(tax_numbers)+1
                    taxpaths = []
                    for i, t in enumerate(taxonomy):
                        taxpaths.append(tax_numbers['|'.join(taxonomy[:(i+1)])])
                    taxpath = '|'.join([str(tp) for tp in taxpaths])
                    # @@TAXID	RANK	TAXPATH	TAXPATHSN	PERCENTAGE
                    print('\t'.join([
                        str(taxpaths[-1]),
                        levels[level],
                        taxpath,
                        taxpathsn,
                        str(percent*100.)
                    ]), file=out)


if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument('--input-condensed-table', help='input condensed table', required=True)
    parent_parser.add_argument('--output-biobox', help='output biobox', required=True)
    parent_parser.add_argument('--template-biobox', help='template biobox')
    parent_parser.add_argument('--no-fill', help='Use .coverage not get_full_coverage()', action="store_true")

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')

    write_biobox(args.input_condensed_table, args.output_biobox, args.template_biobox, args.no_fill)