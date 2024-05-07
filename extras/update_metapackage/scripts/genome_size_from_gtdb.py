#!/usr/bin/env python3

###############################################################################
#
#    Copyright (C) 2023 Ben Woodcroft
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
__copyright__ = "Copyright 2023"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3"
__maintainer__ = "Ben Woodcroft"
__email__ = "benjwoodcroft near gmail.com"
__status__ = "Development"

import argparse
import logging
import sys
import os
import polars as pl

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','..','..')] + sys.path

from singlem.genome_size import GenomeSizes

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument(
        '--checkm2-grep',
        help="grep to header remove checkm2 quality file(s). If not specified use checkm2 stats from metadata files")
    parent_parser.add_argument(
        '--gtdb-bac-metadata',
        help="GTDB metadata file for bacteria", required=True)
    parent_parser.add_argument(
        '--gtdb-ar-metadata',
        help="GTDB metadata file for archaea", required=True)

    args = parent_parser.parse_args()

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

    # Read checkm2 stats if needed
    if args.checkm2_grep:
        checkm = pl.read_csv(args.checkm2_grep, has_header=False, separator="\t")
        checkm.columns = ['Name','Completeness','Contamination','Completeness_Model_Used','Additional_Notes']
        checkm = checkm.with_columns(pl.col("Name").str.replace('_protein','').alias("accession"))
        logging.info("Read in {} checkm2 stats".format(len(checkm)))

    # Read GTDB metadata
    columns_of_interest = [
        "accession", 'gtdb_representative', "gtdb_taxonomy", 'genome_size'
    ]
    if not args.checkm2_grep:
        columns_of_interest += ['checkm2_completeness', 'checkm2_contamination']
    gtdb = pl.concat([
        pl.read_csv(args.gtdb_bac_metadata, separator="\t", infer_schema_length=10000000).select(columns_of_interest),
        pl.read_csv(args.gtdb_ar_metadata, separator="\t", infer_schema_length=10000000).select(columns_of_interest)
    ])
    logging.info("Read in {} GTDB entries".format(len(gtdb)))

    # We do not calculate the checkm stats for all genomes in GTDB, just the reps, to save effort.
    gtdb_reps = gtdb.filter(pl.col("gtdb_representative") == "t")
    logging.info("Read in {} GTDB reps".format(len(gtdb_reps)))

    # Calculate the mean completeness and contamination for each GTDB rep
    if args.checkm2_grep:
        gc_reps = gtdb_reps.join(checkm, on="accession", how="inner")
        if len(gc_reps) != len(checkm):
            logging.warning("Different number of checkm2 stats and GTDB reps detected, something is amiss.")
        if len(gc_reps) != len(gtdb_reps):
            raise Exception("It appears there are not enough checkm2 stats to cover all the reps, found only {} checkm2 stats for {} reps".format(len(gc_reps), len(gtdb_reps)))
        logging.info("Joined {} GTDB reps with checkm2 stats".format(len(gc_reps)))
    else:
        gc_reps = gtdb_reps.with_columns(
            pl.col('checkm2_completeness').cast(pl.Float32).alias('Completeness'),
            pl.col('checkm2_contamination').cast(pl.Float32).alias('Contamination'),
        )

    def correct_size(struct: dict):
        return GenomeSizes.corrected_genome_size(struct['genome_size'], struct['Completeness']/100., struct['Contamination']/100.)

    gc_reps = gc_reps.with_columns(
        pl.struct(['genome_size', 'Completeness', 'Contamination']).map_elements(correct_size).alias('checkm2_adjusted_genome_size')
    )

    final = GenomeSizes.calculate_rank_genome_sizes(
        list(gc_reps['gtdb_taxonomy']),
        list(gc_reps['checkm2_adjusted_genome_size'])
    )

    final.write_csv("/dev/stdout", separator="\t")
