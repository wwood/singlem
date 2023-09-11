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

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path

if __name__ == '__main__':
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    #parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=repeatm.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parent_parser.add_argument(
        '--checkm2-grep',
        help="grep to header remove checkm2 quality file(s)", required=True)
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

    # Read checkm2 stats
    checkm = pl.read_csv(args.checkm2_grep, has_header=False, separator="\t")
    checkm.columns = ['Name','Completeness','Contamination','Completeness_Model_Used','Additional_Notes']
    checkm = checkm.with_columns(pl.col("Name").str.replace('_protein','').alias("accession"))
    logging.info("Read in {} checkm2 stats".format(len(checkm)))

    # Read GTDB metadata
    columns_of_interest = [
        "accession", "gtdb_taxonomy", "gtdb_representative", 'genome_size'
    ]
    gtdb = pl.concat([
        pl.read_csv(args.gtdb_bac_metadata, separator="\t", infer_schema_length=10000000).select(columns_of_interest),
        pl.read_csv(args.gtdb_ar_metadata, separator="\t", infer_schema_length=10000000).select(columns_of_interest)
    ])
    logging.info("Read in {} GTDB reps".format(len(gtdb)))

    # We do not calculate the checkm stats for all genomes in GTDB, just the reps, to save effort.
    gtdb_reps = gtdb.filter(pl.col("gtdb_representative") == "t")
    gtdb_non_reps = gtdb.filter(pl.col("gtdb_representative") == "f")
    if len(gtdb_reps) + len(gtdb_non_reps) != len(gtdb):
        raise Exception("Unexpected GTDB metadata")

    # Calculate the mean completeness and contamination for each GTDB rep
    gc_reps = gtdb_reps.join(checkm, on="accession", how="inner")
    if len(gc_reps) != len(checkm):
        raise Exception()
    gc_reps = gc_reps.with_columns(
        (pl.col('genome_size')/(1+pl.col('Contamination')/100.)/(pl.col('Completeness')/100.)).alias('checkm2_adjusted_genome_size')
    )

    # Trust the genome sizes for the rest as-is
    gc_non_reps = gtdb_non_reps.with_columns(
        pl.col('genome_size').alias('checkm2_adjusted_genome_size').cast(pl.Float64)
    )

    gc = pl.concat([
        gc_reps.select(['accession', 'checkm2_adjusted_genome_size', 'gtdb_taxonomy']),
        gc_non_reps.select(['accession', 'checkm2_adjusted_genome_size', 'gtdb_taxonomy']),
    ])

    gc = gc.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(0).alias("kingdom"))
    gc = gc.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(1).alias("phylum"))
    gc = gc.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(2).alias("class"))
    gc = gc.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(3).alias("order"))
    gc = gc.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(4).alias("family"))
    gc = gc.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(5).alias("genus"))
    gc = gc.with_columns(pl.col("gtdb_taxonomy").str.split(';').list.get(6).alias("species"))
    logging.info("GTDB taxonomy parsed.")

    levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    final = None

    for level in levels:
        taxon_means = gc.select([level, 'checkm2_adjusted_genome_size']).group_by(level).mean()
        taxon_means.columns = ['rank', 'genome_size']
        taxon_mins = gc.select([level, 'checkm2_adjusted_genome_size']).group_by(level).min()
        taxon_mins.columns = ['rank', 'min_genome_size']
        taxon_maxs = gc.select([level, 'checkm2_adjusted_genome_size']).group_by(level).max()
        taxon_maxs.columns = ['rank', 'max_genome_size']
        taxon_99pct = gc.select([level, 'checkm2_adjusted_genome_size']).group_by(level).quantile(0.99, interpolation='linear')
        taxon_99pct.columns = ['rank', '99pct_genome_size']

        new_df = taxon_means.join(taxon_mins, on='rank', how="inner").join(taxon_maxs, on='rank', how="inner").join(taxon_99pct, on='rank', how="inner")

        if final is None:
            final = new_df
        else:
            final = pl.concat([
                final,
                new_df
            ])
    logging.info("Calculated genome size statistics for {} taxa".format(len(final)))

    final.write_csv("/dev/stdout", separator="\t")
