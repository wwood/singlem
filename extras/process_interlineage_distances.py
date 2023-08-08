#!/usr/bin/env python3

"""
Author: Samuel Aroney
Process output from inter_lineage_distance_calculator.py
"""

import sys
import argparse
import logging
import polars as pl

def main(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--debug', help='output debug information', action="store_true")
    parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser.add_argument("--interlineage-distances", help="File containing interlineage distances", required=True)
    parser.add_argument("--output", help="Output json file with package-domain cutoffs", required=True)

    args = parser.parse_args(arguments)

    # Setup logging
    if args.debug:
        loglevel = logging.DEBUG
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.INFO
    logging.basicConfig(level=loglevel, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

    trim = 0.1
    trimmed = (pl.col("distance") >= pl.col("distance").quantile(trim)) | (pl.col("distance") <= pl.col("distance").quantile(1 - trim))
    distances = (
        pl.read_csv(args.interlineage_distances, separator="\t")
        .groupby("rank", "singlem_package", "domain", "taxon")
        .agg(distance = pl.col("distance").filter(trimmed).mean())
        .groupby("rank", "singlem_package", "domain")
        .agg(distance = pl.col("distance").filter(trimmed).mean())
        .pivot(values = "distance", columns = "rank", index = ["singlem_package", "domain"], aggregate_function=None)
        .select(
            "singlem_package", "domain",
            phylum_cutoff = pl.sum_horizontal("d", "p") / 2,
            class_cutoff = pl.sum_horizontal("p", "c") / 2,
            order_cutoff = pl.sum_horizontal("c", "o") / 2,
            family_cutoff = pl.sum_horizontal("o", "f") / 2,
            genus_cutoff = pl.sum_horizontal("f", "g") / 2,
        )
        .sort("singlem_package", "domain")
    )
    distances.write_json(args.output, pretty=True)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
