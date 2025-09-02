import os
import sys
import argparse
import logging

import polars as pl
import pandas as pd

input_dir = snakemake.params.input_dir
input_files = [f for f in os.listdir(input_dir)]
metapackage = snakemake.input.metapackage
included_ids = set()
id_to_tax = {}

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

# read in viral taxonomy file
with open(snakemake.input.viral_taxonomy) as r:
    for line in r:
        id_, tax = line.strip().split("\t")
        id_to_tax[id_] = tax

# get set of possible ids, weird bug where a taxon is present from the transcript search but not the metapackage readname store
# would need metapackage v6 to access directly from the readname store so we go through each taxonomy pickle instead
for spkg in os.listdir(metapackage):
    if spkg.endswith(".spkg"):
        logging.info(f"Reading {spkg}")
        tax_hash = pd.read_pickle(os.path.join(metapackage, spkg, "taxonomy_hash.pickle"))
        for key in tax_hash:
            included_ids.add(key.split('~')[1]) 

#for each input file, read in first line as header, change taxonomy column in all rows to viral taxonomy based on filename, concatenate to single output file
output = pl.DataFrame()
i = 0
j = 0
for f in input_files:
    if not f.endswith(".tsv"):
        continue
    id_ = f.rsplit(".", 2)[0]
    if id_ not in included_ids:
        j += 1
        with open(os.path.join(input_dir, f)) as r:
            if len(r.readlines()) > 1:
                logging.warning(f"Skipping {id_} as it is not in the metapackage")
                i += 1
        continue
    taxonomy = id_to_tax[id_]
    logging.info(f"Reading {f}")
    df = pl.read_csv(os.path.join(input_dir, f), separator="\t")
    df = df.with_columns([pl.lit(taxonomy).alias("taxonomy")])
    output = pl.concat([output, df])
logging.info(f"Skipped {j} files, of which {i} were not empty.")

logging.info(f"Writing output to {snakemake.output[0]}")
output.write_csv(snakemake.output[0], separator="\t")
