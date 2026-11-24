"""Map supplemented genome IDs to their FASTA paths.

The ID is derived from the path the same way 'singlem supplement' derived it:
the basename with its last extension removed (biolib_lite remove_extension), so
'.../AMXMAG_0092.fa.gz' -> 'AMXMAG_0092.fa'.
"""

import logging
import os

import polars as pl

logging.basicConfig(
    filename=snakemake.log[0], level=logging.INFO,
    format='%(asctime)s %(levelname)s: %(message)s')


def genome_id(path):
    f = os.path.splitext(os.path.basename(path))[0]
    return f[:-1] if f.endswith('.') else f


ids = pl.read_csv(snakemake.input.ids, has_header=False, new_columns=['genome'])
paths = (
    pl.read_csv(snakemake.input.genome_fasta_paths, separator='\t', has_header=False,
                new_columns=['path'], columns=[0])
    .with_columns(genome=pl.col('path').map_elements(genome_id, return_dtype=pl.String))
)
logging.info("%d supplemented genome IDs, %d candidate FASTA paths", ids.height, paths.height)

duplicated = paths.filter(pl.col('genome').is_duplicated()).filter(pl.col('genome').is_in(ids['genome']))
if duplicated.height > 0:
    raise Exception("Supplemented genome IDs match more than one FASTA path, e.g. {}".format(
        duplicated.head(10).to_dicts()))

joined = ids.join(paths, on='genome', how='left')
missing = joined.filter(pl.col('path').is_null())
if missing.height > 0:
    raise Exception("{} supplemented genomes have no FASTA path, e.g. {}".format(
        missing.height, missing['genome'].head(10).to_list()))

absent = [p for p in joined['path'] if not os.path.exists(p)]
if absent:
    raise Exception("{} FASTA paths do not exist, e.g. {}".format(len(absent), absent[:10]))

joined.select('path').write_csv(snakemake.output.paths, include_header=False)
logging.info("Wrote %d genome paths", joined.height)
