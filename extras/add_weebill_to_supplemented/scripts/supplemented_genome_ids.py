"""List the genome IDs 'singlem supplement' added to a metapackage.

supplement names its reference sequences '<marker>~<genome>‡<gene>', where
<genome> is the genome FASTA basename with its last extension removed. Original
GTDB sequences ('<marker>~GB_GCA_xxx.1') have no '‡', so they are not listed.
"""

import json
import logging
import os
import re

logging.basicConfig(
    filename=snakemake.log[0], level=logging.INFO,
    format='%(asctime)s %(levelname)s: %(message)s')

# Must agree with singlem's _SUPPLEMENTED_GENOME_REGEX. The spkg prefix never
# contains '~', but the genome name may.
GENOME_REGEX = re.compile(r'^[^~]*~(.+?)‡')
VERSIONED_ACCESSION_REGEX = re.compile(r'GC[AF]_\d+\.\d+')

mpkg = snakemake.input.mpkg
with open(os.path.join(mpkg, 'CONTENTS.json')) as f:
    contents = json.load(f)
store = os.path.join(mpkg, contents['sqlite_db_path_key'])
logging.info("Reading read names from %s", store)

query = "SELECT read_name FROM read_name_taxonomy"
if store.endswith('.duckdb'):
    import duckdb
    con = duckdb.connect(store, read_only=True)
    cursor = con.execute(query)
else:
    import sqlite3
    con = sqlite3.connect('file:{}?mode=ro'.format(store), uri=True)
    cursor = con.execute(query)

genomes = set()
num_reads = 0
num_supplemented = 0
while True:
    rows = cursor.fetchmany(100000)
    if not rows:
        break
    for (read_name,) in rows:
        num_reads += 1
        match = GENOME_REGEX.search(read_name)
        if match:
            num_supplemented += 1
            genomes.add(match.group(1))
con.close()

logging.info("%d read names, %d from supplemented genomes, %d distinct supplemented genomes",
             num_reads, num_supplemented, len(genomes))
if not genomes:
    raise Exception("No supplemented genomes found in {}".format(store))

# SingleM matches weebill genomes to taxonomy by a versioned accession first, so
# a supplemented genome whose name carries one will only match if that accession
# is in the read-name store.
with_accession = sorted(g for g in genomes if VERSIONED_ACCESSION_REGEX.search(g))
if with_accession:
    logging.warning("%d supplemented genome IDs contain a versioned GCA_/GCF_ accession, e.g. %s",
                    len(with_accession), with_accession[:5])

with open(snakemake.output.ids, 'w') as out:
    for g in sorted(genomes):
        out.write(g + '\n')
