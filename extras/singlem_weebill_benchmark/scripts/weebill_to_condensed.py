#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict

from singlem.metapackage import Metapackage
from singlem.metapackage_read_name_store import genome_key_from_genome_file


parser = argparse.ArgumentParser()
parser.add_argument("--profile", required=True)
parser.add_argument("--metapackage", required=True)
parser.add_argument("--sample", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

with open(args.profile) as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
genome_keys = {genome_key_from_genome_file(row["Genome_file"]) for row in rows}
taxonomy = Metapackage.acquire(args.metapackage).genome_key_to_taxonomy(genome_keys)

coverage_by_taxonomy = defaultdict(float)
for row in rows:
    genome_key = genome_key_from_genome_file(row["Genome_file"])
    if genome_key in taxonomy:
        coverage_by_taxonomy[taxonomy[genome_key]] += float(row["Eff_cov"])

if not coverage_by_taxonomy:
    raise SystemExit("No weebill genomes could be mapped to r207 taxonomy")

with open(args.output, "w") as handle:
    handle.write("sample\tcoverage\ttaxonomy\n")
    for taxon, coverage in sorted(coverage_by_taxonomy.items()):
        handle.write(f"{args.sample}\t{coverage}\t{taxon}\n")
