#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict

from singlem.metapackage import Metapackage
from singlem.weebill import WeebillProfiler


parser = argparse.ArgumentParser()
parser.add_argument("--profile", required=True)
parser.add_argument("--metapackage", required=True)
parser.add_argument("--sample", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()

profiler = WeebillProfiler()
with open(args.profile) as handle:
    rows = list(csv.DictReader(handle, delimiter="\t"))
accessions = {profiler._extract_accession(row["Genome_file"]) for row in rows}
accessions.discard(None)
taxonomy = Metapackage.acquire(args.metapackage).genome_accession_to_taxonomy(accessions)

coverage_by_taxonomy = defaultdict(float)
for row in rows:
    accession = profiler._extract_accession(row["Genome_file"])
    if accession in taxonomy:
        coverage_by_taxonomy[taxonomy[accession]] += float(row["Eff_cov"])

if not coverage_by_taxonomy:
    raise SystemExit("No weebill genomes could be mapped to r207 taxonomy")

with open(args.output, "w") as handle:
    handle.write("sample\tcoverage\ttaxonomy\n")
    for taxon, coverage in sorted(coverage_by_taxonomy.items()):
        handle.write(f"{args.sample}\t{coverage}\t{taxon}\n")
