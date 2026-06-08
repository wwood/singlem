#!/usr/bin/env python3
"""Annotate a native sylph profile TSV with GTDB taxonomy.

`sylph profile` reports a `Genome_file` column (the genome in the sylph
database) and an `Eff_cov` column, but not a taxonomy string. SingleM's
`condense --sylph-profile` expects a pre-annotated TSV carrying a GTDB taxonomy
column. This script bridges the two: it extracts the genome accession
(GCA_/GCF_...) from `Genome_file`, looks it up in one or more GTDB taxonomy TSVs
(`accession<TAB>taxonomy`, with the usual GB_/RS_ accession prefixes), and writes
a TSV with columns `Sample_file`, `taxonomy`, `Eff_cov`.
"""

import argparse
import csv
import re
import sys

ACCESSION_RE = re.compile(r'(GC[AF]_\d+\.\d+)')


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sylph-tsv', required=True, help='native sylph profile output TSV')
    parser.add_argument('--taxonomy', nargs='+', required=True,
        help='one or more GTDB taxonomy TSVs (accession<TAB>taxonomy)')
    parser.add_argument('--output', required=True, help='annotated TSV to write')
    return parser.parse_args()


def extract_accession(genome_file):
    match = ACCESSION_RE.search(genome_file)
    return match.group(1) if match else None


def main():
    args = parse_args()

    # First pass: read the sylph rows and collect the accessions we need so the
    # (large) taxonomy files can be streamed once, keeping only relevant lines.
    rows = []
    needed = set()
    with open(args.sylph_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        if reader.fieldnames is None or 'Genome_file' not in reader.fieldnames or 'Eff_cov' not in reader.fieldnames:
            sys.exit("sylph TSV must contain Genome_file and Eff_cov columns; found {}".format(reader.fieldnames))
        for row in reader:
            accession = extract_accession(row['Genome_file'])
            rows.append((row.get('Sample_file', ''), accession, row['Eff_cov'], row['Genome_file']))
            if accession is not None:
                needed.add(accession)

    accession_to_taxonomy = {}
    for taxonomy_path in args.taxonomy:
        with open(taxonomy_path) as f:
            for line in f:
                parts = line.rstrip('\n').split('\t')
                if len(parts) < 2:
                    continue
                accession = parts[0]
                if accession.startswith('GB_') or accession.startswith('RS_'):
                    accession = accession[3:]
                if accession in needed:
                    accession_to_taxonomy[accession] = parts[1]

    num_written = 0
    num_skipped = 0
    with open(args.output, 'w') as out:
        out.write('Sample_file\ttaxonomy\tEff_cov\n')
        for sample, accession, eff_cov, genome_file in rows:
            taxonomy = accession_to_taxonomy.get(accession)
            if taxonomy is None:
                sys.stderr.write("No GTDB taxonomy for {} (accession {}), skipping\n".format(genome_file, accession))
                num_skipped += 1
                continue
            out.write('\t'.join([sample, taxonomy, eff_cov]) + '\n')
            num_written += 1
    sys.stderr.write("Wrote {} annotated sylph rows, skipped {}\n".format(num_written, num_skipped))


if __name__ == '__main__':
    main()
