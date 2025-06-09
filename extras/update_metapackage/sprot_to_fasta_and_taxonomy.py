#!/usr/bin/env python

import sys, gzip, argparse
from Bio import SwissProt
 
parser = argparse.ArgumentParser()
parser.add_argument('-s','--uniprot_sprot_dat_gz', help='path to uniprot_sprot.dat.gz', required=True)
parser.add_argument('-f', help='output fasta file', required=True)
parser.add_argument('-t', help='output taxonomy file', required=True)
args = parser.parse_args()

kingdom_counts = {}
with open(args.f,'w') as output_fasta_fh:
    with open(args.t,'w') as output_taxonomy_fh:
        for record in SwissProt.parse(gzip.open(args.uniprot_sprot_dat_gz)):
            tax = record.organism_classification[0]
            try:
                kingdom_counts[tax] += 1
            except KeyError:
                kingdom_counts[tax] = 1
            if tax == 'Eukaryota':
                new_classification = "d__%s; %s" % (
                    record.organism_classification[0],
                    record.organism_classification[1])
                output_taxonomy_fh.write(
                    "\t".join([record.entry_name, new_classification])+
                    "\n")
                output_fasta_fh.write(">%s\n" % record.entry_name)
                output_fasta_fh.write("%s\n" % record.sequence)
 
for kingdom, count in kingdom_counts.items():
    sys.stderr.write("\t".join([kingdom,str(count)])+"\n")
