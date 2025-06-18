#############################
### get_matches_no_dup.py ###
#############################

# Search hmmsearch output for matching HMM IDs

import argparse
import csv
import re
import os
import Bio.SearchIO.HmmerIO as HmmerIO
import logging
from collections import Counter
from collections import defaultdict

parser = argparse.ArgumentParser(description='Search HMMER output for matching HMM IDs.')
parser.add_argument('--hmmsearch-file', type=str, metavar='<HMMSEARCH TBLOUT>', help='path to pfam output file')
parser.add_argument('--hmm-list', type=str, metavar='<HMMS>', help='path HMM list')
parser.add_argument('--output', type=str, metavar='<OUTPUT>', help='path to output file')

args = parser.parse_args()
HMM_id_list = getattr(args, 'hmm_list')
hmmsearch_input = getattr(args, 'hmmsearch_file')
output_path = getattr(args, 'output')

# Read in HMM IDs
logging.info("Reading in HMM IDs...")
hmms_to_phrogs = {}
with open(HMM_id_list, 'r') as hmm_list_file:
    for line in hmm_list_file.readlines()[1:]:
        phrog, hmm = line.split()[:2]
        hmms_to_phrogs[hmm] = phrog
logging.info("Creating match list...")
match_list = []
hmm_hit_scores = defaultdict(list)
hmm_count = 0
hit_count = 0
with open(hmmsearch_input, 'r') as hmmsearch_file:
    for qresult in HmmerIO.hmmer3_tab.Hmmer3TabParser(hmmsearch_file):
        hmm_count += 1
        hmm = hmms_to_phrogs[qresult.id]
        for hit in qresult.hits:
            match_list.append((hit.id, hmm, hit.bitscore))
            hmm_hit_scores[hit.id].append(hit.bitscore)
            hit_count += 1
logging.info("Found {} hits for {} HMMs".format(hit_count, hmm_count))

gene_counter = Counter([seq_id for seq_id,__,__ in match_list])
logging.info("Removing duplicate hits...")
derep_list = [[seq_id, hmm_id] for seq_id, hmm_id, score in match_list if 
                gene_counter[seq_id] == 1 or score == max(hmm_hit_scores[seq_id])]

# Remove HMMs with multiple gene hits
logging.info("Removing sequences with multiple gene hits...")
hmm_counter = Counter([hmm_id for __, hmm_id in derep_list])
output_list = [[seq_id, hmm_id] for seq_id, hmm_id in derep_list if hmm_counter[hmm_id] == 1]

output_dict = defaultdict(list)
for seq_id, hmm_id in output_list:
    output_dict[seq_id].append(hmm_id)
output_list = [[seq_id, ','.join(hmm_ids)] for seq_id, hmm_ids in output_dict.items()]

with open(output_path, 'w') as output_file:
    output_writer = csv.writer(output_file, delimiter='\t')
    output_writer.writerows(output_list)
logging.info("Done.")