#############################
### get_matches_no_dup.py ###
#############################
# Author: Samuel Aroney
# Search pfam/tigrfam search output for matching HMM IDs

import argparse
import csv
import re
from collections import Counter
from collections import defaultdict


parser = argparse.ArgumentParser(description='Search pfam/tigrfam search output for matching HMM IDs.')
parser.add_argument('--pfam-search', type=str, metavar='<PFAM TBLOUT>', help='path to pfam output file')
parser.add_argument('--tigrfam-search', type=str, metavar='<TIGRFAM TBLOUT>', help='path to tigrfam output file')
parser.add_argument('--hmm-list', type=str, metavar='<REQ HMMS>', help='path to required HMM list')
parser.add_argument('--output', type=str, metavar='<OUTPUT>', help='path to fam output file')

args = parser.parse_args()
pfam_search_path = getattr(args, 'pfam_search')
tigrfam_search_path = getattr(args, 'tigrfam_search')
HMM_id_list = getattr(args, 'hmm_list')
NEW_VERSION_HMM_COLUMN = "HMM"
output_path = getattr(args, 'output')


PFAM_HMM_ID_COLUMN = 5
PFAM_HMM_HIT_COLUMN = 11
TIGRFAM_HMM_ID_COLUMN = 3
TIGRFAM_HMM_HIT_COLUMN = 5


with open(HMM_id_list) as hmm_file:
    hmms = csv.DictReader(hmm_file, delimiter="\t")
    HMM_set = set(line[NEW_VERSION_HMM_COLUMN] for line in hmms)


def get_match_from_file(fam_search_path, fam_type, match_list):
    with open(fam_search_path) as fam_file:
        input = csv.reader(fam_file, delimiter="\t")

        if fam_type.lower() == "pfam":
            FAM_HMM_ID_COLUMN = PFAM_HMM_ID_COLUMN
            FAM_HMM_HIT_COLUMN = PFAM_HMM_HIT_COLUMN
        elif fam_type.lower() == "tigrfam":
            FAM_HMM_ID_COLUMN = TIGRFAM_HMM_ID_COLUMN
            FAM_HMM_HIT_COLUMN = TIGRFAM_HMM_HIT_COLUMN

        for line in input:
            if len(line)>0 and not line[0].startswith("#"):
                line_split = re.split('\s{1,}', line[0])
                if line_split[FAM_HMM_ID_COLUMN] in HMM_set:
                    match_list.append([line_split[i] for i in [0, FAM_HMM_ID_COLUMN, FAM_HMM_HIT_COLUMN]])
    
    return match_list


with open(output_path, 'w') as output_file:
    output = csv.writer(output_file, delimiter="\t")

    match_list = []
    match_list = get_match_from_file(pfam_search_path, "pfam", match_list)
    match_list = get_match_from_file(tigrfam_search_path, "tigrfam", match_list)

    # For genes with multiple HMM hits, choose highest scoring hit
    hmm_hit_scores = defaultdict(list)
    for seq_id,_,score in match_list:
        hmm_hit_scores[seq_id].append(score)

    gene_counter = Counter([seq_id for seq_id,_,_ in match_list])

    derep_list = [[seq_id,hmm_id] for seq_id,hmm_id,score in match_list if
                    gene_counter[seq_id] == 1 or score == max(hmm_hit_scores[seq_id])]
    
    # Remove HMMs with multiple gene hits
    hmm_counter = Counter([hmm_id for _,hmm_id in derep_list])
    output_list = [[seq_id,hmm_id] for seq_id,hmm_id in derep_list if hmm_counter[hmm_id] == 1]


    output.writerows(output_list)

