################################
### get_matches_microbial.py ###
################################

# Search hmmsearch output for matching HMM IDs, skip provirus hits

import argparse
import csv
import Bio.SearchIO.HmmerIO as HmmerIO
import logging
from collections import Counter
from collections import defaultdict

parser = argparse.ArgumentParser(description='Search HMMER output for matching HMM IDs.')
parser.add_argument('--hmmsearch-file', type=str, metavar='<HMMSEARCH TBLOUT>', help='path to pfam output file')
parser.add_argument('--hmm-list', type=str, metavar='<HMMS>', help='path HMM list')
parser.add_argument('--output', type=str, metavar='<OUTPUT>', help='path to output file')
parser.add_argument('--genomad-db', type=str, metavar='<GENOMAD DB>', help='path to genomad database', required=False)
parser.add_argument('--log', type=str, metavar='<LOG>', help='path to log file')

logging.basicConfig(
    filename=getattr(args, 'log'),
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
)

args = parser.parse_args()
hmms_and_names = getattr(args, 'hmm_list')
hmmsearch_input = getattr(args, 'hmmsearch_file')
output_path = getattr(args, 'output')
genomad_db = getattr(args, 'genomad_db')

logging.info("Reading in provirus db...")
proviruses = {}
contamination = set()
if genomad_db:
    with open(genomad_db, 'r') as genomad_file:
        for line in csv.DictReader(genomad_file, delimiter='\t'):
            virus_name = line["seq_name"]
            if line["topology"] == "Provirus":
                virus_name = virus_name.split('|')[0]
                coordinates = tuple(map(int, line["coordinates"].split('-')))
                if virus_name not in proviruses:
                    proviruses[virus_name] = []
                proviruses[virus_name].append(coordinates)
            else:
                contamination.add(virus_name)

logging.info("Found {} proviruses in db".format(len(proviruses)))
logging.info("Found {} viral contamination sequences in db".format(len(contamination)))

# Read in HMM IDs
logging.info("Reading in HMM IDs...")
hmms_to_phrogs = {}
with open(hmms_and_names, 'r') as hmm_list_file:
    for line in hmm_list_file.readlines()[1:]:
        phrog, hmm = line.split()[:2]
        hmms_to_phrogs[hmm] = phrog
logging.info("Creating match list...")
match_list = []
hmm_hit_scores = defaultdict(list)
hmm_count = 0
hit_count = 0
hits_in_provirus = 0
with open(hmmsearch_input, 'r') as hmmsearch_file:
    for qresult in HmmerIO.hmmer3_tab.Hmmer3TabParser(hmmsearch_file):
        hmm = hmms_to_phrogs[qresult.id]
        hmm_count += 1
        for hit in qresult.hits:
            genome = hit.id.rsplit('_', 1)[0]
            if genome in contamination:
                hits_in_provirus += 1
                continue
            if genome in proviruses:
                start = int(hit.description.replace('# ', '').split()[0])
                end = int(hit.description.replace('# ', '').split()[1])
                in_provirus = False
                # skip hmmsearch hit if it overlaps with provirus
                for provirus in proviruses[genome]:
                    provirus_start = provirus[0]
                    provirus_end = provirus[1]
                    if start <= provirus_end and end >= provirus_start:
                        in_provirus = True
                        hits_in_provirus += 1
                        break
                if not in_provirus:
                    match_list.append((hit.id, hmm, hit.bitscore))
                    hmm_hit_scores[hit.id].append(hit.bitscore)
                    hit_count += 1
            else:
                match_list.append((hit.id, hmm, hit.bitscore))
                hmm_hit_scores[hit.id].append(hit.bitscore)
                hit_count += 1
logging.info("Found {} hits for {} HMMs, skipping {} hits found in proviruses db".format(hit_count, hmm_count, hits_in_provirus))

#### We actually want to keep multiple gene hits for decoy sequences since we want the decoys no matter what, no dereplication
output_list = [[seq_id, hmm_id] for seq_id, hmm_id, __ in match_list]

output_dict = defaultdict(list)
for seq_id, hmm_id in output_list:
    output_dict[seq_id].append(hmm_id)
output_list = [[seq_id, ','.join(hmm_ids)] for seq_id, hmm_ids in output_dict.items()]

with open(output_path, 'w') as output_file:
    output_writer = csv.writer(output_file, delimiter='\t')
    output_writer.writerows(output_list)
logging.info("Done.")