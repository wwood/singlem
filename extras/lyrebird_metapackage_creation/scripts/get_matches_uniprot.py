### UNUSED, dsDNA phages don't seem to have any eukaryotic hits but this may be useful for eukaryotic viruses in the future

import os
import logging
import pathlib
import csv
import Bio.SearchIO.HmmerIO as HmmerIO
from collections import Counter
from collections import defaultdict

hmms_and_names = snakemake.params.hmms_and_names
outfile = snakemake.params.outfile
logs_dir = snakemake.params.logs_dir
num_threads = snakemake.threads
hmmsearch_input = snakemake.params.hmmsearch_file

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

logging.info("Reading in HMM IDs...")
hmms_to_vogs = {}
with open(hmms_and_names, 'r') as hmm_list_file:
    for line in hmm_list_file.readlines()[1:]:
        vog, hmm = line.split()[:2]
        hmms_to_vogs[hmm] = vog

logging.info("Creating match list...")
match_list = []
hmm_hit_scores = defaultdict(list)
hmm_count = 0
hit_count = 0
with open(hmmsearch_input, 'r') as hmmsearch_file:
    for qresult in HmmerIO.hmmer3_tab.Hmmer3TabParser(hmmsearch_file):
        hmm = hmms_to_vogs[qresult.id]
        hmm_count += 1
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
logging.info("Removing HMMs with multiple gene hits...")
hmm_counter = Counter([hmm_id for __, hmm_id in derep_list])
output_list = [[seq_id, hmm_id] for seq_id, hmm_id in derep_list if hmm_counter[hmm_id] == 1]

output_dict = defaultdict(list)
for seq_id, hmm_id in output_list:
    output_dict[seq_id].append(hmm_id)
output_list = [[seq_id, ','.join(hmm_ids)] for seq_id, hmm_ids in output_dict.items()]

pathlib.Path(os.path.dirname(outfile)).mkdir(parents=True, exist_ok=True)
with open(outfile, 'w+') as output_file:
    output_writer = csv.writer(output_file, delimiter='\t')
    output_writer.writerows(output_list)
logging.info("Done.")

# touch snakemake.output[0]
with open(snakemake.output[0], 'w') as _: pass