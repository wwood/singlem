import os
import extern
import logging
import pathlib
from tqdm.contrib.concurrent import process_map

def process_a_genome(params):
    hmmsearch_path, output_tsv, log = params
    logging.debug(f"Processing {hmmsearch_path}")
    
    pathlib.Path(os.path.dirname(output_tsv)).mkdir(parents=True, exist_ok=True)
    pathlib.Path(os.path.dirname(log)).mkdir(parents=True, exist_ok=True)

    cmd = f"python scripts/get_matches_no_dup.py --hmmsearch-file {hmmsearch_path} --hmm-list {hmms_and_names} --output {output_tsv} &> {log}"
    extern.run(cmd)

def process_a_microbe(params):
    hmmsearch_path, output_tsv, log, proviruses = params
    logging.debug(f"Processing {hmmsearch_path}")

    pathlib.Path(os.path.dirname(output_tsv)).mkdir(parents=True, exist_ok=True)
    pathlib.Path(os.path.dirname(log)).mkdir(parents=True, exist_ok=True)

    cmd = f"python scripts/get_matches_microbial.py --hmmsearch-file {hmmsearch_path} --hmm-list {hmms_and_names} --output {output_tsv} --genomad-db {proviruses} &> {log}"
    extern.run(cmd)

hmmsearch_directory = snakemake.params.hmmsearch_directory
hmms_and_names = snakemake.params.hmms_and_names
output_dir = snakemake.params.output_dir
logs_dir = snakemake.params.logs_dir
num_threads = snakemake.threads
proviruses = snakemake.params.proviruses

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')
pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

param_sets = []

for hmmsearch_file in os.listdir(hmmsearch_directory):
    hmmsearch_path = os.path.join(hmmsearch_directory, hmmsearch_file)
    output_tsv = os.path.join(output_dir, hmmsearch_file.replace(".txt", ".tsv"))
    log = os.path.join(logs_dir, hmmsearch_file.replace(".txt", "_matching.log"))
    if proviruses:
        param_sets.append((hmmsearch_path, output_tsv, log, proviruses))
    else:
        param_sets.append((hmmsearch_path, output_tsv, log))

logging.info(os.path.basename(__file__) + ": Processing {} genomes with {} threads".format(len(param_sets), num_threads))
if proviruses:
    logging.info("Using provirus db to filter out provirus hits")
    process_map(process_a_microbe, param_sets, max_workers=num_threads, chunksize=1)
else:
    process_map(process_a_genome, param_sets, max_workers=num_threads, chunksize=1)

logging.info('done')

# touch snakemake.output[0]
with open(snakemake.output[0], 'w') as _: pass