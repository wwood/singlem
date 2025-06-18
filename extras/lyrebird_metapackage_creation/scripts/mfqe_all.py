import os
import logging
import pathlib
from tqdm.contrib.concurrent import process_map

def process_a_genome(params):
    fasta, matches, output, log = params
    logging.debug(f"Processing {fasta}")

    pathlib.Path(os.path.dirname(output)).mkdir(parents=True, exist_ok=True)
    pathlib.Path(os.path.dirname(log)).mkdir(parents=True, exist_ok=True)

    cmd = f"cut -f1 {matches} | mfqe  --input-fasta {fasta} --sequence-name-lists /dev/stdin --output-fasta-files {output} --output-uncompressed &> {log}"
    os.system(cmd)

protein_filepaths = snakemake.params.protein_filepaths
match_directory = snakemake.params.match_directory
output_dir = snakemake.params.output_dir
num_threads = snakemake.threads
logs_dir = snakemake.params.logs_dir

logging.basicConfig(filename=snakemake.log[0],level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

count = 0
basename_to_filepath = {}
with open(protein_filepaths) as r:
    for line in r:
        filepath = line.strip()
        basename = filepath.split('/')[-1].rsplit('.', 1)[0]
        basename_to_filepath[basename] = filepath
        count += 1

logging.info(os.path.basename(__file__) + ": Processing {} genome protein files".format(count))

param_sets = []
for match in os.listdir(match_directory):
    fasta = basename_to_filepath[os.path.splitext(match)[0]]
    matches = os.path.join(match_directory, match)
    output = os.path.join(output_dir, os.path.splitext(match)[0] + ".faa")
    log = os.path.join(logs_dir, os.path.splitext(match)[0] + ".log")
    param_sets.append((fasta, matches, output, log))
logging.info(os.path.basename(__file__)+ f": Processing {len(param_sets)} genomes with {num_threads} threads")
process_map(process_a_genome, param_sets, max_workers=num_threads, chunksize=1)

logging.info('done')

with open(snakemake.output[0], 'w') as _: pass