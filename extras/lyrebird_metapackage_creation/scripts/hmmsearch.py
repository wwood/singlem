import os
import extern

import logging
import pathlib
from tqdm.contrib.concurrent import thread_map

hmm = snakemake.input.hmm
genome_proteins = [prot_filepath.strip('\n') for prot_filepath in open(snakemake.input.genome_proteins)]
output_dir = snakemake.params.output_dir
script_dir = snakemake.output.script_dir
num_threads = snakemake.threads

def process_a_chunk(param_set):
    i = param_set[0][2]
    script_file = script_dir + f"/hmmsearch_{i}.sh"
    if os.path.exists(f'{script_file}.done'):
        logging.info(f"Skipping {script_file} as it already exists")
        return
    with open(script_file, 'w+') as f:
        for params in param_set:
            genome_filepath, hmmsearch_output, __ = params
            f.write(f"hmmsearch -E 0.00001 --cpu 1 --tblout {hmmsearch_output} {hmm} {genome_filepath} > /dev/null\n")
    cmd = f"mqsub -t 32 -m 16 --no-email --hours 4 --name hmmsearch_{i} --segregated-log-files -- 'cat {script_file} | parallel -j32' && touch {script_file}.done"
    extern.run(cmd)

pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
pathlib.Path(script_dir).mkdir(parents=True, exist_ok=True)

logging.basicConfig(filename=snakemake.log[0],level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

param_list = []
param_set = []
i = 0
for genome_filepath in genome_proteins:
    hmmsearch_output = os.path.join(output_dir, os.path.splitext(os.path.basename(genome_filepath))[0] + ".txt")
    param_set.append((genome_filepath, hmmsearch_output, i))
    i += 1
    if len(param_set) > len(genome_proteins) / num_threads:
        param_list.append(param_set)
        param_set = []
if len(param_set) > 0:
    param_list.append(param_set)

logging.info(os.path.basename(__file__) + f": Processing {len(param_list)} chunks with {num_threads} threads")
thread_map(process_a_chunk, param_list, max_workers=num_threads, chunksize=1)

logging.info('done')

with open(snakemake.output.touch, 'w') as _: pass