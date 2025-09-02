import os
from os.path import join
import logging
import pathlib
from tqdm.contrib.concurrent import thread_map
import extern
import tempfile
import atexit
from concurrent.futures import ThreadPoolExecutor, as_completed

transcript_files = snakemake.input.transcript_fastas
script_directory = snakemake.output.script_dir
output_directory = snakemake.output.dir
log_directory = snakemake.params.logs
num_threads = snakemake.threads

pathlib.Path(output_directory).mkdir(parents=True, exist_ok=True)
pathlib.Path(script_directory).mkdir(parents=True, exist_ok=True)
pathlib.Path(log_directory).mkdir(parents=True, exist_ok=True)

metapackage = snakemake.input.metapackage

def process_a_chunk(param_set):
    i = param_set[0][3]
    script_file = join(script_directory, f"lyrebird_transcripts_{i}.sh")
    with open(script_file, 'w+') as f:
        for params in param_set:
            input_path, output_path, log_path, __ = params
            f.write(f"singlem pipe -1 {input_path} --metapackage {metapackage} --otu-table {output_path} --no-assign-taxonomy &> {log_path}\n")
        logging.info(f"Finished writing script {script_file}")
    cmd = f"mqsub -t 16 -m 16 --no-email --hours 4 --name lyrebird_transcripts_{i} --segregated-log-files -- 'cat {script_file} | parallel -j16 --will-cite'"
    extern.run(cmd)

logging.basicConfig(
    filename=os.path.join(snakemake.log[0]),
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
)
logging.info("Running lyrebird pipe to generate transcript OTUs from input transcript list: {}".format(transcript_files))

param_list = []
param_set = []
i = 0

for input_file in [f for f in open(transcript_files, 'r').readlines()]:
    input_path = input_file.strip()
    if input_path.endswith('.fna'):
        output_path = join(output_directory, os.path.basename(input_path).replace('.fna', '.otu_table.tsv'))
        log_path = join(log_directory, f"{os.path.basename(input_path).replace('.fna', '.log')}")
        i += 1
        if len(param_set) >= 1000:
            param_list.append(param_set)
            param_set = []
        param_set.append((input_path, output_path, log_path, i))
if len(param_set) > 0:
    param_list.append(param_set)
logging.info(f"Processing {len(param_list)} chunks with {num_threads} threads")

with ThreadPoolExecutor(max_workers=num_threads) as executor:
    futures = [executor.submit(process_a_chunk, chunk) for chunk in param_list]
    try:
        for future in as_completed(futures):
            # re-raise any exception that happened in worker
            future.result()
    except Exception:
        # try to cancel remaining work and propagate the error to main process
        for f in futures:
            try:
                f.cancel()
            except Exception:
                pass
        raise

logging.info("Finished processing lyrebird transcripts.")