import os
import logging
import pathlib
import extern
from tqdm.contrib.concurrent import process_map

def process_a_genome(params):
    matches_fasta, mfqe, output, hmms_and_names, taxfiles,  log = params
    logging.debug("Processing genome: " + genome)

    pathlib.Path(os.path.dirname(output)).mkdir(parents=True, exist_ok=True)
    pathlib.Path(os.path.dirname(log)).mkdir(parents=True, exist_ok=True)

    cmd = f"python scripts/transpose_hmms_with_sequences.py --input-fasta {matches_fasta} --taxonomy {' '.join(taxfiles)} --hmm-seq {mfqe} --hmm-spkg {hmms_and_names} --output {output} &> {log}"
    extern.run(cmd)

def async_lustre_cleanup(target_dir):
    if not os.path.exists(target_dir):
        return
    raise Exception("target_dir {} already exists, please delete".format(target_dir))

protein_filepaths = [prot_filepath.strip('\n') for prot_filepath in open(snakemake.params.protein_filepaths)]
matches_dir = snakemake.params.matches_dir
mfqe_dir = snakemake.params.mfqe_dir
taxfiles = snakemake.input.taxfiles
output_dir = snakemake.params.output_dir

hmms_and_names = snakemake.params.hmms_and_names

logs_dir = snakemake.params.logs_dir
num_threads = snakemake.threads

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

logging.info(os.path.basename(__file__) + ": Processing {} genomes with {} threads".format(len(protein_filepaths), num_threads))

# delete output directory if it exists, or you will have problems with appending duplicates to existing files
if os.path.exists(output_dir):
    async_lustre_cleanup(output_dir)

pathlib.Path(os.path.dirname(output_dir)).mkdir(parents=True, exist_ok=True)

param_sets = []
for filepath in protein_filepaths:
    genome = os.path.basename(filepath).rsplit(".", 1)[0]
    matches_fasta = os.path.join(mfqe_dir, genome + '.faa')
    fam = os.path.join(matches_dir, genome + ".tsv")
    output = os.path.join(output_dir, genome)
    log = os.path.join(logs_dir, f"{genome}_transpose.log")

    if not os.path.exists(fam):
        logging.warning(f"Missing fam file for {genome}")
        continue
    param_sets.append((matches_fasta, fam, output, hmms_and_names, taxfiles, log))

process_map(process_a_genome, param_sets, max_workers=num_threads, chunksize=1)

logging.info('done')

with open(snakemake.output[0], 'w') as _: pass