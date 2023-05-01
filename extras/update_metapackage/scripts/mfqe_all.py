
import os
import logging
import pathlib
import extern
from tqdm.contrib.concurrent import process_map

def process_a_genome(params):
    # pfam_search, tigrfam_search, hmms_and_names, output, log = params
    fasta, fam, output, log = params
    logging.debug("Processing genome: " + genome)

    pathlib.Path(os.path.dirname(output)).mkdir(parents=True, exist_ok=True)
    pathlib.Path(os.path.dirname(log)).mkdir(parents=True, exist_ok=True)

    cmd = "cut -f1 {} |mfqe ".format(fam) + \
        "--input-fasta {} ".format(fasta) + \
        "--sequence-name-lists /dev/stdin " + \
        "--output-fasta-files {} ".format(output) + \
        "--output-uncompressed " + \
        "&> {}".format(log)
    extern.run(cmd)


genomes = snakemake.params.genome_ids
fam_directory = snakemake.params.fam_directory
# pfam_search_directory = snakemake.params.pfam_search_directory
# tigrfam_search_directory = snakemake.params.tigrfam_search_directory
# hmms_and_names = snakemake.params.hmms_and_names
output_dir = snakemake.params.output_dir
logs_dir = snakemake.params.logs_dir
num_threads = snakemake.threads

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

logging.info(os.path.basename(__file__) + ": Processing {} genomes with {} threads".format(len(genomes), num_threads))

param_sets = []
for genome in genomes:
    fasta = os.path.join(snakemake.config['gtdb_protein_faa_reps'], genome)
    fam = os.path.join(fam_directory, genome + ".fam")
    output = os.path.join(output_dir, genome)
    log = os.path.join(logs_dir, f"{genome}_mfqe.log")
    # pfam_search, tigrfam_search, hmms_and_names, output, log
    param_sets.append((fasta, fam, output, log))

process_map(process_a_genome, param_sets, max_workers=num_threads, chunksize=1)

logging.info('done')

# touch snakemake.output[0]
with open(snakemake.output[0], 'w') as _: pass
