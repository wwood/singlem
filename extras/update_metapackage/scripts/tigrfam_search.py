
import os
import logging
import pathlib
import extern
from tqdm.contrib.concurrent import process_map

def process_a_genome(params):
    genome_faa, output, tigrfams, log = params
    logging.debug("Processing genome: " + genome)

    pathlib.Path(os.path.dirname(output_tsv)).mkdir(parents=True, exist_ok=True)
    pathlib.Path(os.path.dirname(log)).mkdir(parents=True, exist_ok=True)

    cmd = "hmmsearch " \
        "-o /dev/null " \
        "--tblout {} ".format(output) + \
        "--noali " \
        "--notextw " \
        "--cut_nc " \
        "--cpu 1 " \
        "{} ".format(tigrfams) + \
        "{} ".format(genome_faa) + \
        "&> {}".format(log)
    extern.run(cmd)

tigrfams = snakemake.config['tigrfams']
genomes = snakemake.params.genome_ids
fasta_base = snakemake.config['gtdb_protein_faa_reps']
output_dir = snakemake.params.output_dir
num_threads = snakemake.threads

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

logging.info(os.path.basename(__file__) + ": Processing {} genomes with {} threads".format(len(genomes), num_threads))

param_sets = []
for genome in genomes:
    genome_fasta = os.path.join(fasta_base, genome)
    output_tsv = os.path.join(output_dir, f"{genome}.tsv")
    log = os.path.join(output_dir, f"{genome}.log")
    param_sets.append((genome_fasta, output_tsv, tigrfams, log))

process_map(process_a_genome, param_sets, max_workers=num_threads, chunksize=1)

logging.info('done')

# touch snakemake.output[0]
with open(snakemake.output[0], 'w') as _: pass
