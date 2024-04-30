
import os
import logging
import pathlib
import extern
from tqdm.contrib.concurrent import process_map
import tempfile
import subprocess

def process_a_genome(params):
    genome_faa, output, tigrfams, log = params
    logging.debug("Processing genome: " + genome)

    pathlib.Path(os.path.dirname(output_tsv)).mkdir(parents=True, exist_ok=True)
    pathlib.Path(os.path.dirname(log)).mkdir(parents=True, exist_ok=True)

    genome_is_compressed = 'compressed_genome_data' in snakemake.config and snakemake.config['compressed_genome_data']
    if genome_is_compressed:
        # Temporary file to store the uncompressed genome
        original_genome_faa = genome_faa
        genome_tmp = tempfile.NamedTemporaryFile(prefix='singlem-tigrfam-scan-', suffix='.faa')
        genome_faa = genome_tmp.name
        subprocess.check_call(['bash','-c', 'zcat {} > {}'.format(original_genome_faa, genome_tmp.name)])

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

    if genome_is_compressed:
        genome_tmp.close()

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
