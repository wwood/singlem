
# Get around reuqirement for GTDBTK_DATA_PATH
import os
os.environ['GTDBTK_DATA_PATH'] = '/work/microbiome/db/gtdb/gtdb_release207_v2'

from gtdbtk.external.pypfam.Scan.PfamScan import PfamScan
import logging
import pathlib
from tqdm.contrib.concurrent import process_map
import tempfile
import subprocess

def process_a_genome(params):
    genome_input, output_tsv, pfams = params
    logging.debug("Processing genome: " + genome)
    pathlib.Path(os.path.dirname(output_tsv)).mkdir(parents=True, exist_ok=True)

    if 'compressed_genome_data' in snakemake.config and snakemake.config['compressed_genome_data']:
        # Temporary file to store the uncompressed genome
        with tempfile.NamedTemporaryFile(prefix='singlem-pfam-scan-', suffix='.faa') as genome_tmp:
            subprocess.check_call(['bash','-c', 'zcat {} > {}'.format(genome_input, genome_tmp.name)])
            pfam_scan = PfamScan(cpu=1, fasta=genome_tmp.name, dir=pfams)
            pfam_scan.search()
            pfam_scan.write_results(output_tsv, None, None, None, None)
    else:
        pfam_scan = PfamScan(cpu=1, fasta=genome_input, dir=pfams)
        pfam_scan.search()
        pfam_scan.write_results(output_tsv, None, None, None, None)


pfams = snakemake.params.pfams
genomes = snakemake.params.genome_ids
fasta_base = snakemake.config['gtdb_protein_faa_reps']
output_dir = snakemake.params.output_dir
num_threads = snakemake.threads

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

logging.info(os.path.basename(__file__) + ": Processing {} genomes with {} threads".format(len(genomes), num_threads))

param_sets = []
for genome in genomes:
    genome_fasta = os.path.join(fasta_base, genome)
    output_tsv = os.path.join(output_dir, f"hmmsearch/pfam/{genome}.tsv")
    param_sets.append((genome_fasta, output_tsv, pfams))

process_map(process_a_genome, param_sets, max_workers=num_threads, chunksize=1)

logging.info('done')

# touch snakemake.output[0]
with open(snakemake.output[0], 'w') as _: pass
