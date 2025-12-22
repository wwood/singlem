
# Get around reuqirement for GTDBTK_DATA_PATH
import os
os.environ['GTDBTK_DATA_PATH'] = '/work/microbiome/db/gtdb/gtdb_release207_v2'

from gtdbtk.external.pypfam.Scan.PfamScan import PfamScan
import logging
import pathlib
from tqdm.contrib.concurrent import process_map
import tempfile
import subprocess
import argparse

def process_a_genome(params):
    genome_input, output_tsv, pfams, compressed = params
    logging.debug("Processing genome: " + os.path.basename(genome_input))
    pathlib.Path(os.path.dirname(output_tsv)).mkdir(parents=True, exist_ok=True)

    if compressed:
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Pfam search on genomes")
    parser.add_argument('--pfams', required=True, help='Path to Pfam database directory')
    parser.add_argument('--genome-ids-file', required=True, help='File with genome IDs (one per line)')
    parser.add_argument('--fasta-base', required=True, help='Base directory for FASTA files')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--compressed', action='store_true', help='Genomes are compressed')
    
    args = parser.parse_args()
    
    pfams = args.pfams
    with open(args.genome_ids_file, 'r') as f:
        genomes = [line.strip() for line in f if line.strip()]
    fasta_base = args.fasta_base
    output_dir = args.output_dir
    num_threads = args.threads
    compressed = args.compressed

    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

    logging.info(os.path.basename(__file__) + ": Processing {} genomes with {} threads".format(len(genomes), num_threads))

    param_sets = []
    for genome in genomes:
        genome_fasta = os.path.join(fasta_base, genome)
        output_tsv = os.path.join(output_dir, f"{genome}.tsv")
        param_sets.append((genome_fasta, output_tsv, pfams, compressed))

    process_map(process_a_genome, param_sets, max_workers=num_threads, chunksize=1)

    logging.info('done')
