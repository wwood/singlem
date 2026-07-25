import os
import logging
import pathlib
import extern
from tqdm.contrib.concurrent import process_map
import tempfile
import subprocess
import argparse

def process_a_genome(params):
    genome_faa, output, tigrfams, log, compressed = params
    logging.debug("Processing genome: " + os.path.basename(genome_faa))

    pathlib.Path(os.path.dirname(output)).mkdir(parents=True, exist_ok=True)
    pathlib.Path(os.path.dirname(log)).mkdir(parents=True, exist_ok=True)

    if compressed:
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

    if compressed:
        genome_tmp.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run TIGRFAMs search on genomes")
    parser.add_argument('--tigrfams', required=True, help='Path to TIGRFAMs HMM file')
    parser.add_argument('--genome-ids-file', required=True, help='File with genome IDs (one per line)')
    parser.add_argument('--fasta-base', required=True, help='Base directory for FASTA files')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads')
    parser.add_argument('--compressed', action='store_true', help='Genomes are compressed')
    
    args = parser.parse_args()
    
    tigrfams = args.tigrfams
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
        log = os.path.join(output_dir, f"{genome}.log")
        param_sets.append((genome_fasta, output_tsv, tigrfams, log, compressed))

    process_map(process_a_genome, param_sets, max_workers=num_threads, chunksize=1)

    logging.info('done')
