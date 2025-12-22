import os
import logging
import pathlib
import extern
import argparse
from tqdm.contrib.concurrent import process_map

def process_a_genome(params):
    fasta, fam, output, log, genome_is_compressed = params
    logging.debug("Processing genome: " + os.path.basename(fasta))

    pathlib.Path(os.path.dirname(output)).mkdir(parents=True, exist_ok=True)
    pathlib.Path(os.path.dirname(log)).mkdir(parents=True, exist_ok=True)

    if genome_is_compressed:
        input_fasta = '<(zcat {})'.format(fasta)
    else:
        input_fasta = fasta

    cmd = "cut -f1 {} |mfqe ".format(fam) + \
        "--input-fasta {} ".format(input_fasta) + \
        "--sequence-name-lists /dev/stdin " + \
        "--output-fasta-files {} ".format(output) + \
        "--output-uncompressed " + \
        "&> {}".format(log)
    extern.run(cmd)


def main():
    parser = argparse.ArgumentParser(description='Process genomes with mfqe')
    parser.add_argument('--genome-ids-file', required=True, help='File containing list of genome IDs to process')
    parser.add_argument('--fam-directory', required=True, help='Directory containing .fam files')
    parser.add_argument('--gtdb-protein-faa-reps', required=True, help='Directory containing FASTA files')
    parser.add_argument('--output-dir', required=True, help='Output directory for processed files')
    parser.add_argument('--logs-dir', required=True, help='Directory for log files')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
    parser.add_argument('--compressed-genome-data', action='store_true', help='Input genomes are compressed')
    
    args = parser.parse_args()

    with open(args.genome_ids_file, 'r') as f:
        genomes = [line.strip() for line in f.readlines()]
    fam_directory = args.fam_directory
    gtdb_protein_faa_reps = args.gtdb_protein_faa_reps
    output_dir = args.output_dir
    logs_dir = args.logs_dir
    num_threads = args.threads
    genome_is_compressed = args.compressed_genome_data

    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

    logging.info(os.path.basename(__file__) + ": Processing {} genomes with {} threads".format(len(genomes), num_threads))

    param_sets = []
    for genome in genomes:
        fasta = os.path.join(gtdb_protein_faa_reps, genome)
        fam = os.path.join(fam_directory, genome + ".fam")
        output = os.path.join(output_dir, genome)
        log = os.path.join(logs_dir, f"{genome}_mfqe.log")
        param_sets.append((fasta, fam, output, log, genome_is_compressed))

    process_map(process_a_genome, param_sets, max_workers=num_threads, chunksize=1)

    logging.info('done')


if __name__ == '__main__':
    main()

"""
Example usage in Snakemake shell directive:

shell:
    '''
    python {input.script} \\
        --genome-ids {params.genome_ids} \\
        --fam-directory {params.fam_directory} \\
        --gtdb-protein-faa-reps {config[gtdb_protein_faa_reps]} \\
        --output-dir {params.output_dir} \\
        --logs-dir {params.logs_dir} \\
        --output-file {output[0]} \\
        --threads {threads} \\
        {params.compressed_flag}
    '''

Where compressed_flag would be set in params as:
    compressed_flag = "--compressed-genome-data" if config.get('compressed_genome_data', False) else ""
"""
