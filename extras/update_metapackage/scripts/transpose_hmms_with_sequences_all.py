import os
import logging
import pathlib
import extern
import argparse
from tqdm.contrib.concurrent import process_map


def process_a_genome(params):
    matches_fasta, fam, output, hmms_and_names, gtdb_bac_tax, gtdb_arc_tax, log = params
    logging.debug("Processing genome: " + os.path.basename(matches_fasta))

    pathlib.Path(os.path.dirname(output)).mkdir(parents=True, exist_ok=True)
    pathlib.Path(os.path.dirname(log)).mkdir(parents=True, exist_ok=True)

    cmd = "python scripts/transpose_hmms_with_sequences.py " + \
        "--input-fasta {} ".format(matches_fasta) + \
        "--bacterial-taxonomy {} ".format(gtdb_bac_tax) + \
        "--archaeal-taxonomy {} ".format(gtdb_arc_tax) + \
        "--hmm-seq {} ".format(fam) + \
        "--hmm-spkg {} ".format(hmms_and_names) + \
        "--output {} ".format(output) + \
        "&> {} ".format(log)
    extern.run(cmd)


def main():
    parser = argparse.ArgumentParser(description='Transpose HMMs with sequences for all genomes')
    parser.add_argument('--genome-ids-list', required=True, help='File containing list of genome IDs to process')
    parser.add_argument('--matches-dir', required=True, help='Directory containing matches files')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--gtdb-bac-tax', required=True, help='GTDB bacterial taxonomy file')
    parser.add_argument('--gtdb-arc-tax', required=True, help='GTDB archaeal taxonomy file')
    parser.add_argument('--hmms-and-names', required=True, help='HMMs and names file')
    parser.add_argument('--logs-dir', required=True, help='Directory for log files')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
    
    args = parser.parse_args()

    with open(args.genome_ids_list, 'r') as f:
        genomes = [line.strip() for line in f if line.strip()]
    matches_dir = args.matches_dir
    output_dir = args.output_dir
    gtdb_bac_tax = args.gtdb_bac_tax
    gtdb_arc_tax = args.gtdb_arc_tax
    hmms_and_names = args.hmms_and_names
    logs_dir = args.logs_dir
    num_threads = args.threads

    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

    logging.info(os.path.basename(__file__) + ": Processing {} genomes with {} threads".format(len(genomes), num_threads))

    pathlib.Path(os.path.dirname(output_dir)).mkdir(parents=True, exist_ok=True)

    param_sets = []
    for genome in genomes:
        matches_fasta = os.path.join(matches_dir, genome)
        fam = os.path.join(matches_dir, genome + ".fam")
        output = os.path.join(output_dir, genome)
        log = os.path.join(logs_dir, f"{genome}_transpose.log")
        param_sets.append((matches_fasta, fam, output, hmms_and_names, gtdb_bac_tax, gtdb_arc_tax, log))

    process_map(process_a_genome, param_sets, max_workers=num_threads, chunksize=1)

    logging.info('done')


if __name__ == '__main__':
    main()

