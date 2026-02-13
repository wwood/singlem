import os
import logging
import pathlib
import extern
import argparse
from tqdm.contrib.concurrent import process_map

def process_a_genome(params):
    pfam_search, tigrfam_search, hmms_and_names, output, log = params
    logging.debug("Processing genome: " + os.path.basename(pfam_search).replace('.tsv', ''))

    pathlib.Path(os.path.dirname(output)).mkdir(parents=True, exist_ok=True)
    pathlib.Path(os.path.dirname(log)).mkdir(parents=True, exist_ok=True)

    cmd = "python scripts/get_matches_no_dup.py " \
        "--pfam-search {} ".format(pfam_search) + \
        "--tigrfam-search {} ".format(tigrfam_search) + \
        "--hmm-list {} ".format(hmms_and_names) + \
        "--output {} ".format(output) + \
        "&> {}".format(log)
    extern.run(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process genomes for HMM matches')
    parser.add_argument('--genome-ids-file', required=True, help='File containing newline-separated genome IDs')
    parser.add_argument('--pfam-search-directory', required=True, help='Directory containing Pfam search results')
    parser.add_argument('--tigrfam-search-directory', required=True, help='Directory containing TIGRFam search results')
    parser.add_argument('--hmms-and-names', required=True, help='HMM list file')
    parser.add_argument('--output-dir', required=True, help='Output directory for results')
    parser.add_argument('--logs-dir', required=True, help='Directory for log files')
    parser.add_argument('--threads', type=int, default=1, help='Number of threads to use')
    
    args = parser.parse_args()
    
    with open(args.genome_ids_file, 'r') as f:
        genomes = [line.strip() for line in f if line.strip()]
    pfam_search_directory = args.pfam_search_directory
    tigrfam_search_directory = args.tigrfam_search_directory
    hmms_and_names = args.hmms_and_names
    output_dir = args.output_dir
    logs_dir = args.logs_dir
    num_threads = args.threads

    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

    logging.info(os.path.basename(__file__) + ": Processing {} genomes with {} threads".format(len(genomes), num_threads))

    param_sets = []
    for genome in genomes:
        pfam_search = os.path.join(pfam_search_directory, f"{genome}.tsv")
        tigrfam_search = os.path.join(tigrfam_search_directory, f"{genome}.tsv")
        output_tsv = os.path.join(output_dir, f"{genome}.fam")
        log = os.path.join(logs_dir, f"{genome}_matching.log")
        # pfam_search, tigrfam_search, hmms_and_names, output, log
        param_sets.append((pfam_search, tigrfam_search, hmms_and_names, output_tsv, log))

    process_map(process_a_genome, param_sets, max_workers=num_threads, chunksize=1)

    logging.info('done')

