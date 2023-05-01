
import os
import logging
import pathlib
import extern
from tqdm.contrib.concurrent import process_map


def process_a_genome(params):
    matches_fasta, fam, output, hmms_and_names, gtdb_bac_tax, gtdb_arc_tax, log = params
    logging.debug("Processing genome: " + genome)

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


genomes = snakemake.params.genome_ids
# fam_directory = snakemake.params.fam_directory
matches_dir = snakemake.params.matches_dir
# pfam_search_directory = snakemake.params.pfam_search_directory
# tigrfam_search_directory = snakemake.params.tigrfam_search_directory
# hmms_and_names = snakemake.params.hmms_and_names
output_dir = snakemake.params.output_dir

gtdb_bac_tax = snakemake.config["gtdb_bac_tax"]
gtdb_arc_tax = snakemake.config["gtdb_arc_tax"]
hmms_and_names = snakemake.params.hmms_and_names

logs_dir = snakemake.params.logs_dir
num_threads = snakemake.threads

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

logging.info(os.path.basename(__file__) + ": Processing {} genomes with {} threads".format(len(genomes), num_threads))

pathlib.Path(os.path.dirname(output_dir)).mkdir(parents=True, exist_ok=True)

param_sets = []
for genome in genomes:
    # fasta = output_dir + "/hmmsearch/matches/{genome}",
    # matches = output_dir + "/hmmsearch/matches/{genome}.fam",
    matches_fasta = os.path.join(matches_dir, genome)
    fam = os.path.join(matches_dir, genome + ".fam")
    output = os.path.join(output_dir, genome)
    log = os.path.join(logs_dir, f"{genome}_transpose.log")
    # pfam_search, tigrfam_search, hmms_and_names, output, log
    param_sets.append((matches_fasta, fam, output, hmms_and_names, gtdb_bac_tax, gtdb_arc_tax, log))

process_map(process_a_genome, param_sets, max_workers=num_threads, chunksize=1)

logging.info('done')

# touch snakemake.output[0]
with open(snakemake.output[0], 'w') as _: pass
