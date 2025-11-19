########################################
### transpose_hmms_with_sequences.py ###
########################################
# Author: Rossen Zhao (blatantly cribbed from SingleM script of same name by Samuel Aroney)
# Trawl through extracted sequences matching target HMMs
# Extract these into separate files for each HMM

# input [GENOME ID].faa
# >[SEQUENCE ID]<space>[RANDOM TEXT]
# GENE SEQUENCE...

# input [GENOME ID].txt:
# [SEQUENCE ID]<tab>[HMM ID]

# input taxonomy tsv:
# [GENOME ID]<tab>[taxonomy]

# intended output HMM file [PACKAGE NAME].faa:
# >[GENOME ID]-[PACKAGE NAME]
# GENE SEQUENCE...

# intended output taxonomy file [PACKAGE NAME]_taxonomy.tsv:
# [GENOME ID]-[PACKAGE NAME]<tab>[taxonomy]

import argparse
import logging
import csv
import os
from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="Extract sequences into separate HMM files."
)
parser.add_argument(
    "--input-fasta", type=str, metavar="<INPUT FASTA>", help="path to sequence file"
)
parser.add_argument(
    "--taxonomy", type=str, metavar="<TAX>", nargs='+', help="path/s to taxonomy file"
)
parser.add_argument(
    "--hmm-seq", type=str, metavar="<HMM SEQ LIST>", help="list of HMMs with sequences"
)
parser.add_argument(
    "--hmm-spkg",
    type=str,
    metavar="<HMM SPKG LIST>",
    help="list of HMMs with SingleM package",
)
parser.add_argument(
    "--output", type=str, metavar="<OUTPUT>", help="path to output directory"
)

args = parser.parse_args()
input_path = getattr(args, "input_fasta")
genome_id = os.path.basename(input_path).rsplit(".", 1)[0].replace("_protein", "")
taxfiles = getattr(args, "taxonomy")
HMM_seq_list = getattr(args, "hmm_seq")
HMM_id_list = getattr(args, "hmm_spkg")
NEW_VERSION_NAME_COLUMN = "gene"
NEW_VERSION_HMM_COLUMN = "spkg_name"
output_dir = getattr(args, "output")
os.makedirs(output_dir, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%m/%d/%Y %I:%M:%S %p",
)

logging.info("Creating HMM package-match and sequence-match dictionaries")
with open(HMM_seq_list) as file:
    hmms = csv.reader(file, delimiter="\t")
    Seq_HMM = {line[0]: line[1] for line in hmms}

with open(HMM_id_list) as file:
    hmms = csv.DictReader(file, delimiter="\t")
    HMM_spkg = {}
    hmm_to_spkg = {} # hacky fix for using the wrong column for filenames
    for line in hmms:
        if line[NEW_VERSION_NAME_COLUMN] not in HMM_spkg:
            HMM_spkg[line[NEW_VERSION_NAME_COLUMN]] = line[NEW_VERSION_HMM_COLUMN]
            hmm_to_spkg[line[NEW_VERSION_HMM_COLUMN]] = line[NEW_VERSION_NAME_COLUMN]
        else:
            logging.warning(f"Duplicate HMM name: {line[NEW_VERSION_NAME_COLUMN]}")
            logging.warning(f"Duplicate HMM name: {line[NEW_VERSION_HMM_COLUMN]}")

HMM_output = {spkg: {} for spkg in set(HMM_spkg.values())}

taxonomy = ""
logging.info("Finding matching taxonomy")
for taxfile in taxfiles:
    with open(taxfile) as file:
        tax = csv.reader(file, delimiter="\t")
        for line in tax:
            if line[0] == genome_id:
                taxonomy = line[1]
                break
        else:
            continue

logging.info(f"Load fasta file from {input_path}")
for sequence in SeqIO.parse(input_path, "fasta"):
    HMM_ID_list = Seq_HMM[sequence.id].split(',')
    for HMM_ID in HMM_ID_list:
        if HMM_ID not in HMM_spkg:
            logging.warning(f"Could not find {HMM_ID} in HMM_spkg")
            continue
        spkg = HMM_spkg[HMM_ID]
        HMM_output[spkg] = sequence.seq
        logging.debug(f"Found {HMM_ID} in {spkg}")

example_output_file = os.path.join(output_dir, list(HMM_output.keys())[0] + ".faa")
logging.info(f"Save matching sequences in spkg files e.g. {example_output_file}")
example_taxonomy_output_file = os.path.join(
    output_dir, list(HMM_output.keys())[0] + "_taxonomy.tsv"
)
logging.info(
    f"Save matching taxonomy in spkg files e.g. {example_taxonomy_output_file}"
)
for spkg in HMM_output.keys():
    spkg_name = hmm_to_spkg[spkg]
    if HMM_output[spkg] != {}:
        with open(os.path.join(output_dir, spkg_name + ".faa"), "a") as output_file:
            output_file.write(f">{genome_id}\n")
            output_file.write(f"{HMM_output[spkg]}\n")

        with open(os.path.join(output_dir, spkg_name + "_taxonomy.tsv"), "a") as output_file:
            output_file.write(f"{genome_id}\t{taxonomy}\n")

logging.info("Process complete")
