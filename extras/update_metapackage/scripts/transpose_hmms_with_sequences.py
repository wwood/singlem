########################################
### transpose_hmms_with_sequences.py ###
########################################
# Author: Samuel Aroney
# Trawl through extracted sequences matching target HMMs
# Extract these into separate files for each HMM

# input [GENOME ID]_protein.fam.faa:
# >[SEQUENCE ID]<space>[RANDOM TEXT]
# GENE SEQUENCE...

# input [GENOME ID]_protein.fam:
# [SEQUENCE ID]<tab>[HMM ID]

# input bacterial/archaeal taxonomy tsv:
# [GENOME ID]<tab>[taxonomy]

# intended output HMM file [PACKAGE NAME].faa:
# >[GENOME ID]-[PACKAGE NAME]
# GENE SEQUENCE...

# intended output taxonomy file [PACKAGE NAME]_taxonomy.tsv:
# [GENOME ID]-[PACKAGE NAME]<tab>[taxonomy]



# Manual testing
# python transpose_hmms_with_sequences.py \
#     --input-fasta ./test/GB_GCA_000091165.1_protein.fam.faa \
#     --bacterial-taxonomy ~/m/db/gtdb/gtdb_release202/bac120_taxonomy_r202.tsv \
#     --archaeal-taxonomy ~/m/db/gtdb/gtdb_release202/ar122_taxonomy_r202.tsv \
#     --hmm-seq ./test/GB_GCA_000091165.1_protein.fam \
#     --hmm-spkg ./hmms_and_names \
#     --output ./test_output/

# Manual testing parallel
# find test/ |
# grep fam$ |
# parallel \
# python transpose_hmms_with_sequences.py \
#     --input-fasta {}.faa \
#     --bacterial-taxonomy ~/m/db/gtdb/gtdb_release202/bac120_taxonomy_r202.tsv \
#     --archaeal-taxonomy ~/m/db/gtdb/gtdb_release202/ar122_taxonomy_r202.tsv \
#     --hmm-seq {} \
#     --hmm-spkg ./hmms_and_names \
#     --output ./test_output/



import argparse
import logging
import csv
import os
import re
from Bio import SeqIO


parser = argparse.ArgumentParser(description='Extract sequences into separate HMM files.')
parser.add_argument('--input-fasta', type=str, metavar='<INPUT FASTA>', help='path to sequence file')
parser.add_argument('--bacterial-taxonomy', type=str, metavar='<BAC TAX>', help='path to bacterial taxonomy file')
parser.add_argument('--archaeal-taxonomy', type=str, metavar='<ARCH TAX>', help='path to archaeal taxonomy file')
parser.add_argument('--hmm-seq', type=str, metavar='<HMM SEQ LIST>', help='list of HMMs with sequences')
parser.add_argument('--hmm-spkg', type=str, metavar='<HMM SPKG LIST>', help='list of HMMs with SingleM package')
parser.add_argument('--output', type=str, metavar='<OUTPUT>', help='path to output directory')

args = parser.parse_args()
input_path = getattr(args, 'input_fasta')
genome_id = re.findall(r'(.*)\.faa', os.path.basename(input_path))[0].removesuffix("_protein")
bac_taxonomy = getattr(args, 'bacterial_taxonomy')
arch_taxonomy = getattr(args, 'archaeal_taxonomy')
HMM_seq_list = getattr(args, 'hmm_seq')
HMM_id_list = getattr(args, 'hmm_spkg')
NEW_VERSION_HMM_COLUMN = "HMM"
NEW_VERSION_NAME_COLUMN = "name"
output_dir = getattr(args, 'output')
os.makedirs(output_dir, exist_ok=True)

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')


logging.info("Creating HMM package-match and sequence-match dictionaries")
with open(HMM_seq_list) as file:
    hmms = csv.reader(file, delimiter="\t")
    Seq_HMM = {line[0]:line[1] for line in hmms}

with open(HMM_id_list) as file:
    hmms = csv.DictReader(file, delimiter="\t")
    HMM_spkg = {line[NEW_VERSION_HMM_COLUMN]:line[NEW_VERSION_NAME_COLUMN] for line in hmms}

HMM_output = {spkg:{} for spkg in set(HMM_spkg.values())}


logging.info("Finding matching bacterial or archaeal taxonomy")
with open(bac_taxonomy) as file:
    bac_tax = csv.reader(file, delimiter="\t")
    taxonomy = [line[1] for line in bac_tax if line[0] == genome_id]
    
if len(taxonomy) == 0:
    with open(arch_taxonomy) as file:
        arch_tax = csv.reader(file, delimiter="\t")
        taxonomy = [line[1] for line in arch_tax if line[0] == genome_id]


logging.info(f"Load fasta file from {input_path}")
for sequence in SeqIO.parse(input_path, "fasta"):
    HMM_ID = Seq_HMM[sequence.id]
    spkg = HMM_spkg[HMM_ID]
    HMM_output[spkg] = sequence.seq
    

example_output_file = os.path.join(output_dir, list(HMM_output.keys())[0] + '.faa')
logging.info(f"Save matching sequences in spkg files e.g. {example_output_file}")
example_taxonomy_output_file = os.path.join(output_dir, list(HMM_output.keys())[0] + '_taxonomy.tsv')
logging.info(f"Save matching taxonomy in spkg files e.g. {example_taxonomy_output_file}")
for spkg in HMM_output.keys():
    if HMM_output[spkg] != {}:
        with open(os.path.join(output_dir, spkg + ".faa"), 'a') as output_file:
            # print(f">{genome_id}")
            # print(HMM_output[spkg])
            output_file.write(f">{genome_id}\n")
            output_file.write(f"{HMM_output[spkg]}\n")

        with open(os.path.join(output_dir, spkg + "_taxonomy.tsv"), 'a') as output_file:
            # print(f"{genome_id}\t{taxonomy[0]}")
            output_file.write(f"{genome_id}\t{taxonomy[0]}\n")

logging.info("Process complete")




