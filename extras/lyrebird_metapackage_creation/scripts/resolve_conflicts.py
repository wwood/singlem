### Identifies individual proteins that match to multiple HMMs, and chooses the HMM that covers the greatest number of species

import os
import logging
import pandas as pd

import pathlib

match_directory = snakemake.params.match_directory
output_dir = snakemake.params.output_dir

pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

logging.basicConfig(filename=snakemake.log[0],level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

proteins_to_hmm = {}
hmm_to_proteins = {}
hmm_to_species = {}
species_list = []

resolved_genes = set()
logging.info(f"Reading hmmsearch matches from {match_directory}")
for filename in os.listdir(match_directory):
    if filename.endswith(".tsv"):
        with open(os.path.join(match_directory, filename)) as f:
            species_list.append(filename.rsplit(".", 1)[0])
            for line in f:
                protein, hmms = line.strip().split('\t')
                for hmm in hmms.split(','):
                    if protein not in proteins_to_hmm:
                        proteins_to_hmm[protein] = [hmm]
                    else:
                        proteins_to_hmm[protein].append(hmm)
                    if hmm not in hmm_to_proteins:
                        hmm_to_proteins[hmm] = [protein]
                    else:
                        hmm_to_proteins[hmm].append(protein)
                    if hmm not in hmm_to_species:
                        hmm_to_species[hmm] = set()
                        hmm_to_species[hmm].add(protein.rsplit('_', 1)[0])
                    else:
                        hmm_to_species[hmm].add(protein.rsplit('_', 1)[0])

logging.info(f"Found {len(proteins_to_hmm)} proteins and {len(hmm_to_proteins)} HMMs in {len(species_list)} species")

# pick for copy number < 1.05
for hmm in hmm_to_species:
    if len(hmm_to_proteins[hmm])/len(hmm_to_species[hmm]) > 1.05:
        hmm_to_proteins.pop(hmm, None)
    elif len(hmm_to_proteins[hmm]) < 5:
        hmm_to_proteins.pop(hmm, None)
logging.info(f"Kept {len(hmm_to_proteins)} HMMs with copy number <= 1.05 and >= 5 proteins")
del(hmm_to_species)

picked_hmms = set()
for hmm in hmm_to_proteins:
    picked_hmms.add(hmm)

# reverse roundrobin
iterations = 1
while True:
    break_outer = True
    for protein in proteins_to_hmm:
        check_hmms = [hmm for hmm in proteins_to_hmm[protein] if hmm in picked_hmms]
        if len(check_hmms) <= 1:
            continue
        break_outer = False
        hmms_sorted = sorted(check_hmms, key=lambda x: len(hmm_to_proteins[x]), reverse=True)
        for hmm in hmms_sorted[1:]:
            picked_hmms.remove(hmm)
    if break_outer:
        break
    iterations += 1

logging.info(f"Found {len(picked_hmms)} HMMs after {iterations} iterations")

species_to_prot_pair = {}

for hmm in picked_hmms:
    for protein in hmm_to_proteins[hmm]:
        species = protein.rsplit("_", 1)[0]
        if species not in species_to_prot_pair:
            species_to_prot_pair[species] = [(protein, hmm)]
        else:
            species_to_prot_pair[species].append((protein, hmm))

logging.info(f"Found {len(species_to_prot_pair)} species with at least one protein, with {len(species_list)-len(species_to_prot_pair)} unable to be resolved.")

for species in species_to_prot_pair:
    with open(os.path.join(output_dir, f"{species}.tsv"), "w") as f:
        for protein, hmm in species_to_prot_pair[species]:
            f.write(protein + "\t" + hmm + "\n")

logging.info(f"Finished writing {len(species_to_prot_pair)} species files to {output_dir}")

old_hmms_and_names = pd.read_csv(snakemake.params.hmms_and_names, sep='\t', header=None, names=['gene', 'spkg_name', 'filepath', 'description'])
hmms_and_names = old_hmms_and_names[old_hmms_and_names['gene'].isin(picked_hmms)]
hmms_and_names.to_csv(snakemake.output.hmms_and_names_noconflict, sep='\t', header=True, index=False)

with open(snakemake.output[0], 'w') as _: pass