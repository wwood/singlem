import os
import logging
import pathlib
import extern

# dups are probably if there is more fasta entries than there are taxon files
# append -1, -2, -3 to the fasta entry id if there are duplicates, then add the duplicated ids to the taxon file
# taxon file is in tsv format of [GENOME ID]\t[taxonomy]

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

fasta = snakemake.input.seqs
taxon_file = snakemake.input.taxonomy

output_fasta = snakemake.output.seqs
output_taxon_file = snakemake.output.taxonomy

logging.info(f"Renaming off-target duplicates in {fasta} and {taxon_file}")
# Read in the fasta file and create a mapping of sequence IDs to their counts
num_seqs = 0
seq_id_counts = {}
with open(fasta, 'r') as fasta_file:
    for line in fasta_file:
        if line.startswith('>'):
            num_seqs += 1
            seq_id = line[1:].strip().split()[0]  # Get the sequence ID
            seq_id_counts[seq_id] = seq_id_counts.get(seq_id, 0) + 1
logging.info(f"Found {len(seq_id_counts)} unique sequence IDs in the fasta file.")
logging.info(f"Total number of sequences in fasta: {num_seqs}")
# Create a mapping of sequence IDs to their new names
seq_id_mapping = {}
for seq_id, count in seq_id_counts.items():
    if count > 1:
        logging.info(f"Found {count} duplicates for sequence ID {seq_id}, renaming them.")
        for i in range(1, count + 1):
            new_seq_id = f"{seq_id}-{i}"
            seq_id_mapping[seq_id] = new_seq_id
    else:
        seq_id_mapping[seq_id] = seq_id
logging.info(f"Renamed {len(seq_id_mapping)} sequence IDs in the fasta file.")
# Write the renamed fasta file
with open(output_fasta, 'w') as fasta_out:
    with open(fasta, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith('>'):
                seq_id = line[1:].strip().split()[0]  # Get the sequence ID
                new_seq_id = seq_id_mapping[seq_id]
                fasta_out.write(f">{new_seq_id}\n")
            else:
                fasta_out.write(line)
logging.info(f"Wrote renamed fasta file to {output_fasta}")
# Write the updated taxonomy file
with open(output_taxon_file, 'w') as taxon_out:
    with open(taxon_file, 'r') as taxon_file:
        for line in taxon_file:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                seq_id = parts[0]
                taxonomy = parts[1]
                new_seq_id = seq_id_mapping.get(seq_id, seq_id)
                taxon_out.write(f"{new_seq_id}\t{taxonomy}\n")
logging.info(f"Wrote updated taxonomy file to {output_taxon_file}")
logging.info("Renaming off-target duplicates completed successfully.")

with open(snakemake.output.done, 'w') as _: pass