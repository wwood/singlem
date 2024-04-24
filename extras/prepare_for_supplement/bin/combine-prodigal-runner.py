#!/usr/bin/env python3

import re
import os

# Tab-separated file of genome_fasta<TAB>transcript_fasta<TAB>protein_fasta
# [default: undefined, call genes using Prodigal]

output_dir = snakemake.params['output_directory']
prodigal_runner_output_directory = output_dir + "/prodigal-runner"
output_file = snakemake.output['gene_definitions']
genomes_input = snakemake.params['genomes_input']

r = re.compile(r'^\d+$')

genome_to_paths = {}
def file_to_genome(filename):
    r = re.compile(r'^(.*)\..+')
    if r.match(filename):
        return r.match(filename).group(1)
    else:
        raise Exception("Error: file name {} does not match expected format".format(filename))

# Cache genome to original fasta path
genome_to_fasta = {}
for fasta in genomes_input:
    genome = file_to_genome(os.path.basename(fasta))
    genome_to_fasta[genome] = fasta

# Get paths of transcripts and proteins
for directory in os.listdir(prodigal_runner_output_directory):
    if r.match(directory):
        for file in os.listdir(prodigal_runner_output_directory + "/" + directory):
            genome = file_to_genome(file)
            if genome not in genome_to_paths:
                genome_to_paths[genome] = {}

            if file.endswith('.faa'):
                protein_fasta = prodigal_runner_output_directory + "/" + directory + "/" + file
                if 'protein' in genome_to_paths[genome]:
                    raise Exception("Error: multiple protein fasta files found for genome {}".format(genome))
                genome_to_paths[genome]['protein'] = protein_fasta
            elif file.endswith('.fna'):
                transcripts_fasta = prodigal_runner_output_directory + "/" + directory + "/" + file
                if 'transcripts' in genome_to_paths[genome]:
                    raise Exception("Error: multiple transcript fasta files found for genome {}".format(genome))
                genome_to_paths[genome]['transcripts'] = transcripts_fasta
            elif file.endswith('.gff'):
                pass
            else:
                raise Exception("Error: unexpected file found in directory {}".format(file))

# Write out the gene definitions file
num_genomes = len(genome_to_paths)
with open(output_file, 'w') as out:
    # genome_fasta', 'transcript_fasta', 'protein_fasta
    out.write("\t".join(['genome_fasta', 'transcript_fasta', 'protein_fasta']))
    out.write("\n")

    for genome in genome_to_paths:
        if 'protein' not in genome_to_paths[genome]:
            raise Exception("Error: no protein fasta file found for genome {}".format(genome))
        if 'transcripts' not in genome_to_paths[genome]:
            raise Exception("Error: no transcript fasta file found for genome {}".format(genome))

        out.write("\t".join([
            genome_to_fasta[genome],
            genome_to_paths[genome]['transcripts'],
            genome_to_paths[genome]['protein'],
        ]))
        out.write("\n")
        
print("Wrote {} gene definitions to {}".format(num_genomes, output_file))