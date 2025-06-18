import os
import logging
import pandas as pd

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')

spkgs = snakemake.input.spkgs
match_directory = snakemake.input.match_directory
resolved_trees_list = snakemake.input.resolved_trees_list

resolved_gene_list = set()
logging.info(os.path.basename(__file__) + ": Processing resolved trees")
with open(resolved_trees_list) as f:
    for line in f.readlines()[1:]:
        gene = line.strip().split("\t")[0]
        resolved_gene_list.add(gene)

gene_to_spkg = {}
for spkg in spkgs:
    gene = os.path.basename(spkg).rsplit(".", 1)[0]
    if gene not in resolved_gene_list:
        continue
    gene_to_spkg[gene] = spkg

logging.info(os.path.basename(__file__) + ": Processing {} genes".format(len(gene_to_spkg)))

gene_to_species = {}
species_to_gene = {}
species_coverage = {}
picked_genes = set()

for filename in os.listdir(match_directory):
    species = filename.rsplit(".", 1)[0]
    with open(os.path.join(match_directory, filename)) as f:
        gene_in_species = False
        for line in f:
            __, gene = line.strip().split("\t")
            if gene not in gene_to_spkg:
                continue
            if gene not in gene_to_species:
                gene_to_species[gene] = set()
            gene_to_species[gene].add(species)
            if species not in species_to_gene:
                species_to_gene[species] = set()
            species_to_gene[species].add(gene)
            gene_in_species = True
        if gene_in_species:
            species_coverage[species] = 0

num_iterations = 4 #3

logging.info(os.path.basename(__file__) + f": Found {len(species_coverage)} of {len(os.listdir(match_directory))} total species covered by at least one gene")
logging.info(os.path.basename(__file__) + f": Starting greedy search to select {num_iterations} genes per species if possible.")

iterations = 1
while iterations <= num_iterations:
    for species in species_coverage:
        if species_coverage[species] > iterations:
            continue
        sorted_genes = sorted([(len(gene_to_species[gene]), gene) for gene in species_to_gene[species]], reverse=True)
        for __, gene in sorted_genes:
            if gene not in picked_genes:
                picked_genes.add(gene)
                for s in gene_to_species[gene]:
                    species_coverage[s] += 1
                break
    logging.info(os.path.basename(__file__) + ": Iteration {}: Picked {} genes".format(iterations, len(picked_genes)))
    iterations += 1

logging.info(os.path.basename(__file__) + ": Picked a total of {} genes from greedy search".format(len(picked_genes)))

old_hmms_and_names = pd.read_csv(snakemake.params.hmms_and_names, sep='\t', header=None, names=['gene', 'spkg_name', 'filepath', 'description'])
hmms_and_names = old_hmms_and_names[old_hmms_and_names['gene'].isin(picked_genes)]
logging.info(os.path.basename(__file__) + ": Writing to " + snakemake.output.hmms_and_names_roundrobin)
hmms_and_names.to_csv(snakemake.output.hmms_and_names_roundrobin, sep='\t', header=True, index=False)

outfile = snakemake.output.coverages
logging.info(os.path.basename(__file__) + ": Writing species coverages to " + outfile)
with open(outfile, "w") as f:
    for species in sorted(species_coverage.keys()):
        f.write("{}\t{}\n".format(species, species_coverage[species]))

logging.info('done')
with open(snakemake.output.done, 'w') as _: pass