################################################################################
## Script for getting mean genome sizes per rank from the GTDB
## Raphael Eisenhofer June 2022
################################################################################

## Load libraries
library(tidyverse)

## Import arguments
## N.B. arg[1] = gtdb_bac_metadata, arg[2] = gtdb_arc_metadata, arg[3] outputdir
args = commandArgs(trailingOnly=TRUE)

## Import data
gtdb_bac_metadata <- read_tsv(args[1])
gtdb_arc_metadata <- read_tsv(args[2])

## Bind bac/arc metadata, #OPTIONAL: filter by minimnum/max completeness/contamination
## Also split out taxonomy strink into columns by taxonomic rank
gtdb_combined_metadata <- rbind(gtdb_arc_metadata, gtdb_bac_metadata) %>%
#  filter(checkm_completeness > 80 & checkm_contamination < 10) %>%
#  mutate(gtdb_taxonomy = str_replace_all(gtdb_taxonomy, ".__", "")) %>%
  separate(., col = gtdb_taxonomy,
           sep = ";",
           into = c("domain", "phylum", "class", "order", "family", "genus", "species")
           )

################################################################################
## Mean genome sizes per taxonomic rank
# Note for family we take the mean of the genus means, not the mean of the
# species that belong to a family. This helps normalise for species that are
# overabundant in the database and that would otherwise inflate estimates -- e.g. E.coli!

genome_size_by_rank <- gtdb_combined_metadata %>%
  group_by(domain, phylum, class, order, family, genus, species) %>%
  summarise(species_mean = mean(genome_size)) %>%
  ungroup()

species_means <- genome_size_by_rank %>%
  group_by(domain, phylum, class, order, family, genus, species) %>%
  summarise(genome_size = mean(species_mean)) %>%
  ungroup() %>%
  select(species, genome_size) %>%
  rename("rank" = species)

genus_means <- genome_size_by_rank %>%
  group_by(domain, phylum, class, order, family, genus) %>%
  summarise(genome_size = mean(species_mean)) %>%
  ungroup() %>%
  select(genus, genome_size) %>%
  rename("rank" = genus)

genome_size_genus_means <- genome_size_by_rank %>%
  group_by(domain, phylum, class, order, family, genus) %>%
  summarise(genome_size = mean(species_mean)) %>%
  ungroup() 

family_means <- genome_size_genus_means %>%
  group_by(domain, phylum, class, order, family) %>%
  summarise(genome_size = mean(genome_size)) %>%
  ungroup() %>%
  select(family, genome_size) %>%
  rename("rank" = family)

order_means <- genome_size_genus_means %>%
  group_by(domain, phylum, class, order) %>%
  summarise(genome_size = mean(genome_size)) %>%
  ungroup() %>%
  select(order, genome_size) %>%
  rename("rank" = order)

class_means <- genome_size_genus_means %>%
  group_by(domain, phylum, class) %>%
  summarise(genome_size = mean(genome_size)) %>%
  ungroup() %>%
  select(class, genome_size) %>%
  rename("rank" = class)

phylum_means <- genome_size_genus_means %>%
  group_by(domain, phylum) %>%
  summarise(genome_size = mean(genome_size)) %>%
  ungroup() %>%
  select(phylum, genome_size) %>%
  rename("rank" = phylum)

domain_means <- genome_size_genus_means %>%
  group_by(domain) %>%
  summarise(genome_size = mean(genome_size)) %>%
  ungroup() %>%
  select(domain, genome_size) %>%
  rename("rank" = domain)

#Combine
full_genome_size_table <- rbind(species_means, genus_means, family_means,
                                order_means, class_means, phylum_means,
                                domain_means)

write_tsv(full_genome_size_table, paste0(args[3], "/", "gtdb_mean_genome_sizes.tsv"))
