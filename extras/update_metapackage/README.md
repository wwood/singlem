See top of Snakefile for some instructions on how to run this pipeline.

## Uniprot processing

Download latest Uniprot swissprot annotations

```bash
# Main swissprot annotations - into sprot_to_fasta_and_taxonomy
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz

# Running sprot_to_fasta_and_taxonomy.py:
python sprot_to_fasta_and_taxonomy.py -s uniprot_sprot.dat.gz -f uniprot_sprot.fa -t uniprot_sprot_taxonomy.tsv

# Output
# Viruses    17451
# Eukaryota  198986
# Bacteria   336739
# Archaea    19794
```

## SRA test run comparison

```r
library(tidyverse)
samples <- tribble(
    ~sample, ~metapackage,
    "SRR1926147", "old",
    "SRR1926147", "new",
    "DRR083195", "old",
    "DRR083195", "new",
    "SRR9224014", "old",
    "SRR9224014", "new",
)

# Check for increases in lower-level relative abundance
tax_fraction <- samples %>%
    mutate(
        tax_fraction_path = str_c("sra/", sample, "_", metapackage, ".taxonomic_coverage.tsv"),
        tax_fraction = map(tax_fraction_path, ~ read_tsv(.x)),
    ) %>%
    unnest(tax_fraction, names_sep="-") %>%
    select(-tax_fraction_path) %>%
    select(sample, metapackage, level = `tax_fraction-level`, rel_abund = `tax_fraction-relative abundance (%)`) %>%
    pivot_wider(names_from = metapackage, values_from = rel_abund)

# Check for similar SMF and average genome size values
smf <- samples %>%
    mutate(
        smf_path = str_c("sra/", sample, "_", metapackage, ".smf.tsv"),
        smf = map(smf_path, ~ read_tsv(.x)),
    ) %>%
    unnest(smf, names_sep="-") %>%
    select(-smf_path) %>%
    select(sample, metapackage, read_fraction = `smf-read_fraction`, average_size = `smf-average_bacterial_archaeal_genome_size`) %>%
    pivot_wider(names_from = metapackage, values_from = c(read_fraction, average_size))
```

## Run renew

```bash
# in mess/198_R226_renew_mach2
notify /work/microbiome/msingle/mess/195_sandpiper1/novel_genome_finder/bin/find_novel_genomes_by_renew.py renew \
    --sample-to-archive-tsv /work/microbiome/msingle/mess/198_R226_renew_mach2/final_unannotated_list.tsv \
    --all-samples  --output-directory renew \
    --new-metapackage /work/microbiome/msingle/mess/196_metapackage_r226/update_metapackage_gtdb_transcripts/metapackage/S5.4.0.GTDB_r226.metapackage_20250331.smpkg \
    --run-through-mqsub
```
