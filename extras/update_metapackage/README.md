See top of Snakefile for some instructions on how to run this pipeline.

## Uniprot processing

Download latest Uniprot swissprot annotations

```bash
# Main swissprot annotations - into sprot_to_fasta_and_taxonomy
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz

# Running sprot_to_fasta_and_taxonomy.py:
pixi run -e update-metapackage \
    python sprot_to_fasta_and_taxonomy.py -s uniprot_sprot.dat.gz -f uniprot_sprot.fa -t uniprot_sprot_taxonomy.tsv

# Output
# Viruses   17497
# Eukaryota 200289
# Bacteria  336999
# Archaea   19842
```

## Rename preliminary taxonomy files

```python
# When taxonomy files have genome names like: G001261915 instead of RS_GCF_015356595.1
import polars as pl
import os

genome_files = (
    pl.DataFrame({"filename": [f for r, d, files in os.walk("protein_faa_reps") for f in files if f.endswith(".faa.gz")]})
    .select(
        genome_id = pl.col("filename").str.replace(r"_protein\.faa\.gz$", ""),
        g_id = pl.col("filename").str.replace(r"\.\d+_protein\.faa\.gz$", "").str.replace(r"^[RG][SB]_GC[FA]_", "G"),
    )
)

fixed_ar = (
    pl.read_csv("gtdb_r232_ar53_draft_taxonomy.tsv", separator="\t", has_header=False, new_columns=["g_id", "taxonomy"])
    .join(genome_files, on="g_id", how="left")
    .select("genome_id", "taxonomy")
)

fixed_ar.filter(pl.col("genome_id").is_null())
fixed_ar.write_csv("ar53_taxonomy_r232.tsv", separator="\t", include_header=False)

fixed_bac = (
    pl.read_csv("gtdb_r232_bac120_draft_taxonomy.tsv", separator="\t", has_header=False, new_columns=["g_id", "taxonomy"])
    .join(genome_files, on="g_id", how="left")
    .select("genome_id", "taxonomy")
)

fixed_bac.filter(pl.col("genome_id").is_null())
fixed_bac.write_csv("bac120_taxonomy_r232.tsv", separator="\t", include_header=False)
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

## Compare outputs

```r
library(tidyverse)

old <- read_csv("/work/microbiome/msingle/mess/174_R220_renew/processing_20240531/per_acc_summary.csv") %>%
    select(sample, old_species_coverage = species_coverage, old_read_fraction = read_fraction)
new <- read_csv("/work/microbiome/msingle/mess/198_R226_renew_mach2/processing_20250409/per_acc_summary.csv") %>%
    select(sample, new_species_coverage = species_coverage, new_read_fraction = read_fraction)

comb <- old %>%
    inner_join(new)

comb %>%
    ggplot(aes(x = old_species_coverage, y = new_species_coverage)) +
    geom_point(alpha = 0.1) +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(x = "Old known species fraction", y = "New known species fraction") +
    theme_bw() +
    theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5)
    )
ggsave("old_vs_new_species_fraction.png", width = 8, height = 6, dpi = 300)

comb %>%
    ggplot(aes(x = old_read_fraction / 100, y = new_read_fraction / 100)) +
    geom_point(alpha = 0.1) +
    geom_abline(slope = 1, intercept = 0, colour = "red") +
    xlim(0, 1) +
    ylim(0, 1) +
    labs(x = "Old microbial fraction", y = "New microbial fraction") +
    theme_bw() +
    theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5)
    )
ggsave("old_vs_new_prokaryotic_fraction.png", width = 8, height = 6, dpi = 300)
```
