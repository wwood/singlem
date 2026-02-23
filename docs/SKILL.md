# SingleM Taxonomic Profiling Skill

## Overview

SingleM is a tool for profiling shotgun metagenomes (short- and long-read) by targeting 20 amino acid "window" sequences within single-copy marker genes. It generates GTDB-based taxonomic profiles and is particularly strong at handling novel lineages.

The primary subcommand for taxonomic profiling is `singlem pipe`.

Skill corresponds to SingleM v0.20.3.

---

## Installation

### Conda (recommended)
```bash
conda create -c conda-forge -c bioconda --override-channels \
  --name singlem 'singlem>=0.20.3'
conda activate singlem

# Download reference data (metapackage) — required after conda install
singlem data --output-directory /path/to/metapackage
```

### Docker (includes reference data — no separate data download needed)
```bash
docker pull wwood/singlem:0.20.3
# Run pipe directly:
docker run -v `pwd`:`pwd` wwood/singlem:0.20.3 pipe \
  --sequences `pwd`/my.fastq.gz -p `pwd`/my.profile.csv --threads 4
```

### Singularity/Apptainer
```bash
singularity pull docker://wwood/singlem:0.20.3
singularity run -B `pwd`:`pwd` singlem_0.20.3.sif pipe \
  --sequences `pwd`/my.fastq.gz -p `pwd`/my.profile.csv --threads 4
```

---

## Core Concepts

- **OTU table**: The intermediate output of `singlem pipe`. Contains per-marker-gene OTU sequences with their coverage/abundance across samples.
- **Taxonomic profile** (condensed profile): The final output summarising community composition. Generated from the OTU table via the `condense` algorithm, which uses trimmed means and expectation maximisation across 59 marker genes.
- **Coverage**: The expected per-base coverage of a genome with that OTU sequence. Derived from `num_hits`. The default minimum coverage to report in a taxonomic profile is 0.35× for reads, 0.1× for genomes.
- **GTDB taxonomy**: SingleM uses GTDB taxonomy strings (e.g. `Root; d__Bacteria; p__Proteobacteria; ...`).

---

## Generating a Taxonomic Profile

### Basic usage — paired-end short reads
```bash
singlem pipe \
  --forward sample_R1.fastq.gz \
  --reverse sample_R2.fastq.gz \
  --taxonomic-profile sample.profile.tsv \
  --threads 8
```

### Single-end or unpaired reads
```bash
singlem pipe \
  --sequences sample.fastq.gz \
  -p sample.profile.tsv \
  --threads 8
```

(`-p` is the short form of `--taxonomic-profile`)

### Long reads (Nanopore ≥R10.4.1 or PacBio HiFi)
```bash
singlem pipe \
  --sequences sample_nanopore.fastq.gz \
  -p sample.profile.tsv \
  --threads 8
```
> Long reads use the same interface; SingleM auto-detects read length.

### Multiple samples — combined in one run
```bash
singlem pipe \
  --forward S1_R1.fq.gz S2_R1.fq.gz \
  --reverse S1_R2.fq.gz S2_R2.fq.gz \
  --otu-table all_samples.otu_table.csv \
  --taxonomic-profile all_samples.profile.tsv \
  --threads 16
```
> For >100 samples, run each individually and combine OTU tables with `singlem summarise`.

### Genome / assembly input
```bash
# Single genome
singlem pipe \
  --genome-fasta-files genome.fna \
  -p genome.profile.tsv

# Many genomes from a directory
singlem pipe \
  --genome-fasta-directory /path/to/genomes/ \
  --genome-fasta-extension fna \
  -p genomes.profile.tsv \
  --threads 16

# From a file listing genome paths
singlem pipe \
  --genome-fasta-list genomes.txt \
  -p genomes.profile.tsv \
  --threads 16
```
> Genome mode uses different defaults: higher `--min-taxon-coverage` (0.1) and `--min-orf-length` (300 bp).

---

## Output Options

### Also save an OTU table (`--otu-table`)
Useful for alpha/beta diversity metrics, ordination, and inspecting raw data (e.g. which marker genes fired, which OTU sequences were found). Compatible with `singlem summarise` and `singlem appraise`.
```bash
singlem pipe \
  --forward sample_R1.fastq.gz \
  --reverse sample_R2.fastq.gz \
  --otu-table sample.otu_table.csv \
  --taxonomic-profile sample.profile.tsv \
  --threads 8
```

### Save an archive OTU table (`--archive-otu-table`) — recommended for long-term archiving
The archive OTU table stores additional information (full sequence context, alignment data) needed to regenerate results without re-running the pipeline. It is the right format for two important downstream modes:

- **`singlem condense`** — re-derive the taxonomic profile from the archive OTU table (e.g. with different `--min-taxon-coverage` settings) without re-running `pipe`
- **`singlem renew`** — re-assign taxonomy against an updated metapackage without re-running `pipe`

```bash
singlem pipe \
  --forward sample_R1.fastq.gz \
  --reverse sample_R2.fastq.gz \
  --archive-otu-table sample.archive.otu_table.json.gz \
  --taxonomic-profile sample.profile.tsv \
  --threads 8

# Later: re-derive profile with different coverage threshold
singlem condense \
  --input-archive-otu-tables sample.archive.otu_table.json.gz \
  --taxonomic-profile sample_recondensed.profile.tsv \
  --min-taxon-coverage 0.1

# Later: re-assign taxonomy with a newer metapackage
singlem renew \
  --archive-otu-table sample.archive.otu_table.json.gz \
  --taxonomic-profile sample_updated.profile.tsv \
  --metapackage /path/to/new_metapackage
```

---

## Key Options

| Option | Description | Default |
|--------|-------------|---------|
| `--forward` / `-1` / `--reads` / `--sequences` | Forward or unpaired reads (FASTA/FASTQ, gzipped ok) | required |
| `--reverse` / `-2` | Reverse reads for paired-end | — |
| `--taxonomic-profile` / `-p` | Output taxonomic profile (TSV) | not set |
| `--otu-table` | Output OTU table (CSV) | not set |
| `--threads` | Number of CPU threads | 1 |
| `--metapackage` | Path to reference metapackage | default system metapackage |
| `--min-taxon-coverage` | Min coverage to report in profile | 0.35 (reads), 0.1 (genomes) |
| `--assignment-method` | Taxonomy assignment algorithm for OTUs | `smafa_naive_then_diamond` |
| `--genome-fasta-files` | Input genome FASTA(s) | — |
| `--genome-fasta-directory` / `-d` | Directory of genome FASTAs | — |
| `--genome-fasta-extension` | Extension for genome FASTAs | `fna` |
| `--genome-fasta-list` | File listing genome paths | — |

---

## Output Format

### Taxonomic profile (`-p` / `--taxonomic-profile`)
Tab-separated file with columns:
```
sample    coverage    taxonomy
sample1   5.23        Root; d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Pseudomonadales; f__Pseudomonadaceae; g__Pseudomonas; s__Pseudomonas aeruginosa
sample1   3.10        Root; d__Bacteria; p__Firmicutes_A; ...
```

### OTU table (`--otu-table`)
CSV with columns: `gene`, `sample`, `sequence`, `num_hits`, `coverage`, `taxonomy`

---

## Important Caveats

1. **Use raw reads, not quality-trimmed reads.** Quality trimming (e.g. Trimmomatic) can shorten reads below 100 bp, making them unusable. Adapter trimming is fine but unnecessary.
2. **Do not use assembled contigs as read input.** Use `--genome-fasta-files` for assemblies/MAGs; `--sequences`/`--forward` is for raw reads only.
3. **Reference data required.** After conda install, run `singlem data` before using `pipe`. Docker images include reference data.

---

## Downstream Analysis

### Convert profile to other formats (e.g. BIOM, Kraken-style)
```bash
singlem summarise \
  --input-taxonomic-profiles sample.profile.tsv \
  --output-taxonomic-profile-krona sample.krona.html
```

### Estimate fraction of reads that are bacterial/archaeal (prokaryotic) rather than eukaryotic/phage/etc
```bash
singlem pipe \
  --forward sample_R1.fq.gz --reverse sample_R2.fq.gz \
  -p sample.profile.tsv --threads 8

singlem prokaryotic_fraction \
  --forward sample_R1.fq.gz --reverse sample_R2.fq.gz \
  -p sample.profile.tsv \
  > sample.prokaryotic_fraction.tsv
```

### Re-profile with updated reference database (no re-running pipe)
Requires that the original run saved an `--archive-otu-table`.
```bash
singlem renew \
  --archive-otu-table sample.archive.otu_table.json.gz \
  --taxonomic-profile sample_updated.profile.tsv \
  --metapackage /path/to/new_metapackage
```

### Combine OTU tables from multiple separate runs
```bash
singlem summarise \
  --input-otu-tables s1.otu_table.csv s2.otu_table.csv s3.otu_table.csv \
  --output-otu-table combined.otu_table.csv
```

### Assess how much of a metagenomes's prokaryotes have an associated genome/MAG
```bash
singlem pipe --sequences raw.fq.gz --otu-table metagenome.otu_table.csv
singlem pipe --genome-fasta-files my-genomes/*.fasta --otu-table genomes.otu_table.csv
singlem appraise \
  --metagenome-otu-tables metagenome.otu_table.csv \
  --genome-otu-tables genomes.otu_table.csv
```

---

## Phage Profiling (Lyrebird)

For dsDNA phage profiling, use the `lyrebird` command with the same interface:
```bash
# Download lyrebird reference data
lyrebird data --output-directory /path/to/lyrebird_metapackage

lyrebird pipe \
  --forward sample_R1.fq.gz \
  --reverse sample_R2.fq.gz \
  -p sample.phage_profile.tsv \
  --threads 8
```
Lyrebird uses >500 phage marker genes and vConTACT3-based taxonomy (not GTDB).

---

## Quick Reference — Most Common Commands

```bash
# 1. Download reference data (once, after conda install)
singlem data --output-directory ~/singlem_metapackage

# 2. Profile paired-end metagenome (save archive OTU table for future re-use)
singlem pipe \
  --forward sample_R1.fq.gz \
  --reverse sample_R2.fq.gz \
  --archive-otu-table sample.archive.otu_table.json.gz \
  --taxonomic-profile sample.profile.tsv \
  --threads 16

# 3. View profile
cat sample.profile.tsv

# 4. Convert to Krona chart
singlem summarise \
  --input-taxonomic-profiles sample.profile.tsv \
  --output-taxonomic-profile-krona sample.krona.html
```

---

## Citation

If you use SingleM, please cite:

> Ben J. Woodcroft et al. *Comprehensive taxonomic identification of microbial species in metagenomic data using SingleM and Sandpiper.* Nat Biotechnol (2025). https://doi.org/10.1038/s41587-025-02738-1