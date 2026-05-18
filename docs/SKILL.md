---
name: singlem
description: Use SingleM to profile metagenomes and genomes, generate OTU tables, and produce GTDB-based taxonomic profiles.
---

# SingleM Taxonomic Profiling Skill

## Overview

SingleM is a tool for profiling shotgun metagenomes (short- and long-read) by targeting 20 amino acid "window" sequences within single-copy marker genes. It generates GTDB-based taxonomic profiles and is particularly strong at handling novel lineages.

The primary subcommand for taxonomic profiling is `singlem pipe`.

Skill corresponds to SingleM v0.21.2.

---

## Installation

### Conda (recommended)
```bash
conda create -c conda-forge -c bioconda --override-channels \
  --name singlem 'singlem>=0.21.2'
conda activate singlem

# Download reference data (metapackage) — required after conda install
singlem data --output-directory /path/to/metapackage
```

### Docker (includes reference data — no separate data download needed)
```bash
docker pull wwood/singlem:0.21.2
# Run pipe directly:
docker run -v `pwd`:`pwd` wwood/singlem:0.21.2 pipe \
  --sequences `pwd`/my.fastq.gz -p `pwd`/my.profile.csv --threads 4
```

### Singularity/Apptainer
```bash
singularity pull docker://wwood/singlem:0.21.2
singularity run -B `pwd`:`pwd` singlem_0.21.2.sif pipe \
  --sequences `pwd`/my.fastq.gz -p `pwd`/my.profile.csv --threads 4
```

---

## Core Concepts

- **OTU table**: The intermediate output of `singlem pipe`. Contains per-marker-gene OTU sequences with their coverage/abundance across samples.
- **Taxonomic profile** (condensed profile): The final output summarising community composition. Generated from the OTU table via the `condense` algorithm, which uses trimmed means and expectation maximisation across 59 marker genes.
- **Coverage**: The expected per-base coverage of a genome with that OTU sequence. Derived from `num_hits`. The default minimum coverage to report in a taxonomic profile is 0.35× for reads, 0.1× for genomes.
- **GTDB taxonomy**: SingleM uses GTDB taxonomy strings (e.g. `Root; d__Bacteria; p__Proteobacteria; ...`).

---

## SingleM Subcommands at a Glance

SingleM (and its phage-focused sibling Lyrebird) is a suite of subcommands. Most users only need `pipe` (and `data` once, to fetch reference data). The rest support downstream analysis, reference-data management, and package development.

### Main tools

| Subcommand | Purpose |
|------------|---------|
| `singlem pipe` | Main workflow: profile reads/genomes → OTU table + GTDB taxonomic profile |
| `singlem data` | Download / verify the reference metapackage |
| `singlem summarise` | Mechanical transformations of `pipe` results (Krona, species-by-site tables, combining OTU tables, etc.) |
| `singlem renew` | Re-run taxonomy assignment on an existing archive OTU table against a new metapackage |
| `singlem supplement` | Add new genomes to a metapackage to create a custom reference |
| `singlem prokaryotic_fraction` | Estimate the bacterial/archaeal fraction (and average genome size) of a metagenome |
| `singlem appraise` | Assess how much of a metagenome is represented by a set of genomes/assemblies |
| `lyrebird data` | Download / verify the Lyrebird (phage) reference metapackage |
| `lyrebird pipe` | Profile dsDNA phages — same interface as `singlem pipe` |

### Advanced / expert modes

| Subcommand | Purpose |
|------------|---------|
| `singlem condense` | Generate a taxonomic profile from an existing (archive) OTU table |
| `singlem makedb` | Build a searchable database (`.sdb`) from OTU tables |
| `singlem query` | Find sequences in a `makedb` database similar to query OTU sequences |
| `singlem seqs` | Choose the best window position within an HMM (step 1 of building a SingleM package) |
| `singlem create` | Create a SingleM package from a GraftM package + taxonomy (step 2 of package building) |
| `singlem regenerate` | Update an existing SingleM package with new sequences/taxonomy |
| `singlem metapackage` | Create (or `--describe`) a metapackage from individual SingleM packages |
| `lyrebird condense` | `condense` for Lyrebird (non-universal phage markers) |
| `lyrebird renew` | `renew` for Lyrebird archive OTU tables |

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
  --input-archive-otu-table sample.archive.otu_table.json.gz \
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

### Taxonomic profile (`-p` / `--taxonomic-profile`) — SingleM condensed format
Tab-separated file (`.tsv`) with three columns: `sample`, `coverage`, `taxonomy`.
```
sample     coverage  taxonomy
marine0.1  3.64      Root; d__Archaea
marine0.1  0.02      Root; d__Bacteria
marine0.1  0.56      Root; d__Archaea; p__Thermoproteota
marine0.1  0.80      Root; d__Bacteria; p__Desulfobacterota
marine0.1  2.17      Root; d__Bacteria; p__Proteobacteria
```

Key properties of the condensed format:
- **Coverage** is the estimated per-base read coverage attributed *directly* to that taxon — it is **not** inclusive of descendants. For example, `Root; d__Bacteria` (coverage 0.02) does not include the coverage from `p__Desulfobacterota` (0.80) or `p__Proteobacteria` (2.17); those are reported on their own lines.
- Every taxonomic level that has any coverage is listed as its own row. Higher-level rows (e.g. domain) represent reads that could not be assigned more specifically.
- Taxonomy strings follow GTDB conventions: `Root; d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; ...`
- Multiple samples can appear in a single file, distinguished by the `sample` column.

### OTU table (`--otu-table`)
CSV with columns: `gene`, `sample`, `sequence`, `num_hits`, `coverage`, `taxonomy`

---

## Important Caveats

1. **Use raw reads, not quality-trimmed reads.** Quality trimming (e.g. Trimmomatic) can shorten reads below 100 bp, making them unusable. Adapter trimming is fine but unnecessary.
2. **Do not use assembled contigs as read input.** Use `--genome-fasta-files` for assemblies/MAGs; `--sequences`/`--forward` is for raw reads only.
3. **Reference data required.** After conda install, run `singlem data` before using `pipe`. Docker images include reference data.

---

## Downstream Analysis

### Convert condensed profile to other formats (`singlem summarise`)

`singlem summarise` transforms the condensed profile into several more analysis-friendly formats.

#### Krona chart (interactive HTML)
```bash
singlem summarise \
  --input-taxonomic-profile sample.profile.tsv \
  --output-taxonomic-profile-krona sample.krona.html
```
Produces an interactive hierarchical chart viewable in any web browser. Can also be generated directly from `pipe` with `--taxonomic-profile-krona`.

#### Relative abundance species-by-site table
Outputs a taxon-by-sample matrix with relative abundance as percentages. Use `--output-species-by-site-level` to choose the taxonomic rank (`domain`, `phylum`, `class`, `order`, `family`, `genus`, or `species`):
```bash
singlem summarise \
  --input-taxonomic-profile sample.profile.tsv \
  --output-species-by-site-relative-abundance sample.phylum.csv \
  --output-species-by-site-level phylum
```
Example output (one column per sample when multiple samples are present):
```
taxonomy                                marine0.1
unassigned                              50.9
Root; d__Archaea; p__Thermoproteota     7.79
Root; d__Bacteria; p__Desulfobacterota  11.13
Root; d__Bacteria; p__Proteobacteria    30.18
```
To generate tables for all taxonomic levels at once, use a prefix:
```bash
singlem summarise \
  --input-taxonomic-profile sample.profile.tsv \
  --output-species-by-site-relative-abundance-prefix myprefix
# produces: myprefix-domain.tsv, myprefix-phylum.tsv, ..., myprefix-species.tsv
```

#### Long form with extra columns (filled coverage, relative abundance, level)
```bash
singlem summarise \
  --input-taxonomic-profile sample.profile.tsv \
  --output-taxonomic-profile-with-extras sample.with_extras.tsv
```
Adds `full_coverage` (coverage including descendants), `relative_abundance` (%), and `level` columns:
```
sample     coverage  full_coverage  relative_abundance  level   taxonomy
marine0.1  0         7.19           100.0               root    Root
marine0.1  3.64      4.20           58.41               domain  Root; d__Archaea
marine0.1  0.02      2.99           41.59               domain  Root; d__Bacteria
marine0.1  0.56      0.56           7.79                phylum  Root; d__Archaea; p__Thermoproteota
marine0.1  0.80      0.80           11.13               phylum  Root; d__Bacteria; p__Desulfobacterota
marine0.1  2.17      2.17           30.18               phylum  Root; d__Bacteria; p__Proteobacteria
```
Note: `coverage` here is *unfilled* (not including descendants); `full_coverage` is *filled* (sum of a taxon and all its descendants).

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
  --input-archive-otu-table sample.archive.otu_table.json.gz \
  --taxonomic-profile sample_updated.profile.tsv \
  --metapackage /path/to/new_metapackage
```
`renew` also accepts `--assignment-method`, `--threads`, and `--min-taxon-coverage`, just like `pipe`.

### Combine OTU tables from multiple separate runs
```bash
singlem summarise \
  --input-otu-tables s1.otu_table.csv s2.otu_table.csv s3.otu_table.csv \
  --output-otu-table combined.otu_table.csv
```

### Assess how much of a metagenome's prokaryotes have an associated genome/MAG (`singlem appraise`)
`appraise` compares OTU sequences from genomes and/or assemblies against those from the raw metagenome, reporting which lineages are represented and which are missing.
```bash
singlem pipe --sequences raw.fq.gz --otu-table metagenome.otu_table.csv
singlem pipe --genome-fasta-files my-genomes/*.fasta --otu-table genomes.otu_table.csv
singlem appraise \
  --metagenome-otu-tables metagenome.otu_table.csv \
  --genome-otu-tables genomes.otu_table.csv
```
Useful extras:
- `--assembly-otu-tables` — appraise an assembly alongside (or instead of) binned genomes.
- `--imperfect` — match OTU sequences that are similar but not identical (e.g. to credit a genus-level representative); tune with `--sequence-identity`.
- `--plot appraise.svg` — render the appraisal visually (one rectangle per OTU sequence, sized by abundance).
- `--output-binned-otu-table` / `--output-unbinned-otu-table` / `--output-unaccounted-for-otu-table` — write OTU tables of the represented vs. missing populations.

---

## Advanced & Expert Modes

These subcommands support custom reference data and lower-level analyses. Most users never need them.

### Add genomes to a reference metapackage (`singlem supplement`)
Creates a new metapackage that includes your genomes, so future `pipe` runs can identify them. Taxonomy for the new genomes is assigned with GTDB-Tk (installed separately, with a version matching the metapackage's GTDB release) unless supplied via `--taxonomy-file` or `--new-fully-defined-taxonomies`.
```bash
singlem supplement \
  --new-genome-fasta-files genome1.fna genome2.fna \
  --input-metapackage /path/to/metapackage \
  --output-metapackage supplemented.smpkg \
  --checkm2-quality-file checkm2_quality.tsv \
  --dereplicate-with-galah \
  --threads 8
```
A dereplication mode is required: either `--dereplicate-with-galah` (run galah at species level) or `--no-dereplication` (inputs are already dereplicated). A quality-filtering choice is also required: pass CheckM2 results with `--checkm2-quality-file`, or skip with `--no-quality-filter` (and optionally `--no-taxon-genome-lengths` if no CheckM2 file is supplied).

### Build and query a SingleM database (`singlem makedb` / `singlem query`)
Useful for asking "is this OTU sequence (or anything similar) present in samples B, C, D?". `.sdb` is the conventional database extension.
```bash
# Build a database from OTU tables
singlem makedb \
  --otu-tables B.otu_table.csv C.otu_table.csv D.otu_table.csv \
  --db BCD.sdb

# Find database sequences within a given divergence of query OTUs
singlem query \
  --db BCD.sdb \
  --query-otu-table A.otu_table.csv \
  --max-divergence 3
```
`query` can also dump database contents filtered by sample (`--sample-names`), by taxonomy (`--taxonomy Archaea`), or in full (`--dump`).

### Re-derive a profile from an OTU table (`singlem condense`)
`condense` turns an archive OTU table into a taxonomic profile. It is normally invoked implicitly by `pipe`'s `-p` / `--taxonomic-profile`, but can be run standalone — e.g. to recompute a profile with a different `--min-taxon-coverage` without re-running `pipe`. See "Save an archive OTU table" under Output Options for an example.

### Create or inspect a metapackage (`singlem metapackage`)
Assemble individual SingleM packages (`.spkg`) into a metapackage, or inspect an existing one with `--describe`.
```bash
# Describe the contents of an existing metapackage
singlem metapackage --metapackage /path/to/metapackage --describe

# Create a metapackage from individual packages
singlem metapackage \
  --singlem-packages pkg1.spkg pkg2.spkg \
  --metapackage new.smpkg \
  --nucleotide-sdb markers.sdb
```

### Build SingleM packages from scratch (`singlem seqs` → `create` → `regenerate`)
Building a marker package is a multi-step expert workflow:
- **`singlem seqs`** — given an HMM-aligned FASTA, choose the best (most conserved) window position.
- **`singlem create`** — finalise a SingleM package from a GraftM package, a taxonomy file, and the window position from `seqs`.
- **`singlem regenerate`** — update an existing SingleM package with new sequences/taxonomy without rebuilding from scratch.
```bash
# 1. Choose the window position within the HMM
singlem seqs --alignment aligned.fasta --alignment-type aa --hmm marker.hmm

# 2. Create the package using the hmm-position reported by step 1
singlem create \
  --input-graftm-package marker.gpkg \
  --input-taxonomy marker_taxonomy.tsv \
  --hmm-position 25 \
  --target-domains Bacteria Archaea \
  --gene-description "Ribosomal protein S2" \
  --output-singlem-package marker.spkg
```
`--gene-description` is required — it is the free-form text shown by `singlem metapackage --describe`.

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

Lyrebird also provides `condense` and `renew` for archive OTU tables, mirroring their SingleM counterparts but using a Lyrebird metapackage. Save an archive OTU table from `lyrebird pipe` with `--archive-otu-table` to use them:
```bash
# Re-derive a phage profile from an archive OTU table
lyrebird condense \
  --input-archive-otu-table sample.archive.otu_table.json.gz \
  -p sample.phage_profile.tsv

# Re-assign phage taxonomy against an updated Lyrebird metapackage
lyrebird renew \
  --input-archive-otu-table sample.archive.otu_table.json.gz \
  -p sample.updated.phage_profile.tsv \
  --metapackage /path/to/new_lyrebird_metapackage
```

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
