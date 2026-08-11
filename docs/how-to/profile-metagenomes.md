# Profile metagenomes

Generate a microbial taxonomic profile from paired metagenomic reads.

```bash
singlem pipe -1 sample.1.fastq.gz -2 sample.2.fastq.gz \
  -p sample.profile.tsv --archive-otu-table sample.archive.json --threads 4
```

For single-end or long reads:

```bash
singlem pipe -1 sample.fastq.gz \
  -p sample.profile.tsv --archive-otu-table sample.archive.json --threads 4
```

Nanopore R10.4.1 or newer and PacBio HiFi reads are recommended.
