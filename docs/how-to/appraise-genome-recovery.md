# Appraise genome recovery

Measure how much of a metagenomic community is represented by recovered genomes.

```bash
singlem pipe -1 reads.fastq.gz --otu-table metagenome.otu.tsv --threads 4
singlem pipe --genome-fasta-directory genomes --genome-fasta-extension fna \
  --otu-table genomes.otu.tsv --threads 4
singlem appraise --metagenome-otu-tables metagenome.otu.tsv \
  --genome-otu-tables genomes.otu.tsv
```
