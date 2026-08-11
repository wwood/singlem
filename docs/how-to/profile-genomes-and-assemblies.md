# Profile genomes and assemblies

Generate OTU tables from genomes or assembled contigs.

```bash
singlem pipe --genome-fasta-files genomes/genome-a.fna genomes/genome-b.fna \
  --otu-table genomes.otu.tsv --threads 4
```

```bash
singlem pipe --genome-fasta-directory genomes \
  --genome-fasta-extension fna --otu-table genomes.otu.tsv --threads 4
```
