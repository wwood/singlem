# Profile phages

Generate a dsDNA phage taxonomic profile from paired reads.

```bash
lyrebird pipe -1 sample.1.fastq.gz -2 sample.2.fastq.gz \
  -p sample.phage-profile.tsv --archive-otu-table sample.phage-archive.json \
  --threads 4
```
