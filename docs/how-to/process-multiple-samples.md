# Process multiple samples

Profile multiple paired samples and collect them in one profile.

```bash
singlem pipe -1 sample-a.1.fastq.gz sample-b.1.fastq.gz \
  -2 sample-a.2.fastq.gz sample-b.2.fastq.gz \
  -p samples.profile.tsv --threads 4
```

Combine OTU tables produced by separate runs:

```bash
singlem summarise --input-otu-tables sample-a.otu.tsv sample-b.otu.tsv \
  --output-otu-table samples.otu.tsv
```
