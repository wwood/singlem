# Estimate the prokaryotic fraction

Estimate the bacterial and archaeal fraction of a metagenome.

```bash
singlem pipe -1 sample.1.fastq.gz -2 sample.2.fastq.gz \
  -p sample.profile.tsv --threads 4
singlem prokaryotic_fraction -1 sample.1.fastq.gz -2 sample.2.fastq.gz \
  -p sample.profile.tsv --output-tsv sample.spf.tsv
```

Read the `read_fraction` column as a percentage. See [what SPF means](../explanation/spf.md) for its assumptions and limitations.
