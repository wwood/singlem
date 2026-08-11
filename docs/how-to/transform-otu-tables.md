# Transform OTU tables

Combine OTU tables into one table.

```bash
singlem summarise --input-otu-tables sample-a.otu.tsv sample-b.otu.tsv \
  --output-otu-table combined.otu.tsv
```

Cluster closely related sequences:

```bash
singlem summarise --input-otu-tables combined.otu.tsv --cluster \
  --output-otu-table clustered.otu.tsv
```

Rarefy samples to 100 OTU observations:

```bash
singlem summarise --input-otu-tables combined.otu.tsv \
  --rarefied-output-otu-table rarefied.otu.tsv --number-to-choose 100
```
