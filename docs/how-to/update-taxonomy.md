# Update taxonomy

Reannotate an archived OTU table with the current reference taxonomy.

```bash
singlem renew --input-archive-otu-table sample.archive.json \
  -p sample.renewed.profile.tsv --otu-table sample.renewed.otu.tsv --threads 4
```

For a Lyrebird archive:

```bash
lyrebird renew --input-archive-otu-table sample.phage-archive.json \
  -p sample.renewed.phage-profile.tsv --threads 4
```
