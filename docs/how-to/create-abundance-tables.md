# Create abundance tables

Convert taxonomic profiles into a genus-by-sample relative-abundance table.

```bash
singlem summarise --input-taxonomic-profiles sample-a.profile.tsv sample-b.profile.tsv \
  --output-species-by-site-relative-abundance abundance.tsv \
  --output-species-by-site-level genus
```

Add filled coverage and relative abundance to a long-form profile:

```bash
singlem summarise --input-taxonomic-profiles sample-a.profile.tsv \
  --output-taxonomic-profile-with-extras sample-a.extras.tsv
```

See [why profiles use condense](../explanation/condense.md).
