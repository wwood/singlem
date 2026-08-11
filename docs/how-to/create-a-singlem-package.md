# Create a SingleM package

Choose a marker window and build a SingleM package from an existing GraftM package.

```bash
hmmalign marker.hmm marker-sequences.faa > marker.sto
seqmagick convert marker.sto marker.fasta
singlem seqs --alignment marker.fasta
singlem create --input-graftm-package marker.gpkg \
  --input-taxonomy taxonomy.tsv --output-singlem-package marker.spkg \
  --hmm-position 42 --target-domains Bacteria \
  --gene-description "Marker gene"
```

Use the position reported by `seqs` in place of `42`. See [the marker window](../explanation/marker-window.md).
