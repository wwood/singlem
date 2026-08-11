# Add genomes to a metapackage

Create a metapackage containing additional genomes.

```bash
singlem supplement --new-genome-fasta-files genomes/genome-a.fna genomes/genome-b.fna \
  --output-metapackage supplemented.smpkg --threads 4
```

If taxonomy is inferred with GTDB-Tk, use a GTDB-Tk release compatible with the input metapackage.
