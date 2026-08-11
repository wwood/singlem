# Profile a small genome with SingleM

This tutorial downloads a small example genome, runs `pipe`, and reads its taxonomic profile.

Create a working directory and download the example:

```bash
mkdir singlem-tutorial
cd singlem-tutorial
wget https://github.com/wwood/singlem/raw/44e1f81404c12931742259088999290edbb271b3/test/data/methanobacteria/genomes/GCA_000309865.1_genomic.fna
```

Run `pipe`:

```bash
singlem pipe -1 GCA_000309865.1_genomic.fna -p example.profile.tsv --threads 4
```

Read the result:

```bash
column -t -s $'\t' example.profile.tsv
```

Each row reports the sample, estimated coverage, and a taxonomic lineage. The most specific populated lineage identifies the example as a member of the genus `Methanobacterium`.
