# Use a 16S package

Profile reads with a separately downloaded 16S SingleM package.

Download a suitable package from the [SingleM auxiliary package repository](https://github.com/wwood/singlem_extra_packages), then run:

```bash
singlem pipe -1 sample.fastq.gz --singlem-packages 16s.spkg \
  --otu-table sample.16s.otu.tsv --threads 4
```
