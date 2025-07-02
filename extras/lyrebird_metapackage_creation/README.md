Update config.yaml and other relevant files, then run each pipeline as follows:
```
snakemake --snakefile <pipeline> --configfile config.yaml --threads 64 --use-conda
```
In short, run in this order: 
1. 1_initial_filtering.smk
2. 2_concatenate_and_build.smk
3. 3_build_mpkg.smk

To generate the gtdb_proviruses file, run "genomad find-proviruses" on all GTDB genome representatives