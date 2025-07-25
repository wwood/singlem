Update config.yaml and other relevant files, then run each pipeline as follows:
```
snakemake --snakefile <pipeline> --configfile config.yaml --threads 64 --use-conda
```
In short, run in this order: 
1. 1_initial_filtering.smk
2. 2_concatenate_and_build.smk
3. 3_build_mpkg.smk

To generate viral proteins and transcripts, run prodigal-gv on all viral genome representatives.

The viral taxonomy file is generated as a tab-separated genome ID to taxonomy list:
```
{Genome ID} d__{domain};p__{phylum};c__{class};o__{order};f__{family};g__{genus};s__{species}
```

To generate the gtdb_proviruses file, run genomad end-to-end on all GTDB genome representatives, and concatenate the resulting {genome}virus_summary.tsv file. This is required for masking proviral regions and viral contamination for off-target sequences. 