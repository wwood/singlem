**TLDR**: A taxonomic overview of your community can be obtained like so:
```
singlem pipe -1 <fastq_or_fasta1> -2 <fastq_or_fasta2> \
   -p <output.profile.tsv>
```
or, if you have long reads (or single-end short reads):
```
singlem pipe -1 <fastq_or_fasta1> \
   -p <output.profile.tsv>
```
Nanopore reads >= R10.4.1 or PacBio HiFi reads are recommended so reads covering marker gene windows are not missed due to sequencing error.

To further convert the generated taxonomic profile to other formats that might be more convenient, see [`summarise`](/tools/summarise).

## Algorithm details

**Details**: In its most common usage, the SingleM `pipe` subcommand takes as input raw metagenomic reads and outputs a taxonomic profile. Note that taxonomic profiles are generated from OTU tables, they are [not the same thing](/Glossary).

It can also take as input whole genomes (or contigs) via `--genome-fasta-files`, `--genome-fasta-directory` or `--genome-fasta-list`. These options behave the same as providing sequences with `--forward`, except that different defaults for `--min-taxon-coverage` and `--min-orf-length` are used.

`pipe` performs three steps:

1. Find discrete operational taxonomic units (OTUs) from a shotgun metagenome.
2. Assign taxonomy to marker-specific OTU tables.
3. Convert OTU tables into an overall taxonomic profile (otherwise known as a condensed profile). This step happens when the `-p` / `--taxonomic-profile` option of `pipe` is used.

Workflow for the first 2 steps:

![steps 1 and 2](/singlem_pipe_v2.svg)

In the 1st step, reads that encode conserved single copy marker genes are found. SingleM specifically finds reads which cover short highly conserved sections ("*windows*") within those genes. In most species, these windows are 20 amino acids encoded by 60 nucleotides - in rare cases there are inserts or deletions. Sequences covering those small sections are OTU sequences, and these OTU sequences exist independent of taxonomy. By default, SingleM currently uses 35 bacterial and 37 archaeal single copy marker genes.

In the 2nd step, taxonomy is assigned based on comparing the nucleotide sequence of the window to GTDB species representatives' window sequences. If none are similar enough (i.e. within 96.7% identity or 2bp of the 60bp window), then diamond blastx is used instead.

Finally, in the 3rd step, the set of window sequences (i.e. a metagenome's OTU table) is converted into a taxonomic profile, which describes the amount of the microbial community belonging to each species or higher level taxon. This is achieved by considering the OTUs from the 59 different marker genes holistically, using trimmed means and expectation maximisation in a somewhat complicated overall algorithm "*condense*":

![step 3](/singlem_condense_v2.svg)

Please use **raw** metagenomic reads, not quality trimmed reads, if possible. Quality trimming with e.g. [Trimmomatic](https://doi.org/10.1093/bioinformatics/btu170) reads can make them too short for SingleM to use, particularly if they are trimmed to be shorter than 100 bp. Adapter trimming is unlikely to be detrimental, but is not needed.

The [examples section](/tools/pipe#examples) may be of use.

For a more detailed explanation of the SingleM pipeline, see the [SingleM paper](https://doi.org/10.1101/2024.01.30.578060).

