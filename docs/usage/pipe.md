---
title: SingleM pipe
---
# singlem pipe

The SingleM `pipe` subcommand performs three steps:

1. Finding discrete operational taxonomic units (OTUs) from a shotgun metagenome
2. Assign taxonomy to marker-specific OTU tables
3. Convert OTU tables into a overall taxonomic profile

In the first step, reads that encode conserved single copy marker genes are found. SingleM specifically finds reads which cover short (20 amino acid / 60 base pair) highly conserved sections ("windows") within those genes. Sequences covering those small sections are OTU sequences, and these OTU sequences exist independent of taxonomy. By default, SingleM currently uses 35 bacterial and 37 archaeal single copy marker genes.

In the second step, taxonomy is assigned based on comparing the nucleotide sequence of the window to GTDB species representatives' window sequences. If none are similar enough (i.e. within 95% identity or 3bp of the 60), then diamond blastx is used.

A common analysis is to convert a set of window sequences (i.e. a metagenome's OTU table) into a taxonomic community profile. This is achieved in the second step by considering the OTUs from the 59 different marker genes holistically, using trimmed means and expectation maximisation.

An overview of your community can be obtained like so:
```
singlem pipe -1 <fastq_or_fasta1> -2 <fastq_or_fasta2> -p \
   <output.profile.tsv>
```
Please use **raw** metagenomic reads, not quality trimmed reads. Quality trimming with e.g. [Trimmomatic](https://doi.org/10.1093/bioinformatics/btu170) reads often makes them too short for SingleM to use. Adapter trimming is unlikely to be detrimental, but is not needed.

The [examples section](/usage/pipe#examples) may be of use.

COMMON OPTIONS
==============

**-1**, **\--forward**, **\--reads**, **\--sequences** sequence_file [sequence_file \...]

  nucleotide read sequence(s) (forward or unpaired) to be searched.
    Can be FASTA or FASTQ format, GZIP-compressed or not.

**-2**, **\--reverse** sequence_file [sequence_file \...]

  reverse reads to be searched. Can be FASTA or FASTQ format,
    GZIP-compressed or not.

**\--genome-fasta-files** sequence_file [sequence_file \...]

  nucleotide genome sequence(s) to be searched

**\--sra-files** sra_file [sra_file \...]

  \"sra\" format files (usually from NCBI SRA) to be searched

**-p**, **\--taxonomic-profile** FILE

  output a \'condensed\' taxonomic profile for each sample based on
    the OTU table

**\--taxonomic-profile-krona** FILE

  output a \'condensed\' taxonomic profile for each sample based on
    the OTU table

**\--otu-table** filename

  output OTU table

**\--threads** num_threads

  number of CPUS to use [default: 1]

**\--assignment-method** {smafa_naive_then_diamond,scann_naive_then_diamond,annoy_then_diamond,scann_then_diamond,diamond,diamond_example,annoy,pplacer}

  Method of assigning taxonomy to OTUs and taxonomic profiles
    [default: smafa_naive_then_diamond]

| Method                   | Description                                                                                                                                                                                                                                   |
|:-------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| smafa_naive_then_diamond | Search for the most similar window sequences \<= 3bp different using a brute force algorithm (using the smafa implementation) over all window sequences in the database, and if none are found use DIAMOND blastx of all reads from each OTU. |
| scann_naive_then_diamond | Search for the most similar window sequences \<= 3bp different using a brute force algorithm over all window sequences in the database, and if none are found use DIAMOND blastx of all reads from each OTU.                                  |
| annoy_then_diamond       | Same as scann_naive_then_diamond, except search using ANNOY rather than using brute force. Requires a non-standard metapackage.                                                                                                               |
| scann_then_diamond       | Same as scann_naive_then_diamond, except search using SCANN rather than using brute force. Requires a non-standard metapackage.                                                                                                               |
| diamond                  | DIAMOND blastx best hit(s) of all reads from each OTU.                                                                                                                                                                                        |
| diamond_example          | DIAMOND blastx best hit(s) of all reads from each OTU, but report the best hit as a sequence ID instead of a taxonomy.                                                                                                                        |
| annoy                    | Search for the most similar window sequences \<= 3bp different using ANNOY, otherwise no taxonomy is assigned. Requires a non-standard metapackage.                                                                                           |
| pplacer                  | Use pplacer to assign taxonomy of each read in each OTU. Requires a non-standard metapackage.                                                                                                                                                 |

**\--output-extras**

  give extra output for each sequence identified (e.g. the read(s)
    each OTU was generated from) in the output OTU table [default: not
    set]

LESS COMMON OPTIONS
===================

**\--archive-otu-table** filename

  output OTU table in archive format for making DBs etc. [default:
    unused]

**\--output-jplace** filename

  Output a jplace format file for each singlem package to a file
    starting with this string, each with one entry per OTU. Requires
    \'pplacer\' as the \--assignment_method [default: unused]

**\--metapackage** *METAPACKAGE*

  Set of SingleM packages to use [default: use the default set]

**\--singlem-packages** *SINGLEM_PACKAGES* [*SINGLEM_PACKAGES* \...]

  SingleM packages to use [default: use the set from the default
    metapackage]

**\--assignment-singlem-db** *ASSIGNMENT_SINGLEM_DB*

  Use this SingleM DB when assigning taxonomy [default: not set, use
    the default]

**\--diamond-taxonomy-assignment-performance-parameters** *DIAMOND_TAXONOMY_ASSIGNMENT_PERFORMANCE_PARAMETERS*

  Performance-type arguments to use when calling \'diamond blastx\'
    during the taxonomy assignment step. [default: \'\--block-size 0.5
    \--target-indexed -c1\']

**\--evalue** *EVALUE*

  GraftM e-value cutoff [default: the GraftM default]

**\--min-orf-length** length

  When predicting ORFs require this many base pairs uninterrupted by a
    stop codon [default: 72 for reads, 300 for genomes]

**\--restrict-read-length** length

  Only use this many base pairs at the start of each sequence searched
    [default: no restriction]

**\--filter-minimum-protein** length

  Ignore reads aligning in less than this many positions to each
    protein HMM [default: 24]

**\--exclude-off-target-hits**

  Exclude hits that are not in the target domain of each SingleM
    package

**\--min-taxon-coverage** FLOAT

  Minimum coverage to report in a taxonomic profile. [default: 0.35
    for reads, 0.1 for genomes]

**\--working-directory** directory

  use intermediate working directory at a specified location, and do
    not delete it upon completion [default: not set, use a temporary
    directory]

**\--working-directory-dev-shm**

  use an intermediate results temporary working directory in /dev/shm
    rather than the default [default: the usual temporary working
    directory, currently /tmp]

**\--force**

  overwrite working directory if required [default: not set]

**\--filter-minimum-nucleotide** length

  Ignore reads aligning in less than this many positions to each
    nucleotide HMM [default: 72]

**\--include-inserts**

  print the entirety of the sequences in the OTU table, not just the
    aligned nucleotides [default: not set]

**\--known-otu-tables** *KNOWN_OTU_TABLES* [*KNOWN_OTU_TABLES* \...]

  OTU tables previously generated with trusted taxonomies for each
    sequence [default: unused]

**\--no-assign-taxonomy**

  Do not assign any taxonomy except for those already known [default:
    not set]

**\--known-sequence-taxonomy** FILE

  A 2-column \"sequence\<tab\>taxonomy\" file specifying some
    sequences that have known taxonomy [default: unused]

**\--no-diamond-prefilter**

  Do not parse sequence data through DIAMOND blastx using a database
    constructed from the set of singlem packages. Should be used with
    \--hmmsearch-package-assignment. NOTE: ignored for nucleotide
    packages [default: protein packages: use the prefilter, nucleotide
    packages: do not use the prefilter]

**\--diamond-prefilter-performance-parameters** *DIAMOND_PREFILTER_PERFORMANCE_PARAMETERS*

  Performance-type arguments to use when calling \'diamond blastx\'
    during the prefiltering. By default, SingleM should run in \<4GB of
    RAM except in very large (\>100Gbp) metagenomes. [default:
    \'\--block-size 0.5 \--target-indexed -c1\']

**\--hmmsearch-package-assignment**

  Assign each sequence to a SingleM package using HMMSEARCH, and a
    sequence may then be assigned to multiple packages. [default: not
    set]

**\--diamond-prefilter-db** *DIAMOND_PREFILTER_DB*

  Use this DB when running DIAMOND prefilter [default: use the one in
    the metapackage, or generate one from the SingleM packages]

**\--assignment-threads** *ASSIGNMENT_THREADS*

  Use this many processes in parallel while assigning taxonomy
    [default: 1]

**\--sleep-after-mkfifo** *SLEEP_AFTER_MKFIFO*

  Sleep for this many seconds after running os.mkfifo [default:
    None]

OTHER GENERAL OPTIONS
=====================

**\--debug**

  output debug information

**\--version**

  output version information and quit

**\--quiet**

  only output errors

**\--full-help**

  print longer help message

**\--full-help-roff**

  print longer help message in ROFF (manpage) format

AUTHORS
=======

>     Ben J. Woodcroft, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
>     Samuel Aroney, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
>     Rossen Zhao, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology

EXAMPLES
========

Get a taxonomic profile from paired read input:

  **\$ singlem pipe -1 \<fastq_or_fasta1\> -2 \<fastq_or_fasta2\> -p
    \<output.profile.tsv\>**

Get a taxonomic profile Krona diagram from single read input:

  **\$ singlem pipe -i \<fastq_or_fasta\> \--taxonomic-profile-krona
    \<output.profile.html\>**

Gather an OTU table (per marker sequence groupings) from paired reads:

  **\$ singlem pipe -1 \<fastq_or_fasta1\> -2 \<fastq_or_fasta2\>
    \--otu-table \<output.otu_table.tsv\>**
