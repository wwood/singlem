---
title: "singlem pipe"
author: "Ben Woodcroft, Centre for Microbiome Research, Queensland University of Technology"
date: "`r Sys.Date()` (`r system('bin/singlem --version',intern=T)`)"
---
NAME
====

singlem pipe

SYNOPSIS
========

**singlem** pipe [-h] [-1 sequence_file [sequence_file \...]] [-2
sequence_file [sequence_file \...]] [\--genome-fasta-files
sequence_file [sequence_file \...]] [\--sra-files sra_file
[sra_file \...]] [-p FILE] [\--taxonomic-profile-krona FILE]
[\--otu-table filename] [\--threads num_threads]
[\--assignment-method
{naive_then_diamond,annoy_then_diamond,scann_then_diamond,diamond,diamond_example,annoy,pplacer}]
[\--output-extras] [\--archive-otu-table filename]
[\--output-jplace filename] [\--metapackage METAPACKAGE]
[\--singlem-packages SINGLEM_PACKAGES [SINGLEM_PACKAGES \...]]
[\--assignment-singlem-db ASSIGNMENT_SINGLEM_DB]
[\--diamond-taxonomy-assignment-performance-parameters
DIAMOND_TAXONOMY_ASSIGNMENT_PERFORMANCE_PARAMETERS] [\--evalue
EVALUE] [\--min-orf-length length] [\--restrict-read-length length]
[\--filter-minimum-protein length] [\--working-directory directory]
[\--working-directory-dev-shm] [\--force]
[\--filter-minimum-nucleotide length] [\--include-inserts]
[\--known-otu-tables KNOWN_OTU_TABLES [KNOWN_OTU_TABLES \...]]
[\--no-assign-taxonomy] [\--known-sequence-taxonomy FILE]
[\--no-diamond-prefilter]
[\--diamond-prefilter-performance-parameters
DIAMOND_PREFILTER_PERFORMANCE_PARAMETERS]
[\--hmmsearch-package-assignment] [\--diamond-prefilter-db
DIAMOND_PREFILTER_DB] [\--assignment-threads ASSIGNMENT_THREADS]
[\--sleep-after-mkfifo SLEEP_AFTER_MKFIFO] [\--debug] [\--version]
[\--quiet] [\--full-help] [\--full-help-roff]

DESCRIPTION
===========

Generate a taxonomic profile or OTU table from raw sequences

COMMON OPTIONS
==============

**-1**, **\--forward**, **\--reads**, **\--sequences** sequence_file [sequence_file \...]

:   nucleotide read sequence(s) (forward or unpaired) to be searched.
    Can be FASTA or FASTQ format, GZIP-compressed or not.

**-2**, **\--reverse** sequence_file [sequence_file \...]

:   reverse reads to be searched. Can be FASTA or FASTQ format,
    GZIP-compressed or not.

**\--genome-fasta-files** sequence_file [sequence_file \...]

:   nucleotide genome sequence(s) to be searched

**\--sra-files** sra_file [sra_file \...]

:   \"sra\" format files (usually from NCBI SRA) to be searched

**-p**, **\--taxonomic-profile** FILE

:   output a \'condensed\' taxonomic profile for each sample based on
    the OTU table

**\--taxonomic-profile-krona** FILE

:   output a \'condensed\' taxonomic profile for each sample based on
    the OTU table

**\--otu-table** filename

:   output OTU table

**\--threads** num_threads

:   number of CPUS to use [default: 1]

**\--assignment-method** {naive_then_diamond,annoy_then_diamond,scann_then_diamond,diamond,diamond_example,annoy,pplacer}

:   Method of assigning taxonomy to OTUs and taxonomic profiles
    [default: naive_then_diamond]

| Method             | Description                                                                                                                                                                                                  |
|:-------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| naive_then_diamond | Search for the most similar window sequences \<= 3bp different using a brute force algorithm over all window sequences in the database, and if none are found use DIAMOND blastx of all reads from each OTU. |
| annoy_then_diamond | Same as naive_then_diamond, except search using ANNOY rather than using brute force. Requires a non-standard metapackage.                                                                                    |
| scann_then_diamond | Same as naive_then_diamond, except search using SCANN rather than using brute force. Requires a non-standard metapackage.                                                                                    |
| diamond            | DIAMOND blastx best hit(s) of all reads from each OTU.                                                                                                                                                       |
| diamond_example    | DIAMOND blastx best hit(s) of all reads from each OTU, but report the best hit as a sequence ID instead of a taxonomy.                                                                                       |
| annoy              | Search for the most similar window sequences \<= 3bp different using ANNOY, otherwise no taxonomy is assigned. Requires a non-standard metapackage.                                                          |
| pplacer            | Use pplacer to assign taxonomy of each read in each OTU. Requires a non-standard metapackage.                                                                                                                |

**\--output-extras**

:   give extra output for each sequence identified (e.g. the read(s)
    each OTU was generated from) in the output OTU table [default: not
    set]

LESS COMMON OPTIONS
===================

**\--archive-otu-table** filename

:   output OTU table in archive format for making DBs etc. [default:
    unused]

**\--output-jplace** filename

:   Output a jplace format file for each singlem package to a file
    starting with this string, each with one entry per OTU. Requires
    \'pplacer\' as the \--assignment_method [default: unused]

**\--metapackage** *METAPACKAGE*

:   Set of SingleM packages to use [default: use the default set]

**\--singlem-packages** *SINGLEM_PACKAGES* [*SINGLEM_PACKAGES* \...]

:   SingleM packages to use [default: use the set from the default
    metapackage]

**\--assignment-singlem-db** *ASSIGNMENT_SINGLEM_DB*

:   Use this SingleM DB when assigning taxonomy [default: not set, use
    the default]

**\--diamond-taxonomy-assignment-performance-parameters** *DIAMOND_TAXONOMY_ASSIGNMENT_PERFORMANCE_PARAMETERS*

:   Performance-type arguments to use when calling \'diamond blastx\'
    during the taxonomy assignment step. [default: \'\--block-size 0.5
    \--target-indexed -c1\']

**\--evalue** *EVALUE*

:   GraftM e-value cutoff [default: the GraftM default]

**\--min-orf-length** length

:   When predicting ORFs require this many base pairs uninterrupted by a
    stop codon [default: 72 for reads, 300 for genomes]

**\--restrict-read-length** length

:   Only use this many base pairs at the start of each sequence searched
    [default: no restriction]

**\--filter-minimum-protein** length

:   Ignore reads aligning in less than this many positions to each
    protein HMM [default: 24]

**\--working-directory** directory

:   use intermediate working directory at a specified location, and do
    not delete it upon completion [default: not set, use a temporary
    directory]

**\--working-directory-dev-shm**

:   use an intermediate results temporary working directory in /dev/shm
    rather than the default [default: the usual temporary working
    directory, currently /tmp]

**\--force**

:   overwrite working directory if required [default: not set]

**\--filter-minimum-nucleotide** length

:   Ignore reads aligning in less than this many positions to each
    nucleotide HMM [default: 72]

**\--include-inserts**

:   print the entirety of the sequences in the OTU table, not just the
    aligned nucleotides [default: not set]

**\--known-otu-tables** *KNOWN_OTU_TABLES* [*KNOWN_OTU_TABLES* \...]

:   OTU tables previously generated with trusted taxonomies for each
    sequence [default: unused]

**\--no-assign-taxonomy**

:   Do not assign any taxonomy except for those already known [default:
    not set]

**\--known-sequence-taxonomy** FILE

:   A 2-column \"sequence\<tab\>taxonomy\" file specifying some
    sequences that have known taxonomy [default: unused]

**\--no-diamond-prefilter**

:   Do not parse sequence data through DIAMOND blastx using a database
    constructed from the set of singlem packages. Should be used with
    \--hmmsearch-package-assignment. NOTE: ignored for nucleotide
    packages [default: protein packages: use the prefilter, nucleotide
    packages: do not use the prefilter]

**\--diamond-prefilter-performance-parameters** *DIAMOND_PREFILTER_PERFORMANCE_PARAMETERS*

:   Performance-type arguments to use when calling \'diamond blastx\'
    during the prefiltering. By default, SingleM should run in \<4GB of
    RAM except in very large (\>100Gbp) metagenomes. [default:
    \'\--block-size 0.5 \--target-indexed -c1\']

**\--hmmsearch-package-assignment**

:   Assign each sequence to a SingleM package using HMMSEARCH, and a
    sequence may then be assigned to multiple packages. [default: not
    set]

**\--diamond-prefilter-db** *DIAMOND_PREFILTER_DB*

:   Use this DB when running DIAMOND prefilter [default: not set,
    generate one from the SingleM packages]

**\--assignment-threads** *ASSIGNMENT_THREADS*

:   Use this many processes in parallel while assigning taxonomy
    [default: 1]

**\--sleep-after-mkfifo** *SLEEP_AFTER_MKFIFO*

:   Sleep for this many seconds after running os.mkfifo [default:
    None]

OTHER GENERAL OPTIONS
=====================

**\--debug**

:   output debug information

**\--version**

:   output version information and quit

**\--quiet**

:   only output errors

**\--full-help**

:   print longer help message

**\--full-help-roff**

:   print longer help message in ROFF (manpage) format

AUTHORS
=======

>     Ben J. Woodcroft, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
>     Samuel Aroney, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
>     Rossen Zhao, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
