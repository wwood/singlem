---
title: renew
---
# singlem renew

DESCRIPTION
===========

Reannotate an OTU table with an updated taxonomy

INPUT
=====

**\--input-archive-otu-table** *INPUT_ARCHIVE_OTU_TABLE*

  Renew this table

COMMON ARGUMENTS IN SHARED WITH \'PIPE\'
========================================

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

**\--assignment-method** {naive_then_diamond,annoy_then_diamond,scann_then_diamond,diamond,diamond_example,annoy,pplacer}

  Method of assigning taxonomy to OTUs and taxonomic profiles
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

  give extra output for each sequence identified (e.g. the read(s)
    each OTU was generated from) in the output OTU table [default: not
    set]

LESS COMMON ARGUMENTS SHARED WITH \'PIPE\'
==========================================

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
