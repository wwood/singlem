---
title: "singlem summarise"
author: "Ben Woodcroft, Centre for Microbiome Research, Queensland University of Technology"
date: "`r Sys.Date()` (`r system('bin/singlem --version',intern=T)`)"
---
NAME
====

singlem summarise

SYNOPSIS
========

**singlem** summarise [-h] [\--input-otu-tables INPUT_OTU_TABLES
[INPUT_OTU_TABLES \...]] [\--input-otu-tables-list
INPUT_OTU_TABLES_LIST] [\--input-archive-otu-tables
INPUT_ARCHIVE_OTU_TABLES [INPUT_ARCHIVE_OTU_TABLES \...]]
[\--input-gzip-archive-otu-table-list
INPUT_GZIP_ARCHIVE_OTU_TABLE_LIST] [\--input-archive-otu-table-list
INPUT_ARCHIVE_OTU_TABLE_LIST] [\--stream-inputs] [\--cluster]
[\--cluster-id CLUSTER_ID] [\--taxonomy TAXONOMY]
[\--rarefied-output-otu-table RAREFIED_OUTPUT_OTU_TABLE]
[\--number-to-choose NUMBER_TO_CHOOSE] [\--collapse-coupled]
[\--collapse-paired-with-unpaired-archive-otu-table
COLLAPSE_PAIRED_WITH_UNPAIRED_ARCHIVE_OTU_TABLE] [\--output-otu-table
OUTPUT_OTU_TABLE] [\--output-translated-otu-table
OUTPUT_TRANSLATED_OTU_TABLE] [\--output-extras] [\--krona KRONA]
[\--wide-format-otu-table WIDE_FORMAT_OTU_TABLE]
[\--strain-overview-table STRAIN_OVERVIEW_TABLE] [\--unifrac-by-otu
UNIFRAC_BY_OTU] [\--unifrac-by-taxonomy UNIFRAC_BY_TAXONOMY]
[\--clustered-output-otu-table CLUSTERED_OUTPUT_OTU_TABLE]
[\--exclude-off-target-hits] [\--singlem-packages SINGLEM_PACKAGES
[SINGLEM_PACKAGES \...]] [\--metapackage METAPACKAGE]
[\--unaligned-sequences-dump-file UNALIGNED_SEQUENCES_DUMP_FILE]
[\--debug] [\--version] [\--quiet] [\--full-help]
[\--full-help-roff]

DESCRIPTION
===========

Summarise and transform OTU tables.

INPUT
=====

**\--input-otu-tables**, **\--input-otu-table** *INPUT_OTU_TABLES* [*INPUT_OTU_TABLES* \...]

:   Summarise these tables

**\--input-otu-tables-list** *INPUT_OTU_TABLES_LIST*

:   Summarise the OTU table files newline separated in this file

**\--input-archive-otu-tables**, **\--input-archive-otu-table** *INPUT_ARCHIVE_OTU_TABLES* [*INPUT_ARCHIVE_OTU_TABLES* \...]

:   Summarise these tables

**\--input-gzip-archive-otu-table-list** *INPUT_GZIP_ARCHIVE_OTU_TABLE_LIST*

:   Summarise the list of newline-separated gzip-compressed archive OTU
    tables specified in this file

**\--input-archive-otu-table-list** *INPUT_ARCHIVE_OTU_TABLE_LIST*

:   Summarise the archive tables newline separated in this file

**\--stream-inputs**

:   Stream input OTU tables, saving RAM. Only works with
    \--output-otu-table and transformation options do not work [expert
    option].

TRANSFORMATION
==============

**\--cluster**

:   Apply sequence clustering to the OTU table

**\--cluster-id** *CLUSTER_ID*

:   Sequence clustering identity cutoff if \--cluster is used

**\--taxonomy** *TAXONOMY*

:   Restrict analysis to OTUs that have this taxonomy (exact taxonomy or
    more fully resolved)

**\--rarefied-output-otu-table** *RAREFIED_OUTPUT_OTU_TABLE*

:   Output rarefied output OTU table, where each gene and sample
    combination is rarefied

**\--number-to-choose** *NUMBER_TO_CHOOSE*

:   Rarefy using this many sequences. Sample/gene combinations with an
    insufficient number of sequences are ignored with a warning
    [default: maximal number such that all samples have sufficient
    counts]

**\--collapse-coupled**

:   Merge forward and reverse read OTU tables into a unified table.
    Sample names of coupled reads must end in \'1\' and \'2\'
    respectively. Read names are ignored, so that if the forward and
    reverse from a pair contain the same OTU sequence, they will each
    count separately.

**\--collapse-paired-with-unpaired-archive-otu-table** *COLLAPSE_PAIRED_WITH_UNPAIRED_ARCHIVE_OTU_TABLE*

:   For archive OTU tables that have both paired and unpaired
    components, merge these into a single output archive OTU table

OUTPUT
======

**\--output-otu-table** *OUTPUT_OTU_TABLE*

:   Output combined OTU table to this file

**\--output-translated-otu-table** *OUTPUT_TRANSLATED_OTU_TABLE*

:   Output combined OTU table to this file, with seqeunces translated
    into amino acids

**\--output-extras**

:   Output extra information in the standard output OTU table

**\--krona** *KRONA*

:   Name of krona file to generate

**\--wide-format-otu-table** *WIDE_FORMAT_OTU_TABLE*

:   Name of output species by site CSV file

**\--strain-overview-table** *STRAIN_OVERVIEW_TABLE*

:   Name of output strains table to generate

**\--unifrac-by-otu** *UNIFRAC_BY_OTU*

:   Output UniFrac format file where entries are OTU sequences

**\--unifrac-by-taxonomy** *UNIFRAC_BY_TAXONOMY*

:   Output UniFrac format file where entries are taxonomies (generally
    used for phylogeny-driven beta diversity when pipe was run with
    \'\--assignment_method diamond_example\')

**\--clustered-output-otu-table** *CLUSTERED_OUTPUT_OTU_TABLE*

:   Output an OTU table with extra information about the clusters

**\--exclude-off-target-hits**

:   Exclude hits that are not in the target domain of each SingleM
    package

**\--singlem-packages** *SINGLEM_PACKAGES* [*SINGLEM_PACKAGES* \...]

:   Packages used in the creation of the OTU tables

**\--metapackage** *METAPACKAGE*

:   Metapackage used in the creation of the OTU tables

**\--unaligned-sequences-dump-file** *UNALIGNED_SEQUENCES_DUMP_FILE*

:   Output unaligned sequences from in put archive OTU table to this
    file. After each read name \'\~N\' is added which corresponds to the
    order of the read in the archive OTU table, so that no two sequences
    have the same read name. N\>1 can happen e.g. when the input file
    contains paired reads, but \~0 does not necessarily correspond to
    the first read in the original input sequence set, but instead to
    the order in the input archive OTU table.

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
