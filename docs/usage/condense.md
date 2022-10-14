---
title: SingleM condense
---
# singlem condense

DESCRIPTION
===========

Combine OTU tables across different markers into a single taxonomic
profile.

INPUT ARGUMENTS (1+ REQUIRED)
=============================

**\--input-otu-tables**, **\--input-otu-table** *INPUT_OTU_TABLES* [*INPUT_OTU_TABLES* \...]

  Condense these table(s).

**\--input-otu-tables-list** *INPUT_OTU_TABLES_LIST*

  Condense the OTU table files newline separated in this file

**\--input-archive-otu-tables**, **\--input-archive-otu-table** *INPUT_ARCHIVE_OTU_TABLES* [*INPUT_ARCHIVE_OTU_TABLES* \...]

  Condense from these archive tables

**\--input-archive-otu-table-list** *INPUT_ARCHIVE_OTU_TABLE_LIST*

  Condense from the archive tables newline separated in this file

**\--input-gzip-archive-otu-table-list** *INPUT_GZIP_ARCHIVE_OTU_TABLE_LIST*

  Condense from the gzip\'d archive tables newline separated in this
    file

OUTPUT ARGUMENTS (1+ REQUIRED)
==============================

**-p**, **\--taxonomic-profile** filename

  output OTU table

**\--taxonomic-profile-krona** filename

  name of krona file to generate.

**\--output-after-em-otu-table** filename

  output OTU table after expectation maximisation has been applied.
    Note that this table usually contains multiple rows with the same
    window sequence.

OTHER OPTIONS
=============

**\--metapackage** *METAPACKAGE*

  Set of SingleM packages to use [default: use the default set]

**\--min-taxon-coverage** FRACTION

  Set taxons with less coverage to coverage=0. [default: 0.05]

**\--trim-percent** *TRIM_PERCENT*

  percentage of markers to be trimmed for each taxonomy [default:
    10]

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
