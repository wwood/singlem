---
title: makedb
---
# singlem makedb

DESCRIPTION
===========

Create a searchable database from an OTU table

REQUIRED ARGUMENTS
==================

**\--otu-tables**, **\--otu-table** *OTU_TABLES* [*OTU_TABLES* \...]

  Make a db from these OTU tables

**\--otu-tables-list** *OTU_TABLES_LIST*

  Make a db from the OTU table files newline separated in this file

**\--archive-otu-tables**, **\--archive-otu-table** *ARCHIVE_OTU_TABLES* [*ARCHIVE_OTU_TABLES* \...]

  Make a db from these archive tables

**\--archive-otu-table-list** *ARCHIVE_OTU_TABLE_LIST*

  Make a db from the archive tables newline separated in this file

**\--gzip-archive-otu-table-list** *GZIP_ARCHIVE_OTU_TABLE_LIST*

  Make a db from the gzip\'d archive tables newline separated in this
    file

**\--db** *DB*

  Name of database to create e.g. tundra.sdb

OTHER ARGUMENTS
===============

**\--threads** *THREADS*

  Use this many threads where possible [default 1]

**\--sequence-database-methods** {annoy,scann,nmslib,naive,none} [{annoy,scann,nmslib,naive,none} \...]

  Index sequences using these methods. Note that specifying \"naive\"
    means \"scann\" databases will also be built [default
    [\'naive\']]

**\--sequence-database-types** {nucleotide,protein} [{nucleotide,protein} \...]

  Index sequences using these types. [default: [\'nucleotide\']]

**\--pregenerated-otu-sqlite-db** *PREGENERATED_OTU_SQLITE_DB*

  [for internal usage] remake the indices using this input SQLite
    database

**\--num-annoy-nucleotide-trees** *NUM_ANNOY_NUCLEOTIDE_TREES*

  make annoy nucleotide sequence indices with this ntrees [default
    10]

**\--num-annoy-protein-trees** *NUM_ANNOY_PROTEIN_TREES*

  make annoy protein sequence indices with this ntrees [default 10]

**\--tmpdir** *TMPDIR*

  [for internal usage] use this directory internally for working

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
