---
title: SingleM query
---
# singlem query

# DESCRIPTION

Find closely related sequences in a SingleM database.

# OPTIONS

# REQUIRED ARGUMENTS

**\--db** *DB*

  Output from \'makedb\' mode

# DATABASE QUERYING BY OTU SEQUENCE

**\--query-otu-table**, **\--query-otu-tables** file [file \...]

  Query the database with all sequences in this OTU table

**\--query-otu-tables-list** *QUERY_OTU_TABLES_LIST*

  Query the database with all sequences in OTU table files newline
    separated in this file

**\--query-archive-otu-tables** *QUERY_ARCHIVE_OTU_TABLES* [*QUERY_ARCHIVE_OTU_TABLES* \...]

  Query the database with all sequences in these archive tables

**\--query-archive-otu-table-list** *QUERY_ARCHIVE_OTU_TABLE_LIST*

  Query the database with all sequences in archive tables newline
    separated in this file

**\--query-gzip-archive-otu-table-list** *QUERY_GZIP_ARCHIVE_OTU_TABLE_LIST*

  Query the database with all sequences in gzip\'d archive tables
    newline separated in this file

**\--max-nearest-neighbours** *MAX_NEAREST_NEIGHBOURS*

  How many nearest neighbours to report. Each neighbour is a distinct
    sequence from the DB. [default: 20]

**\--max-divergence** INT

  Report sequences less than or equal to this divergence i.e. number
    of different bases/amino acids

**\--search-method** {smafa-naive,nmslib,annoy,scann,scann-naive}

  Algorithm to perform search [default: smafa-naive]

**\--sequence-type** {nucleotide,protein}

  Which sequence types to compare (i.e. protein for blastp, nucleotide
    for blastn) [default: nucleotide]

**\--max-search-nearest-neighbours** *MAX_SEARCH_NEAREST_NEIGHBOURS*

  How many nearest neighbours to search for with approximate nearest
    neighbours. Of these hits, only \--max-nearest-neighbours will
    actually be reported. Ignored for \--search-method naive and
    scann-naive. [default: 100]

**\--threads** *THREADS*

  Use this many threads where possible [default 1]

**\--limit-per-sequence** *LIMIT_PER_SEQUENCE*

  How many entries (samples/genomes from DB with identical sequences)
    to report for each distinct, matched sequence (arbitrarily chosen)
    [default: No limit]

**\--preload-db**

  Cache all DB data in python-land instead of querying for it by SQL
    each time. This is faster particularly for querying many sequences,
    but uses more memory and has a larger start-up time for each marker
    gene.

# OTHER DATABASE EXTRACTION METHODS

**\--sample-names** name [name \...]

  Print all OTUs from these samples

**\--sample-list** path

  Print all OTUs from the samples listed in the file
    (newline-separated)

**\--taxonomy** name

  Print all OTUs assigned a taxonomy including this string e.g.
    \'Archaea\'

**\--dump**

  Print all OTUs in the DB

**\--continue-on-missing-genes**

  Continue if a gene is missing from the DB. Only works with
    smafa/nuclotide search method.

# OTHER GENERAL OPTIONS

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

# AUTHORS

>     Ben J. Woodcroft, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
>     Samuel Aroney, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
>     Raphael Eisenhofer, Centre for Evolutionary Hologenomics, University of Copenhagen, Denmark
>     Rossen Zhao, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
