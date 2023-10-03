---
title: SingleM metapackage
---
# singlem metapackage

DESCRIPTION
===========

Create or describe a metapackage (i.e. set of SingleM packages)

OPTIONS
=======

**\--singlem-packages** *SINGLEM_PACKAGES* [*SINGLEM_PACKAGES* \...]

  Input packages

**\--nucleotide-sdb** *NUCLEOTIDE_SDB*

  Nucleotide SingleM database for initial assignment pass

**\--no-nucleotide-sdb**

  Skip nucleotide SingleM database

**\--taxon-genome-lengths** *TAXON_GENOME_LENGTHS*

  TSV file of genome lengths for each taxon

**\--no-taxon-genome-lengths**

  Skip taxon genome lengths

**\--metapackage** *METAPACKAGE*

  Path to write generated metapackage to

**\--describe**

  Describe a metapackage rather than create it

**\--threads** num_threads

  number of CPUS to use [default: 1]

**\--prefilter-clustering-threshold** fraction

  ID for dereplication of prefilter DB [default: 0.6]

**\--prefilter-diamond-db** DMND

  Dereplicated DIAMOND db for prefilter to use [default: dereplicate
    from input SingleM packages]

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
