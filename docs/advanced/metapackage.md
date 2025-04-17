---
title: SingleM metapackage
---
# singlem metapackage

# DESCRIPTION

Create or describe a metapackage (i.e. set of SingleM packages)

# OPTIONS

**\--metapackage** *METAPACKAGE*

  Path to write generated metapackage to

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

**\--taxonomy-database-name** *TAXONOMY_DATABASE_NAME*

  Name of the taxonomy database to use [default:
    custom_taxonomy_database]

**\--taxonomy-database-version** *TAXONOMY_DATABASE_VERSION*

  Version of the taxonomy database to use [default: unspecified]

**\--diamond-prefilter-performance-parameters** *DIAMOND_PREFILTER_PERFORMANCE_PARAMETERS*

  Performance-type arguments to use when calling \'diamond blastx\'
    during the prefiltering. [default: \'\--block-size 0.5
    \--target-indexed -c1\']

**\--diamond-taxonomy-assignment-performance-parameters** *DIAMOND_TAXONOMY_ASSIGNMENT_PERFORMANCE_PARAMETERS*

  Performance-type arguments to use when calling \'diamond blastx\'
    during the taxonomy assignment. [default: \'\--block-size 0.5
    \--target-indexed -c1\']

**\--describe**

  Describe a metapackage rather than create it

**\--threads** num_threads

  number of CPUS to use [default: 1]

**\--prefilter-clustering-threshold** fraction

  ID for dereplication of prefilter DB [default: 0.6]

**\--prefilter-diamond-db** DMND

  Dereplicated DIAMOND db for prefilter to use [default: dereplicate
    from input SingleM packages]

**\--makeidx-sensitivity-params** PARAMS

  DIAMOND sensitivity parameters to use when indexing the prefilter
    DIAMOND db. [default: None]

**\--calculate-average-num-genes-per-species**

  Calculate the average number of genes per species in the
    metapackage. [default: False]

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
