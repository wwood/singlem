---
title: SingleM supplement
---
# singlem supplement

DESCRIPTION
===========

Create a new metapackage from a vanilla one plus new genomes

OPTIONS
=======

**\--new-genome-fasta-files** *NEW_GENOME_FASTA_FILES* [*NEW_GENOME_FASTA_FILES* \...]

  FASTA files of new genomes

**\--new-genome-fasta-files-list** *NEW_GENOME_FASTA_FILES_LIST* [*NEW_GENOME_FASTA_FILES_LIST* \...]

  File containing FASTA file paths of new genomes

**\--new-taxonomies** *NEW_TAXONOMIES*

  newline separated file containing taxonomies of new genomes
    (path\<TAB\>taxonomy). Must be fully specified to species level. If
    not specified, the taxonomy will be inferred from the new genomes
    using GTDB-tk

**\--input-metapackage** *INPUT_METAPACKAGE*

  metapackage to build upon [default: Use default package]

**\--output-metapackage** *OUTPUT_METAPACKAGE*

  output metapackage

**\--threads** *THREADS*

  parallelisation

**\--pplacer-threads** *PPLACER_THREADS*

  for GTDBtk classify_wf

**\--working-directory** *WORKING_DIRECTORY*

  working directory [default: use a temporary directory]

**\--gtdbtk-output-directory** *GTDBTK_OUTPUT_DIRECTORY*

  use this GTDBtk result. Not used if \--new-taxonomies is used
    [default: not set, run GTDBtk]

**\--taxonomy-file** *TAXONOMY_FILE*

  A 2 column tab-separated file containing each genome\'s taxonomy as
    output by GTDBtk [default: not set, run GTDBtk]

**\--output-taxonomies** *OUTPUT_TAXONOMIES*

  TSV output file of taxonomies of new genomes, whether they are novel
    species or not.

**\--checkm2-quality-file** *CHECKM2_QUALITY_FILE*

  CheckM2 quality file of new genomes

**\--no-quality-filter**

  skip quality filtering

**\--no-taxon-genome-lengths**

  Do not include taxon genome lengths in updated metapackage

**\--no-dereplication**

  Assume genome inputs are already dereplicated

**\--dereplicate-with-galah**

  Run galah to dereplicate genomes at species level

**\--checkm2-min-completeness** *CHECKM2_MIN_COMPLETENESS*

  minimum completeness for CheckM2 [default: 70]

**\--checkm2-max-contamination** *CHECKM2_MAX_CONTAMINATION*

  maximum contamination for CheckM2 [default: 10]

**\--hmmsearch-evalue** *HMMSEARCH_EVALUE*

  evalue for hmmsearch run on proteins to gather markers [default:
    1e-20]

**\--skip-taxonomy-check**

  skip check which ensures that GTDBtk assigned taxonomies are
    concordant with the old metapackage\'s

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
>     Raphael Eisenhofer, Centre for Evolutionary Hologenomics, University of Copenhagen, Denmark
>     Rossen Zhao, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
