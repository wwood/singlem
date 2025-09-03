---
title: SingleM supplement
---
# singlem supplement

**TLDR**: Supplement a metapackage with new genomes like so:
```
singlem supplement --new-genome-fasta-files <genome1.fna> <genome2.fna> \
    --output-metapackage <supplemented.smpkg>
```
GTDB-Tk may be needed to assign taxonomy to the new genomes. If required, install it separately (e.g. `conda install gtdbtk`) and ensure its version matches the GTDB release of the metapackage.

# DESCRIPTION

Create a new metapackage from a vanilla one plus new genomes

# OPTIONS

**\--new-genome-fasta-files** *NEW_GENOME_FASTA_FILES* [*NEW_GENOME_FASTA_FILES* \...]

  FASTA files of new genomes

**\--new-genome-fasta-files-list** *NEW_GENOME_FASTA_FILES_LIST* [*NEW_GENOME_FASTA_FILES_LIST* \...]

  File containing FASTA file paths of new genomes

**\--input-metapackage** *INPUT_METAPACKAGE*

  metapackage to build upon [default: Use default package]

**\--output-metapackage** *OUTPUT_METAPACKAGE*

  output metapackage

**\--threads** *THREADS*

  parallelisation

# TAXONOMY

**\--new-fully-defined-taxonomies** *NEW_FULLY_DEFINED_TAXONOMIES*

  newline separated file containing taxonomies of new genomes
    (path\<TAB\>taxonomy). Must be fully specified to species level. If
    not specified, the taxonomy will be inferred from the new genomes
    using GTDB-tk or read from \--taxonomy-file [default: not set, run
    GTDBtk].

**\--taxonomy-file** *TAXONOMY_FILE*

  A 2 column tab-separated file containing each genome\'s taxonomy as
    output by GTDBtk [default: not set, run GTDBtk]

**\--gtdbtk-output-directory** *GTDBTK_OUTPUT_DIRECTORY*

  use this GTDBtk result. Not used if \--new-taxonomies is used
    [default: not set, run GTDBtk]

**\--pplacer-threads** *PPLACER_THREADS*

  for GTDBtk classify_wf

**\--output-taxonomies** *OUTPUT_TAXONOMIES*

  TSV output file of taxonomies of new genomes, whether they are novel
    species or not.

**\--skip-taxonomy-check**

  skip check which ensures that GTDBtk assigned taxonomies are
    concordant with the old metapackage\'s [default: do the check]

# QUALITY FILTERING OF NEW GENOMES

**\--checkm2-quality-file** *CHECKM2_QUALITY_FILE*

  CheckM2 quality file of new genomes

**\--no-quality-filter**

  skip quality filtering

**\--checkm2-min-completeness** *CHECKM2_MIN_COMPLETENESS*

  minimum completeness for CheckM2 [default: 50]

**\--checkm2-max-contamination** *CHECKM2_MAX_CONTAMINATION*

  maximum contamination for CheckM2 [default: 10]

# DEREPLICATION

**\--no-dereplication**

  Assume genome inputs are already dereplicated

**\--dereplicate-with-galah**

  Run galah to dereplicate genomes at species level

# LESS COMMON OPTIONS

**\--hmmsearch-evalue** *HMMSEARCH_EVALUE*

  evalue for hmmsearch run on proteins to gather markers [default:
    1e-20]

**\--gene-definitions** *GENE_DEFINITIONS*

  Tab-separated file of
    genome_fasta\<TAB\>transcript_fasta\<TAB\>protein_fasta [default:
    undefined, call genes using Prodigal]

**\--working-directory** *WORKING_DIRECTORY*

  working directory [default: use a temporary directory]

**\--no-taxon-genome-lengths**

  Do not include taxon genome lengths in updated metapackage

**\--ignore-taxonomy-database-incompatibility**

  Do not halt when the old metapackage is not the default metapackage.

**\--new-taxonomy-database-name** *NEW_TAXONOMY_DATABASE_NAME*

  Name of the taxonomy database to record in the created metapackage
    [default: custom_taxonomy_database]

**\--new-taxonomy-database-version** *NEW_TAXONOMY_DATABASE_VERSION*

  Version of the taxonomy database to use [default: None]

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
