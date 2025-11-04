---
title: SingleM regenerate
---
# singlem regenerate

# DESCRIPTION

Update a SingleM package with new sequences and taxonomy (expert mode).

# OPTIONS

**\--min-aligned-percent** percent

  remove sequences from the alignment which do not cover this
    percentage of the HMM [default: 10]

**\--window-position** *WINDOW_POSITION*

  change window position of output package [default: do not change]

**\--sequence-prefix** *SEQUENCE_PREFIX*

  add a prefix to sequence names

**\--candidate-decoy-sequences**, **\--euk-sequences** *CANDIDATE_DECOY_SEQUENCES*

  candidate amino acid sequences fasta file to search for decoys

**\--candidate-decoy-taxonomy**, **\--euk-taxonomy** *CANDIDATE_DECOY_TAXONOMY*

  tab-separated sequence ID to taxonomy file of candidate decoy
    sequences

**\--no-candidate-decoy-sequences**, **\--no-further-euks**

  Do not include any euk sequences beyond what is already in the
    current SingleM package

# REQUIRED ARGUMENTS

**\--input-singlem-package** PATH

  input package path

**\--output-singlem-package** PATH

  output package path

**\--input-sequences** *INPUT_SEQUENCES*

  all on-target amino acid sequences fasta file for new package

**\--input-taxonomy** *INPUT_TAXONOMY*

  tab-separated sequence ID to taxonomy file of on-target sequence
    taxonomy

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
