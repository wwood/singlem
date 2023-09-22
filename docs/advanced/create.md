---
title: SingleM create
---
# singlem create

DESCRIPTION
===========

Create a SingleM package.

OPTIONS
=======

**\--input-graftm-package** PATH

  input package

**\--input-taxonomy** PATH

  input taxonomy file in GreenGenes format (2 column tab separated, ID
    then taxonomy with taxonomy separated by \';\' or \'; \'.

**\--output-singlem-package** PATH

  output package

**\--hmm-position** INTEGER

  position in the GraftM alignment HMM where the SingleM window starts

**\--window-size** INTEGER

  length of NUCLEOTIDE residues in the window, counting only those
    that match the HMM [default: 60]

**\--target_domains** *TARGET_DOMAINS* [*TARGET_DOMAINS* \...]

  input domains targeted by this package e.g. \'Archaea\',
    \'Bacteria\', \'Eukaryota\' or \'Viruses\'. Input with multiple
    domains must be space separated.

**\--gene-description** STRING

  input long description of this marker package

**\--force**

  overwrite output path if it already exists [default: false]

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
