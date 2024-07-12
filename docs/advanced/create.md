---
title: SingleM create
---
# singlem create

# DESCRIPTION

Create a SingleM package.

# OPTIONS

**\--input-graftm-package** PATH

  Input GraftM package underlying the new SingleM package. The GraftM
    package is usually made with \'graftM create \--no_tree \--hmm
    \<your.hmm\>\' where \<your.hmm\> is the one provided to \'singlem
    seqs\'.

**\--input-taxonomy** PATH

  Input taxonomy file in GreenGenes format (2 column tab separated, ID
    then taxonomy with taxonomy separated by \';\' or \'; \'.

**\--output-singlem-package** PATH

  Output package path

**\--hmm-position** INTEGER

  Position in the GraftM alignment HMM where the SingleM window
    starts. To choose the best position, use \'singlem seqs\'. Note that
    this position (both the one output by \'seqs\' and the one specified
    here) is a 1-based index, but this positions stored within the
    SingleM package as a 0-based index.

**\--window-size** INTEGER

  Length of NUCLEOTIDE residues in the window, counting only those
    that match the HMM [default: 60]

**\--target-domains** *TARGET_DOMAINS* [*TARGET_DOMAINS* \...]

  Input domains targeted by this package e.g. \'Archaea\',
    \'Bacteria\', \'Eukaryota\' or \'Viruses\'. Input with multiple
    domains must be space separated.

**\--gene-description** STRING

  Input free form text description of this marker package, for use
    with \'singlem metapackage \--describe\'.

**\--force**

  Overwrite output path if it already exists [default: false]

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
