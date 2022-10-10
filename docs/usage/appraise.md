---
title: "singlem appraise"
author: "Ben Woodcroft, Centre for Microbiome Research, Queensland University of Technology"
date: "`r Sys.Date()` (`r system('bin/singlem --version',intern=T)`)"
---
NAME
====

singlem appraise

SYNOPSIS
========

**singlem** appraise [-h] [\--metagenome-otu-tables
METAGENOME_OTU_TABLES [METAGENOME_OTU_TABLES \...]]
[\--metagenome-archive-otu-tables METAGENOME_ARCHIVE_OTU_TABLES
[METAGENOME_ARCHIVE_OTU_TABLES \...]] [\--genome-otu-tables
GENOME_OTU_TABLES [GENOME_OTU_TABLES \...]]
[\--genome-archive-otu-tables GENOME_ARCHIVE_OTU_TABLES
[GENOME_ARCHIVE_OTU_TABLES \...]] [\--assembly-otu-tables
ASSEMBLY_OTU_TABLES [ASSEMBLY_OTU_TABLES \...]]
[\--assembly-archive-otu-tables ASSEMBLY_ARCHIVE_OTU_TABLES
[ASSEMBLY_ARCHIVE_OTU_TABLES \...]] \--metapackage METAPACKAGE
[\--imperfect] [\--sequence-identity SEQUENCE_IDENTITY] [\--plot
PLOT] [\--plot-marker PLOT_MARKER] [\--plot-basename PLOT_BASENAME]
[\--output-binned-otu-table OUTPUT_BINNED_OTU_TABLE]
[\--output-unbinned-otu-table OUTPUT_UNBINNED_OTU_TABLE]
[\--output-assembled-otu-table OUTPUT_ASSEMBLED_OTU_TABLE]
[\--output-unaccounted-for-otu-table OUTPUT_UNACCOUNTED_FOR_OTU_TABLE]
[\--output-found-in] [\--debug] [\--version] [\--quiet]
[\--full-help] [\--full-help-roff]

DESCRIPTION
===========

How much of the metagenome do the genomes or assembly represent?

INPUT OTU TABLE OPTIONS
=======================

**\--metagenome-otu-tables** *METAGENOME_OTU_TABLES* [*METAGENOME_OTU_TABLES* \...]

:   output of \'pipe\' run on metagenomes

**\--metagenome-archive-otu-tables** *METAGENOME_ARCHIVE_OTU_TABLES* [*METAGENOME_ARCHIVE_OTU_TABLES* \...]

:   archive output of \'pipe\' run on metagenomes

**\--genome-otu-tables** *GENOME_OTU_TABLES* [*GENOME_OTU_TABLES* \...]

:   output of \'pipe\' run on genomes

**\--genome-archive-otu-tables** *GENOME_ARCHIVE_OTU_TABLES* [*GENOME_ARCHIVE_OTU_TABLES* \...]

:   archive output of \'pipe\' run on genomes

**\--assembly-otu-tables** *ASSEMBLY_OTU_TABLES* [*ASSEMBLY_OTU_TABLES* \...]

:   output of \'pipe\' run on assembled sequence

**\--assembly-archive-otu-tables** *ASSEMBLY_ARCHIVE_OTU_TABLES* [*ASSEMBLY_ARCHIVE_OTU_TABLES* \...]

:   archive output of \'pipe\' run on assembled sequence

**\--metapackage** *METAPACKAGE*

:   Metapackage used in the creation of the OTU tables

INEXACT APPRAISAL OPTIONS
=========================

**\--imperfect**

:   use sequence searching to account for genomes that are similar to
    those found in the metagenome [default: False]

**\--sequence-identity** *SEQUENCE_IDENTITY*

:   sequence identity cutoff to use if \--imperfect is specified
    [default: \~genus level divergence i.e. 0.86]

PLOTTING-RELATED OPTIONS
========================

**\--plot** *PLOT*

:   Output plot SVG filename (marker chosen automatically unless
    \--plot-marker is also specified)

**\--plot-marker** *PLOT_MARKER*

:   Marker gene to plot OTUs from

**\--plot-basename** *PLOT_BASENAME*

:   Plot visualisation of appraisal results from all markers to this
    basename (one SVG per marker)

OUTPUT SUMMARY OTU TABLES
=========================

**\--output-binned-otu-table** *OUTPUT_BINNED_OTU_TABLE*

:   output OTU table of binned populations

**\--output-unbinned-otu-table** *OUTPUT_UNBINNED_OTU_TABLE*

:   output OTU table of assembled but not binned populations

**\--output-assembled-otu-table** *OUTPUT_ASSEMBLED_OTU_TABLE*

:   output OTU table of all assembled populations

**\--output-unaccounted-for-otu-table** *OUTPUT_UNACCOUNTED_FOR_OTU_TABLE*

:   Output OTU table of populations not accounted for

**\--output-found-in**

:   Output sample name (genome or assembly) the hit was found in

OTHER GENERAL OPTIONS
=====================

**\--debug**

:   output debug information

**\--version**

:   output version information and quit

**\--quiet**

:   only output errors

**\--full-help**

:   print longer help message

**\--full-help-roff**

:   print longer help message in ROFF (manpage) format

AUTHORS
=======

>     Ben J. Woodcroft, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
>     Samuel Aroney, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
>     Rossen Zhao, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
