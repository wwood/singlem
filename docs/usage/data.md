---
title: "singlem data"
author: "Ben Woodcroft, Centre for Microbiome Research, Queensland University of Technology"
date: "`r Sys.Date()` (`r system('bin/singlem --version',intern=T)`)"
---
NAME
====

singlem data

SYNOPSIS
========

**singlem** data [-h] [\--output-directory OUTPUT_DIRECTORY]
[\--verify-only] [\--debug] [\--version] [\--quiet]
[\--full-help] [\--full-help-roff]

DESCRIPTION
===========

Download reference metapackage data

OPTIONS
=======

**\--output-directory** *OUTPUT_DIRECTORY*

:   Output directory [required unless SINGLEM_METAPACKAGE_PATH is
    specified]

**\--verify-only**

:   Check that the data is up to date and each file has the correct
    checksum

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
