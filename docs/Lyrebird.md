---
title: Lyrebird (phage profiling)
---
# Lyrebird (phage profiling)

![Lyrebird](./lyrebird_resized.png)

Lyrebird is a module of SingleM specialized towards the profiling of viruses in metagenomic data. Currently, Lyrebird profiles dsDNA phages belonging to the [Caudoviricetes](https://en.wikipedia.org/wiki/Caudoviricetes).

Lyrebird is similar to standard microbial SingleM, with some main conceptual differences: 

1. Since phage lack universal single copy marker genes, Lyrebird uses an expanded set of 630 marker genes, and does not assume each of them is present in each phage. 
2. Lyrebird uses a more sensitive homology search method to recruit reads to marker gene windows. 
3. Lyrebird does not use GTDB taxonomy, but instead uses a custom taxonomy based on [vConTACT3](https://bitbucket.org/MAVERICLab/vcontact3/src/master/) analysis of a custom collection of reference sequences derived from [ICTV](https://ictv.global/) and [IMG-VR](https://doi.org/10.1093/nar/gkac1037).

In testing, Lyrebird shows good ability to detect both known and novel phage, and to estimate their relative abundances. In metagenomic datasets, many more phage sequences are detected compared to standard contig-centric methods. A manuscript describing Lyrebird and its applications is in preparation.

The following subcommands are available:

* [lyrebird pipe](/tools/lyrebird_pipe) - the main workflow which generates OTU tables and [vConTACT3-derived](https://bitbucket.org/MAVERICLab/vcontact3/src/master/) taxonomic profiles.
* [lyrebird renew](/tools/lyrebird_renew) - Given previously generated results, re-run the pipeline with a new reference sequence/taxonomy database.
* [lyrebird condense](/advanced/lyrebird_condense) - Given an OTU table, summarise the results into a taxonomic profile.

## FAQ
### How does Lyrebird handle prophages?
Lyrebird does not attempt to differentiate between free phage and prophage in the community profiles it generates. Instead,  community profiles contain all phage sequences, regardless of whether they are free or integrated. This is due to the technical limitations of read-centric analysis, where an individual read may be derived from either a free phage or a prophage - determining which is challenging or potentially impossible, at least for short reads.
