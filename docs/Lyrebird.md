---
title: Lyrebird (phage profiling)
---
# Lyrebird (phage profiling)

![Lyrebird](./lyrebird_resized.png)

Lyrebird is a module of SingleM specialized towards the profiling of viruses in metagenomic data. Currently, Lyrebird profiles dsDNA phages belonging to the [Caudoviricetes](https://en.wikipedia.org/wiki/Caudoviricetes).

Lyrebird is similar to standard microbial SingleM, with some main conceptual differences: 

1. Since phage lack universal single copy marker genes, Lyrebird uses an expanded set of >500 marker genes, and does not assume each of them is present in each phage. 
2. Lyrebird uses a more sensitive homology search method to recruit reads to marker gene windows. 
3. Lyrebird does not use GTDB taxonomy, but instead uses a custom taxonomy based on [vConTACT3](https://bitbucket.org/MAVERICLab/vcontact3/src/master/) analysis of a custom collection of reference sequences derived from [ICTV](https://ictv.global/) and [IMG-VR](https://doi.org/10.1093/nar/gkac1037).

In testing, Lyrebird shows good ability to detect both known and novel phage, and to estimate their relative abundances. In metagenomic datasets, many more phage sequences are detected compared to standard contig-centric methods. A manuscript describing Lyrebird and its applications is in preparation.

The following subcommands are available:

* [lyrebird pipe](/tools/lyrebird_pipe) - the main workflow which generates OTU tables and [vConTACT3-derived](https://bitbucket.org/MAVERICLab/vcontact3/src/master/) taxonomic profiles.
* [lyrebird renew](/advanced/lyrebird_renew) - Given previously generated results, re-run the pipeline with a new reference sequence/taxonomy database.
* [lyrebird condense](/advanced/lyrebird_condense) - Given an OTU table, summarise the results into a taxonomic profile.

# Lyrebird database versions

Unlike standard SingleM, Lyrebird does not currently maintain backwards
compatibility in its database. It uses a standard SingleM metapackage format,
but the marker genes, window positions and taxonomy strings are specific to the
Lyrebird database version, given the rapid changes in available phage genomes,
marker gene databases and taxonomy. This means that the `renew` subcommand
cannot be used to update the taxonomy of an archive OTU table generated with a previous version of the Lyrebird database.

## v0.3.1
In SingleM/Lyrebird v0.20.0, the default Lyrebird DB is v0.3.1.
* It improves on v0.3.0 by better removal of off-target sequences i.e. microbial homologues of the marker genes that were previously interpreted as being phage.
* It uses 762 marker genes, derived from PHROG v4.1 applied to phage genomes included in [ICTV Master Species List (MSL)](https://ictv.global/msl) 38 or [IMG-VR](https://doi.org/10.1093/nar/gkac1037) v4.1.
* It was built from 69,281 95% ANI-dereplicated phage genomes with CheckV completeness >80%.
* 69,405 (99%) of the dereplicated set had at least one marker gene.
* Including genomes without any marker genes, each genome had 13.6 marker genes on average.
* Taxonomy was assigned using [vConTACT3](https://bitbucket.org/MAVERICLab/vcontact3/src/master/) v3.0.0.b65.

## v0.3.0
In SingleM/Lyrebird v0.19.0, the default Lyrebird DB was v0.3.0.
* It uses 642 marker genes, derived from PHROG v4.1 applied to phage genomes included in [ICTV Master Species List (MSL)](https://ictv.global/msl) 38 or [IMG-VR](https://doi.org/10.1093/nar/gkac1037) v4.1.
* It was built from 69,281 95% ANI-dereplicated phage genomes with CheckV completeness >80%.
* 65,384 (94%) of the dereplicated set had at least one marker gene.
* Including genomes without any marker genes, each genome had 5.3 marker genes on average.
* Taxonomy was assigned using [vConTACT3](https://bitbucket.org/MAVERICLab/vcontact3/src/master/).

# FAQ
## How does Lyrebird handle prophages?
Lyrebird does not attempt to differentiate between free phage and prophage in the community profiles it generates. Instead,  community profiles contain all phage sequences, regardless of whether they are free or integrated. This is due to the technical limitations of read-centric analysis, where an individual read may be derived from either a free phage or a prophage - determining which is challenging or potentially impossible, at least for short reads.
