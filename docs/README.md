[![Anaconda-Server Badge](https://anaconda.org/bioconda/singlem/badges/version.svg)](https://anaconda.org/bioconda/singlem)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/singlem/badges/downloads.svg)](https://anaconda.org/bioconda/singlem)

# SingleM
Welcome.

SingleM is a tool for profiling shotgun metagenomes. It shows good accuracy in estimating the relative abundances of microbial community members, and has a particular strength in dealing with novel lineages. The method it uses also makes it suitable for some related tasks, such as assessing eukaryotic contamination, finding bias in genome recovery and lineage-targeted MAG recovery.

SingleM has been applied to ~250,000 public metagenomes. The resulting data are available at a companion website [Sandpiper](https://sandpiper.qut.edu.au).

The main idea of SingleM is to profile metagenomes by targeting short 20 amino acid stretches (windows) within single copy marker genes. It finds reads which cover an entire window, and analyses these further. By constraining analysis to these short windows, it becomes possible to know how novel each read is compared to known genomes. Then, using the fact that each analysed gene is (almost always) found exactly once in each genome, the abundance of each lineage can be accurately estimated.

It is currently aimed at the analysis of metagenomes sequenced using Illumina short read technology.

There are several tools (subcommands) which can be used after [installation](/Installation):

* [singlem pipe](/tools/pipe) - the main workflow which generates OTU tables and [GTDB](https://gtdb.ecogenomic.org/) taxonomic profiles. 
* [single summarise](/tools/summarise) - Mechanical transformations of `singlem pipe` results.
* [singlem renew](/tools/renew) - Given previously generated results, re-run the pipeline with a new reference sequence/taxonomy database.
* [singlem supplement](/tools/supplement) - Add new genomes to a reference metapackage.
* [singlem microbial_fraction](/tools/microbial_fraction) - How much of a metagenome is prokaryotic?
* [singlem appraise](/tools/appraise) - How much of a metagenome do the genomes or assembly represent?

And more specialised / expert modes:

* [singlem condense](/advanced/condense) - Given an OTU table, summarise the results into a taxonomic profile.
* [singlem makedb](/advanced/makedb) & [query](/advanced/query)- Create a database of OTU sequences and query it using various sequence similarity methods e.g. [smafa](https://github.com/wwood/smafa).

## Lyrebird
Lyrebird is a module of SingleM specialized towards the profiling of viruses in metagenomic data. Currently, Lyrebird is only able to profile dsDNA phages belonging to the [Caudoviricetes](https://en.wikipedia.org/wiki/Caudoviricetes).

The following subcommands are available for Lyrebird:

* [lyrebird pipe](/tools/lyrebird_pipe) - the main workflow which generates OTU tables and [vConTACT3-derived](https://bitbucket.org/MAVERICLab/vcontact3/src/master/) taxonomic profiles.
* [lyrebird renew](/tools/lyrebird_renew) - Given previously generated results, re-run the pipeline with a new reference sequence/taxonomy database.
* [lyrebird condense](/advanced/lyrebird_condense) - Given an OTU table, summarise the results into a taxonomic profile.

## Help
If you have any questions or comments, raise a [GitHib issue](https://github.com/wwood/singlem/issues) or just send us an [email](https://research.qut.edu.au/cmr/team/ben-woodcroft/).

## License
SingleM is developed by the [Woodcroft lab](https://research.qut.edu.au/cmr/team/ben-woodcroft/) at the [Centre for Microbiome Research](https://research.qut.edu.au/cmr), School of Biomedical Sciences, QUT, with contributions several including [Samuel Aroney](https://github.com/AroneyS) and [Rossen Zhao](https://github.com/rzhao-2) and many others. It is licensed under [GPL3 or later](https://gnu.org/licenses/gpl.html).

The source code is available at [https://github.com/wwood/singlem](https://github.com/wwood/singlem).

## Citation
SingleM and Sandpiper: Robust microbial taxonomic profiles from metagenomic data. Ben J Woodcroft, Samuel T. N. Aroney, Rossen Zhao, Mitchell Cunningham, Joshua A. M. Mitchell, Linda Blackall, Gene W Tyson. bioRxiv 2024.01.30.578060; doi: https://doi.org/10.1101/2024.01.30.578060
