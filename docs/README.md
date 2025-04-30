[![Anaconda-Server Badge](https://anaconda.org/bioconda/singlem/badges/version.svg)](https://anaconda.org/bioconda/singlem)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/singlem/badges/downloads.svg)](https://anaconda.org/bioconda/singlem)

# SingleM
Welcome.

SingleM is a tool for profiling shotgun metagenomes, determining the relative abundance of microbial taxa in a sample. As of version 0.19.0, it can also be used to profile dsDNA phages (see ![Lyrebird logo](./lyrebird_small.png) below).

It shows good accuracy in estimating the relative abundances of community members, and has a particular strength in dealing with novel lineages. The method it uses also makes it suitable for some related tasks, such as assessing eukaryotic contamination, finding bias in genome recovery, and lineage-targeted MAG recovery. It can also be used as the basis for choosing metagenomes which, when coassembled, maximise the recovery of novel MAGs (see the [![BinChicken logo](./binchicken_small.png) documentation](https://aroneys.github.io/binchicken/)).

Microbial SingleM has been applied to ~700,000 public metagenomes. The resulting data are available at the [Sandpiper companion website ![Sandpiper logo](./sandpiper_small.png)](https://sandpiper.qut.edu.au).

The main idea of SingleM is to profile metagenomes by targeting short 20 amino acid stretches ('windows') within single copy marker genes. It finds reads which cover an entire window, and analyses these further. By constraining analysis to these short windows, it becomes possible to know how novel each read is compared to known genomes. Then, using the fact that each analysed gene is (almost always) found exactly once in each genome, the abundance of each lineage can be accurately estimated.

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
![Lyrebird](./lyrebird_resized.png)

Lyrebird is a module of SingleM specialized towards the profiling of viruses in metagenomic data. Currently, Lyrebird profiles dsDNA phages belonging to the [Caudoviricetes](https://en.wikipedia.org/wiki/Caudoviricetes).

Lyrebird is similar to standard microbial SingleM, with two main conceptual differences: 

1. Since phage lack universal single copy marker genes, Lyrebird uses an expanded set of 630 marker genes, and does not assume each of them is present in each phage. 
2. Lyrebird uses a more sensitive homology search method to recruit reads to marker gene windows. 

In testing, Lyrebird shows good ability to detect both known and novel phage, and to estimate their relative abundances. In metagenomic datasets, many more phage sequences are detected compared to standard contig-centric methods. A manuscript describing Lyrebird and its applications is in preparation.

The following subcommands are available:

* [lyrebird pipe](/tools/lyrebird_pipe) - the main workflow which generates OTU tables and [vConTACT3-derived](https://bitbucket.org/MAVERICLab/vcontact3/src/master/) taxonomic profiles.
* [lyrebird renew](/tools/lyrebird_renew) - Given previously generated results, re-run the pipeline with a new reference sequence/taxonomy database.
* [lyrebird condense](/advanced/lyrebird_condense) - Given an OTU table, summarise the results into a taxonomic profile.

## Help
If you have any questions or comments, raise a [GitHib issue](https://github.com/wwood/singlem/issues) or just send us an [email](https://research.qut.edu.au/cmr/team/ben-woodcroft/).

## License
SingleM is developed by the [Woodcroft lab](https://research.qut.edu.au/cmr/team/ben-woodcroft/) at the [Centre for Microbiome Research](https://research.qut.edu.au/cmr), School of Biomedical Sciences, QUT, with contributions several including [Samuel Aroney](https://github.com/AroneyS), [Rossen Zhao](https://github.com/rzhao-2), [Raphael Eisenhofer](https://github.com/EisenRa) and many others. It is licensed under [GPL3 or later](https://gnu.org/licenses/gpl.html).

The source code is available at [https://github.com/wwood/singlem](https://github.com/wwood/singlem).

## Citations
### Profiling microbial communities with SingleM / Sandpiper
Ben J. Woodcroft, Samuel T. N. Aroney, Rossen Zhao, Mitchell Cunningham, Joshua A. M. Mitchell, Linda Blackall, Gene W Tyson. *SingleM and Sandpiper: Robust microbial taxonomic profiles from metagenomic data.* bioRxiv 2024.01.30.578060; doi: [https://doi.org/10.1101/2024.01.30.578060](https://doi.org/10.1101/2024.01.30.578060).
### SingleM microbial_fraction
Raphael Eisenhofer, Antton Alberdi, Ben J. Woodcroft, 2024. *Large-scale estimation of bacterial and archaeal DNA prevalence in metagenomes reveals biome-specific patterns.* bioRxiv, pp.2024-05; [https://doi.org/10.1101/2024.05.16.594470](https://doi.org/10.1101/2024.05.16.594470).
### SingleM-powered coassembly with BinChicken
Samuel T. N. Aroney, Rhys J. Newell, Gene W. Tyson and Ben J. Woodcroft, 2024. Bin Chicken: targeted metagenomic coassembly for the efficient recovery of novel genomes. bioRxiv, pp.2024-11. [https://doi.org/10.1101/2024.11.24.625082](https://doi.org/10.1101/2024.11.24.625082).
### Lyrebird
Rossen Zhao, Gene W. Tyson, Ben J. Woodcroft. *Lyrebird: a tool for profiling dsDNA phage communities in metagenomic data.* (in preparation).
