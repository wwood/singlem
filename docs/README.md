[![Anaconda-Server Badge](https://anaconda.org/bioconda/singlem/badges/version.svg)](https://anaconda.org/bioconda/singlem)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/singlem/badges/downloads.svg)](https://anaconda.org/bioconda/singlem)

<!-- NOTE: This intro should manually be kept in sync between the repo README and the docs README -->

[![Sandpiper logo](./sandpiper_small.png)](https://sandpiper.qut.edu.au)
[![Bin Chicken logo](./binchicken_small.png)](https://aroneys.github.io/binchicken)
[![Lyrebird logo](./lyrebird_small.png)](/Lyrebird)

Welcome.

At heart, SingleM is a tool for profiling shotgun metagenomes. It was originally designed to determine the relative abundance of bacterial and archaeal taxa in a sample. As of version 0.19.0, it can also be used to profile dsDNA phages (see [Lyrebird](/Lyrebird)).

It [shows](https://doi.org/10.1101/2024.01.30.578060) good accuracy in estimating the relative abundances of community members, and has a particular strength in dealing with novel lineages. The method it uses also makes it suitable for some related tasks, such as assessing eukaryotic contamination, finding bias in genome recovery, and lineage-targeted MAG recovery. It can also be used as the basis for choosing metagenomes which, when coassembled, maximise the recovery of novel MAGs (see [Bin Chicken](https://aroneys.github.io/binchicken/)).

Microbial SingleM has been applied to ~700,000 public metagenomes. The resulting data are available at the [Sandpiper companion website](https://sandpiper.qut.edu.au).

The main idea of SingleM is to profile metagenomes by targeting short 20 amino acid stretches ("*windows*") within single copy marker genes. It finds reads which cover an entire window, and analyses these further. By constraining analysis to these short windows, it becomes possible to know how novel each read is compared to known genomes. Then, using the fact that each analysed gene is (almost always) found exactly once in each genome, the abundance of each lineage can be accurately estimated.

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

## Help
If you have any questions or comments, raise a [GitHib issue](https://github.com/wwood/singlem/issues) or just send us an [email](https://research.qut.edu.au/cmr/team/ben-woodcroft/).

## License
SingleM is developed by the [Woodcroft lab](https://research.qut.edu.au/cmr/team/ben-woodcroft/) at the [Centre for Microbiome Research](https://research.qut.edu.au/cmr), School of Biomedical Sciences, QUT, with contributions from many helpful people including [Samuel Aroney](https://github.com/AroneyS), [Rossen Zhao](https://github.com/rzhao-2), [Raphael Eisenhofer](https://github.com/EisenRa). It is licensed under [GPL3 or later](https://gnu.org/licenses/gpl.html).

The source code is available at [https://github.com/wwood/singlem](https://github.com/wwood/singlem).

## Citations
<!-- NOTE: Citations should manually be kept in sync between the repo README and the docs README -->
### Profiling microbial communities with SingleM / Sandpiper
Woodcroft B.J., Aroney T.N.S., Zhao R., Cunningham M., Mitchell J.A.M., Nurdiansyah R., Blackall L., Tyson G.W.. *Comprehensive taxonomic identification of microbial species in metagenomic data using SingleM and Sandpiper.* Nat Biotechnol (2025). https://doi.org/10.1038/s41587-025-02738-1
### SingleM microbial_fraction
Raphael Eisenhofer, Antton Alberdi, Ben J. Woodcroft, 2024. *Large-scale estimation of bacterial and archaeal DNA prevalence in metagenomes reveals biome-specific patterns.* bioRxiv, pp.2024-05; [https://doi.org/10.1101/2024.05.16.594470](https://doi.org/10.1101/2024.05.16.594470).
### SingleM-powered coassembly with Bin Chicken
Samuel T. N. Aroney, Rhys J. Newell, Gene W. Tyson and Ben J. Woodcroft, 2024. *Bin Chicken: targeted metagenomic coassembly for the efficient recovery of novel genomes.* bioRxiv, pp.2024-11. [https://doi.org/10.1101/2024.11.24.625082](https://doi.org/10.1101/2024.11.24.625082).
### Lyrebird
Rossen Zhao, Gene W. Tyson, Ben J. Woodcroft. *Lyrebird: a tool for profiling dsDNA phage communities in metagenomic data.* (in preparation).
