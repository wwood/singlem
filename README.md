

<!-- NOTE: This intro should manually be kept in sync between the repo README and the docs README -->

[![Sandpiper logo](./docs/_include/sandpiper_small.png)](https://sandpiper.qut.edu.au)
[![Bin Chicken logo](./docs/_include/binchicken_small.png)](https://aroneys.github.io/binchicken)
[![Lyrebird logo](./docs/_include/lyrebird_small.png)](/Lyrebird)


[![Current Build](https://github.com/wwood/singlem/actions/workflows/test-singlem.yml/badge.svg)](https://github.com/wwood/singlem/actions)
[![Conda version](https://img.shields.io/conda/vn/bioconda/singlem?label=Conda&color=43b02a)](https://anaconda.org/bioconda/singlem)
[![Conda downloads](https://img.shields.io/conda/dn/bioconda/singlem?label=Downloads&color=43b02a)](https://anaconda.org/bioconda/singlem)
[![Docker version](https://img.shields.io/docker/v/wwood/singlem?label=Docker&color=1D63ED)](https://hub.docker.com/r/wwood/singlem/tags)
[![Docker pulls](https://img.shields.io/docker/pulls/wwood/singlem.svg?label=Pulls&color=1D63ED)](https://hub.docker.com/r/wwood/singlem)
[![PyPI version](https://img.shields.io/pypi/v/singlem.svg?label=PyPI&color=ffd43b)](https://pypi.org/project/singlem/)
[![PyPI downloads](https://img.shields.io/pypi/dm/singlem.svg?label=Downloads&color=ffd43b)](https://pypi.org/project/singlem/)

Welcome.

***💉 Now works on dsDNA phages!***  
***🧬 Now works with long reads!***

At heart, SingleM is a tool for profiling shotgun and long-read metagenomes. It was originally designed to determine the relative abundance of bacterial and archaeal taxa in a sample. As of version 0.19.0, it can also be used to profile dsDNA phages (see [Lyrebird](https://wwood.github.io/singlem/Lyrebird)). Long-read compatibility was added in version 0.19.1, but requires `Nanopore >= R10.4.1` or `PacBio HiFi` reads to ensure high enough base-pair calling accuracy for reliable taxonomic profiling.

It [shows](https://doi.org/10.1101/2024.01.30.578060) good accuracy in estimating the relative abundances of community members, and has a particular strength in dealing with novel lineages. The method it uses also makes it suitable for some related tasks, such as assessing eukaryotic contamination, finding bias in genome recovery, and lineage-targeted MAG recovery. It can also be used as the basis for choosing metagenomes which, when coassembled, maximise the recovery of novel MAGs (see [Bin Chicken](https://aroneys.github.io/binchicken/)).

Microbial SingleM has been applied to ~700,000 public metagenomes. The resulting data are available at the [Sandpiper companion website](https://sandpiper.qut.edu.au).

Documentation can be found at https://wwood.github.io/singlem/

## Citations
<!-- NOTE: Citations should manually be kept in sync between the repo README and the docs README -->
### Profiling microbial communities with SingleM / Sandpiper
Ben J. Woodcroft, Samuel T. N. Aroney, Rossen Zhao, Mitchell Cunningham, Joshua A. M. Mitchell, Rizky Nurdiansyah, Linda Blackall & Gene W. Tyson. *Comprehensive taxonomic identification of microbial species in metagenomic data using SingleM and Sandpiper.* Nat Biotechnol (2025). [https://doi.org/10.1038/s41587-025-02738-1](https://doi.org/10.1038/s41587-025-02738-1).
### SingleM microbial_fraction
Raphael Eisenhofer, Antton Alberdi, Ben J. Woodcroft, 2024. *Large-scale estimation of bacterial and archaeal DNA prevalence in metagenomes reveals biome-specific patterns.* bioRxiv, pp.2024-05; [https://doi.org/10.1101/2024.05.16.594470](https://doi.org/10.1101/2024.05.16.594470).
### SingleM-powered coassembly with Bin Chicken
Samuel T. N. Aroney, Rhys J. Newell, Gene W. Tyson and Ben J. Woodcroft, 2024. *Bin Chicken: targeted metagenomic coassembly for the efficient recovery of novel genomes.* bioRxiv, pp.2024-11. [https://doi.org/10.1101/2024.11.24.625082](https://doi.org/10.1101/2024.11.24.625082).
### Lyrebird
Rossen Zhao, Gene W. Tyson, Ben J. Woodcroft. *Lyrebird: a tool for profiling dsDNA phage communities in metagenomic data.* (in preparation).
