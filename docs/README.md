[![Anaconda-Server Badge](https://anaconda.org/bioconda/singlem/badges/version.svg)](https://anaconda.org/bioconda/singlem)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/singlem/badges/downloads.svg)](https://anaconda.org/bioconda/singlem)

# SingleM
Welcome.

SingleM is a tool for profiling shotgun metagenomes. It shows good accuracy in estimating the relative abundances of microbial community members, and has a particular strength in dealing with novel lineages. The method it uses also makes it suitable for some related tasks, such as assessing eukaryotic contamination, finding bias in genome recovery, computing ecological diversity metrics, and lineage-targeted MAG recovery.

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

## Help
If you have any questions or comments, raise a [GitHib issue](https://github.com/wwood/singlem/issues) or just send us an [email](https://research.qut.edu.au/cmr/team/ben-woodcroft/).

### Glossary

* **Taxonomic profile** - A tab-separated table containing the estimated abundances of GTDB taxons in a metagenome. It is in TSV format with 3 columns, with each row corresponding to a taxon. A taxonomic profile may also be called a **condensed profile**, since it is the output of the `condense` algorithm within the main `pipe` workflow. Columns:
  1. sample name. A taxonomic profile can consist of more than one sample. Usually all the taxons in the first sample are listed, and then the taxons in the second sample, and so on.
  2. coverage of that taxon. This is an approximation of the total read coverage of all genomes from this taxon. However, note that this coverage does not include the coverage of sub-taxons. For instance, the coverage of a species is not included in the coverage shown for its genus.
  3. taxonomy string of the taxon
```
sample	coverage	taxonomy
ERR1914274	0	Root
ERR1914274	3.16	Root; d__Bacteria
ERR1914274	0	Root; d__Bacteria; p__Pseudomonadota
ERR1914274	0.06	Root; d__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria
ERR1914274	0	Root; d__Bacteria; p__Bacillota_A
ERR1914274	0.61	Root; d__Bacteria; p__Bacillota_A; c__Clostridia
ERR1914274	0	Root; d__Bacteria; p__Bacteroidota
ERR1914274	0.39	Root; d__Bacteria; p__Bacteroidota; c__Bacteroidia
ERR1914274	0	Root; d__Bacteria; p__Bacillota
...
```

* **OTU table** - A table containing window sequences per metagenome/contig and marker gene. It may be in default form (a TSV with 6 columns, like below), or an extended form with more detail in further columns. The default OTU table output from [pipe](/tools/pipe), [renew](/tools/renew) and [summarise](/tools/summarise) subcommands has 6 columns, with one sequence per row. The extended form OTU table and archive OTU tables have further information (see below). Columns of a default OTU table:
  1. marker name
  2. sample name
  3. sequence of the OTU
  4. number of reads detected from that OTU
  5. estimated coverage of a genome from this OTU
  6. "median" taxonomic classification of each of the reads in the OTU i.e. the most specific taxonomy that 50%+ of the reads agree with.
```
gene    sample  sequence        num_hits        coverage        taxonomy
4.21.ribosomal_protein_S19_rpsS my_sequences  TGGTCGCGCCGTTCGACGGTCACTCCGGACTTCATCGGCCTACAGTTCGCCGTGCACATC    1       1.64    Root; d__Bacteria; p__Proteobacteria; c__Deltaproteobacteria; o__Desulfuromonadales
4.21.ribosomal_protein_S19_rpsS my_sequences  TGGTCGCGGCGCTCAACCATTCTGCCCGAGTTCGTCGGCCACACCGTGGCCGTTCACAAC    1       1.64    Root; d__Bacteria; p__Acidobacteria; c__Solibacteres; o__Solibacterales; f__Solibacteraceae; g__Candidatus_Solibacter; s__Candidatus_Solibacter_usitatus
```

* **OTU table (extended form)** The extended OTU table form generated with the `--output-extras` option to the [pipe](/tools/pipe), [renew](/tools/renew) and [summarise](/tools/summarise) subcommands, has all the columns of a regular OTU table, but with several additional columns which contain more information about each OTU:
  1. read_names - the names of the reads which encode the OTU sequence
  2. nucleotides_aligned - the number of nucleotides which aligned to the window (usually 60, but can be more or less if there are gaps or inserts)
  3. taxonomy_by_known? - whether the taxonomy of the OTU was determined by known genomes (TRUE) or by the reads themselves (FALSE). Currently this is a disused column and is always marked FALSE.
  4. read_unaligned_sequences - the raw sequences of the reads which encode the OTU sequence
  5. equal_best_hit_taxonomies - the taxonomies of the best hits to the OTU sequence, if there are multiple equally good hits. This is a JSON array of strings.


* **Archive OTU table** - Similar to an extended form OTU table, but in JSON form for machine readability and with formatting version recorded. The [renew](/tools/renew) subcommand which re-analyses a dataset requires this format of OTU table rather than the default tab-separated OTU table format. The canonical file extension for SingleM packages is `.json`.
* **SingleM package (spkg)** - Reference data for one particular marker gene and its window position. The canonical file extension for SingleM packages is `.spkg`.
* **SingleM metapackage** - A collection of SingleM packages, with additional indices. The canonical file extension for SingleM metapackages is `.smpkg`.
* **SingleM database** - An OTU table which has been converted to SQLite3 format and sequence similarity search indexes. Canonically SingleM databases are named with the `.sdb` extension, but this is not enforced. SingleM databases are created with the [makedb](/advanced/makedb) subcommand, and queried with the [query](/advanced/query) subcommand.

### FAQ
#### Can you target the 16S rRNA gene instead of the default set of single copy marker genes with SingleM?
Yes. By default, SingleM builds OTU tables from protein genes rather than 16S because this in general gives more strain-level resolution due to redundancy in the genetic code. If you are really keen on using 16S, then you can use SingleM with a 16S SingleM package (spkg). There is a [repository of auxiliary packages](https://github.com/wwood/singlem_extra_packages) at which includes a 16S package that is suitable for this purpose. The resolution won't be as high taxonomically, and there are issues around copy number variation, but it could be useful to use 16S for various reasons e.g. linking it to an amplicon study or using the GreenGenes taxonomy. For now there's no 16S spkg that gets installed by default, you have to use the `--singlem-packages` flag in `pipe` mode pointing to a separately downloaded package - see [https://github.com/wwood/singlem_extra_packages](https://github.com/wwood/singlem_extra_packages). Searching for 16S reads is also much slower than searching for protein-encoding reads.

#### How should SingleM be run on multiple samples?
There are two ways. It is possible to specify multiple input files to the `singlem pipe` subcommand directly by space separating them. Alternatively `singlem pipe` can be run on each sample and OTU tables combined using `singlem summarise`. The results should be identical, though there are some performance trade-offs. For large numbers of samples (>100) it is probably preferable to run each sample individually or in smaller groups.

#### What is the difference between the num_hits and coverage columns in the OTU table generated by the pipe mode?
`num_hits` is the number of reads found from the sample in that OTU. The
`coverage` is the expected coverage of a genome with that OTU sequence i.e. the
average number of bases covering each position in a genome after read mapping.
This is calculated from `num_hits`. In particular, `num_hits` is the 'kmer
coverage' formula used by genome assembly programs, and so `coverage` is
calculated according to the following formula, adapted from the one given in
the Velvet assembler's
[manual](https://raw.githubusercontent.com/dzerbino/velvet/master/Manual.pdf):

```
coverage = num_hits * L / (L - k + 1)
```

Where `L` is the length of a read and `k` is the length of the OTU sequence including inserts and gaps (usually `60` bp).


## License
SingleM is developed by the [Woodcroft lab](https://research.qut.edu.au/cmr/team/ben-woodcroft/) at the [Centre for Microbiome Research](https://research.qut.edu.au/cmr), School of Biomedical Sciences, QUT, with contributions several including [Samuel Aroney](https://github.com/AroneyS) and [Rossen Zhao](https://github.com/rzhao-2) and many others. It is licensed under [GPL3 or later](https://gnu.org/licenses/gpl.html).

The source code is available at [https://github.com/wwood/singlem](https://github.com/wwood/singlem).

## Citation
SingleM and Sandpiper: Robust microbial taxonomic profiles from metagenomic data. Ben J Woodcroft, Samuel T. N. Aroney, Rossen Zhao, Mitchell Cunningham, Joshua A. M. Mitchell, Linda Blackall, Gene W Tyson. bioRxiv 2024.01.30.578060; doi: https://doi.org/10.1101/2024.01.30.578060
