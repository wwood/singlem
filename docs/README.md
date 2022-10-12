[![Anaconda-Server Badge](https://anaconda.org/bioconda/singlem/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)

<img src="joel_logo.png" alt="SingleM logo" width="600"/>

# SingleM
Welcome.

SingleM is a tool for profiling shotgun metagenomes. It has a particular strength in detecting microbial lineages which are not in reference databases. The method it uses also makes it suitable for some other tasks, such as assessing eukaryotic contamination, finding bias in genome recovery and computing ecological diversity metrics.

SingleM has been applied to public metagenomes. These data are available at a companion website [sandpiper](https://sandpiper.qut.edu.au).

There are several tools, after [installation](/Installation):

* [singlem pipe](/usage/pipe) - the main workflow (`singlem pipe`) which generates OTU tables and [GTDB](https://gtdb.ecogenomic.org/) taxonomic profiles. 





## Installation



## Generating taxonomic profiles

An overview of your community can be obtained like so:
```
singlem pipe -1 <fastq_or_fasta1> -2 <fastq_or_fasta2> -p <output.profile.tsv>
```
Please use **raw** metagenome reads, not quality trimmed reads. Quality trimming with e.g. [Trimmomatic](https://doi.org/10.1093/bioinformatics/btu170) reads often makes them too short for SingleM to use. Adapter trimming is unlikely to be detrimental, but is not needed.



The output table consists of columns:
```
gene    sample  sequence        num_hits        coverage        taxonomy
4.21.ribosomal_protein_S19_rpsS my_sequences  TGGTCGCGCCGTTCGACGGTCACTCCGGACTTCATCGGCCTACAGTTCGCCGTGCACATC    1       1.64    Root; d__Bacteria; p__Proteobacteria; c__Deltaproteobacteria; o__Desulfuromonadales
4.21.ribosomal_protein_S19_rpsS my_sequences  TGGTCGCGGCGCTCAACCATTCTGCCCGAGTTCGTCGGCCACACCGTGGCCGTTCACAAC    1       1.64    Root; d__Bacteria; p__Acidobacteria; c__Solibacteres; o__Solibacterales; f__Solibacteraceae; g__Candidatus_Solibacter; s__Candidatus_Solibacter_usitatus
```
1. marker name
2. sample name
3. sequence of the OTU
4. number of reads detected from that OTU
5. estimated coverage of a genome from this OTU
6. "median" taxonomic classification of each of the reads in the OTU according to [pplacer](http://matsen.fhcrc.org/pplacer/)

Currently SingleM concentrates on 14 single copy marker genes to provide fine-grained differentiation of species that is independent of the copy-number variation issues that hamper 16S analyses. SingleM is reasonably fast and is quite scalable, although there is much room for improvement. On average, each of the 14 genes better differentiates closely related lineages than a typical 16S amplicon-based study.

## Further processing of OTU tables


### Creating and querying SingleM databases
It can be useful in some situations to search for sequences in OTU tables. For instance, you may ask "is the most abundant OTU or anything similar in samples B, C or D?" To answer this question make a SingleM database from sample B, C & D's OTU tables:
```
singlem makedb --otu-tables sample_B.csv sample_C.csv sample_D.csv --db-path sample_BCD.sdb
```
`.sdb` is the conventional file extension for SingleM databases. Then to query this database
```
singlem query --query-sequence TGGTCGCGGCGCTCAACCATTCTGCCCGAGTTCGTCGGCCACACCGTGGCCGTTCACAAC --db sample_BCD.sdb
```


### Appraising assembly and genome recovery efforts 
SingleM can be used to determine how much of a community is represented in an assembly or represented
by a set of genomes.

The assessment is carried out by comparing the set of OTU sequences in the
assembly/genomes to those found from the raw metagenome reads. A similar metric
can be estimated by the fraction of reads mapping to either the assembly or the
genome, but the approach here is more flexible and has several advantages:

1. We can determine which specific lineages are missing
2. We can match OTU sequences imperfectly, so we can e.g. make statements about whether a genus level representative genome has been recovered
3. Since the metric assesses only single copy marker genes, the appraisal is on a per-cell basis, not per-read
4. Some care must be taken, but we can prevent Eukaryotic DNA from skewing the estimate

To assess how well a set of sequences represent a metagenome, first run `pipe`
on both the genomes and the raw reads, and then use `appraise`:
```
singlem pipe --sequences raw.fq.gz --otu-table metagenome.otu_table.csv
singlem pipe --sequences my-genomes/*.fasta --otu_table genomes.otu_table.csv
singlem appraise --metagenome-otu-tables metagenome.otu_table.csv --genome-otu-tables genomes.otu_table.csv
```
One may also accommodate some sequence differences, with `--imperfect`, or
output OTU tables of those sequences that match and those that do not (see
`singlem appraise -h`). Assessing assemblies is similar to assessing genomes,
except that when assessing bins duplicate markers from the same genome are
excluded as likely contamination.

An appraisal can also be represented visually, using `appraise --plot`:

<img src="appraise_plot.png" alt="Image of appraise"/>

Each rectangle represents a single OTU sequence where its size represents its
abundance (the number of reads that OTU represents in the metagenome). Colours
represent 89% OTU clustering of these sequences (~genus level), with the
taxonomy of the most common sequence written out. Here we see that highly
abundant OTUs in SRR5040536 were assembled, but only 1 of the 3 abundant
Gallionellales OTUs has an associated bin. As is common, the highest abundance
lineages did not necessarily assemble and bin successfully. The marker
`4.20.ribosomal_protein_S15P_S13e` was chosen as the representative marker
because it has a representative fraction of OTUs binned, assembled and
unassembled.



## Help
If you have any questions or comments, send a message to the [SupportM mailing list](https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/supportm) or raise a [GitHib issue](https://github.com/wwood/singlem/issues).

### FAQ
#### Can you target the 16S rRNA gene instead of the default set of ribosomal proteins with SingleM?
Yes. By default, SingleM builds OTU tables from ribosomal protein genes rather than 16S because this in general gives more strain-level resolution due to redundancy in the genetic code. If you are really keen on using 16S, then you can use SingleM with a 16S SingleM package (spkg). There is a repository of auxiliary packages at https://github.com/wwood/singlem_extra_packages including a 16S package that is suitable for this purpose. The resolution won't be as high taxonomically, and there are issues around copy number variation, but it could be useful to use 16S for various reasons e.g. linking it to an amplicon study or using the GreenGenes taxonomy. For now there's no 16S spkg that gets installed by default, you have to use the `--singlem_packages` flag in `pipe` mode pointing to a separately downloaded package - see https://github.com/wwood/singlem_extra_packages.

#### How should SingleM be run on multiple samples?
There are two ways. It is possible to specify multiple input files to the `singlem pipe` subcommand directly by space separating them. Alternatively `singlem pipe` can be run on each sample and OTU tables combined using `singlem summarise`. The results should be identical, though there are some performance trade-offs. For large numbers of samples (>100) it is probably preferable to run each sample individually or in smaller groups.

#### What is the difference between the num_hits and coverage columns in the OTU table generated by the pipe mode?
`num_hits` is the number of reads found from the sample in that OTU. The
`coverage` is the expected coverage of a genome with that OTU sequence i.e. the
average number of bases covering each position in a genome after read mapping.
This is calculated from `num_hits`. In particular, `num_hits` is the 'kmer
coverage' formula used by genome assembly programs, and so `coverage` is
calculated accoreding to the following formula, adapted from the one given in
the Velvet assembler's
[manual](https://raw.githubusercontent.com/dzerbino/velvet/master/Manual.pdf):

```
coverage = num_hits * L / (L - k + 1)
```

Where `L` is the length of a read and `k` is the length of the OTU sequence including inserts and gaps (usually `60` bp).


## License
SingleM is written by [Ben Woodcroft](http://ecogenomic.org/personnel/dr-ben-woodcroft) (@wwood) at the [Australian Centre for Ecogenomics (UQ)](http://ecogenomic.org/) and is licensed under [GPL3 or later](https://gnu.org/licenses/gpl.html).
