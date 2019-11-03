<!-- markdown-toc start - Don't edit this section. Run M-x markdown-toc-refresh-toc -->
**Table of Contents**

- [SingleM](#singlem)
    - [Generating an OTU table](#generating-an-otu-table)
    - [Further processing of OTU tables](#further-processing-of-otu-tables)
        - [Summarising OTU tables](#summarising-otu-tables)
        - [Calculating beta diversity between samples](#calculating-beta-diversity-between-samples)
        - [Creating and querying SingleM databases](#creating-and-querying-singlem-databases)
        - [Appraising assembly and genome recovery efforts](#appraising-assembly-and-genome-recovery-efforts)
        - [Installation](#installation)
            - [Installation via GNU Guix](#installation-via-gnu-guix)
            - [Installation via DockerHub](#installation-via-dockerhub)
            - [Installation via PyPI](#installation-via-pypi)
    - [Help](#help)
        - [FAQ](#faq)
            - [Can you target the 16S rRNA gene instead of the default set of ribosomal proteins with SingleM?](#can-you-target-the-16s-rrna-gene-instead-of-the-default-set-of-ribosomal-proteins-with-singlem)
            - [How should SingleM be run on multiple samples?](#how-should-singlem-be-run-on-multiple-samples)
    - [License](#license)

<!-- markdown-toc end -->
# SingleM
Welcome.

SingleM is a tool to find the abundances of discrete operational taxonomic units (OTUs) directly from shotgun metagenome data, without heavy reliance on reference sequence databases. It is able to differentiate closely related species even if those species are from lineages new to science.

Where [GraftM](https://github.com/geronimp/graftM) can give a taxonomic overview of your community e.g. proportion of a community from a particular taxonomic family, SingleM finds sequence-based OTUs from raw, untrimmed metagenomic reads.

This gives you the ability to answer questions such as:

* How many different kinds of TM6 do I have?
* What is the Chao 1 diversity of my sample?
* Are the Acidobacteria in sample 1 very closely related to the Acidobacteria in sample 2?
* Do I have population genomes for the main community members?
* How diverse are the Pelagibacteria relative to the Flavobacteria?
* Has my genome been observed in any samples submitted to the [SRA](http://www.ncbi.nlm.nih.gov/sra)?
* Which community members are more likely to assemble?

## Generating an OTU table

An overview of your community can be obtained like so:
```
singlem pipe --sequences my_sequences.fastq.gz --otu_table otu_table.csv --threads 24
```
Please use **raw** metagenome reads, not quality trimmed reads. Quality trimming with e.g. [Trimmomatic](https://doi.org/10.1093/bioinformatics/btu170) reads often makes them too short for SingleM to use. On the other hand, adapter trimming is unlikely to be detrimental.

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
### Summarising OTU tables by rarefying, clustering, etc.
Once an OTU table has been generated with the `pipe` command, it can be further processed in various ways using `summarise`:

Create a [Krona](https://sourceforge.net/p/krona/) plot of the community. The following command generates a Krona file `my_krona.html` which can be viewed in a web browser:
```
singlem summarise --input_otu_tables otu_table.csv --krona my_krona.html
```

Several OTU tables can be combined into one. Note that this is not necessary if the combined output is to be input again into summarise (or many other commands) - it is possible to just specify multiple input tables with `--input_otu_tables`. To combine:
```
singlem summarise --input_otu_tables otu_table1.csv otu_table2.csv --output_otu_table combined.otu_table.csv
```

Cluster sequences, collapsing them into OTUs with less resolution, but with more robustness against sequencing error:
```
singlem summarise --input_otu_tables otu_table.csv --cluster --clustered_output_otu_table clustered.otu_table.csv
```

Rarefy a set of OTU tables so that each sample contains the same number of OTU sequences:
```
singlem summarise --input_otu_tables otu_table.csv other_samples.otu_table.csv --rarefied_output_otu_table rarefied.otu_table.csv --number_to_choose 100
```

Conversion to [BIOM format](http://biom-format.org) for use with QIIME:
```
singlem summarise --input_otu_tables otu_table.csv other_samples.otu_table.csv --biom_prefix myprefix
```
This generates a BIOM table for each marker gene e.g. `myprefix.4.12.ribosomal_protein_L11_rplK.biom`.

### Calculating beta diversity between samples
As SingleM generates OTUs that are independent of taxonomy, they can be used as input to beta diversity methods known to be appropriate for the analysis of 16S amplicon studies, of which there are many. We recommend [express beta diversity](https://github.com/dparks1134/ExpressBetaDiversity) (EBD) as it implements many different metrics with a unified interface. For instance to calculate Bray-Curtis beta diversity, first convert your OTU table to unifrac file format using `singlem summarise`. Note that this file format does not contain any phylogenetic information, even if the format is called 'unifrac'.
```
singlem summarise --input_otu_table otu_table.csv --unifrac_by_otu otu_table.unifrac
```
The above commands generates 14 different unifrac format files, one for each marker gene used in SingleM. At this point, you need to choose one table to proceed with. Hopefully, the choice matters little, but it might pay to use multiple tables and ensure that the results are consistent.

To calculate beta diversity that does not account for the phylogenetic relationships between the OTU sequences, use the EBD script `convertToEBD.py` to convert the unifrac format into ebd format, and calculate the diversity metric:
```
convertToEBD.py otu_table.unifrac.4.12.ribosomal_protein_L11_rplK.unifrac otu_table.ebd
ExpressBetaDiversity -s otu_table.ebd -c Bray-Curtis
```
Phylogenetic tree-based methods of calculating beta diversity can also be calculated, but `pipe` must be used to generate a new OTU table using the `diamond_example` taxonomy assignment method so that each OTU is assigned to a single leaf in the tree:
```
singlem pipe --sequences my_sequences.fastq.gz --otu_table otu_table.diamond_example.csv --threads 24 --assignment_method diamond_example
```
Then, use the `--unifrac_by_taxonomy` flag to create a unifrac format file indexed by taxonomy identifier:
```
singlem summarise --otu_tables otu_table.diamond_example.csv --unifrac_by_taxonomy otu_table.diamond_example.csv
convertToEBD.py otu_table.diamond_example.unifrac otu_table.diamond_example.ebd
```
Then, finally run `ExpressBetaDiversity` using the `-t` flag.
```
ExpressBetaDiversity -s otu_table.diamond_example.ebd -c Bray-Curtis -t <path_to_tree_in_singlem_package>
```
where `<path_to_tree_in_singlem_package>` is the newick format file in the SingleM package used to find the OTU sequences. This path can be found using `singlem get_tree`.


### Creating and querying SingleM databases
It can be useful in some situations to search for sequences in OTU tables. For instance, you may ask "is the most abundant OTU or anything similar in samples B, C or D?" To answer this question make a SingleM database from sample B, C & D's OTU tables:
```
singlem makedb --otu_tables sample_B.csv sample_C.csv sample_D.csv --db_path sample_BCD.sdb
```
`.sdb` is the conventional file extension for SingleM databases. Then to query this database
```
singlem query --query_sequence TGGTCGCGGCGCTCAACCATTCTGCCCGAGTTCGTCGGCCACACCGTGGCCGTTCACAAC --db sample_BCD.sdb
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
singlem pipe --sequences raw.fq.gz --otu_table metagenome.otu_table.csv
singlem pipe --sequences my_genomes/*.fasta --otu_table genomes.otu_table.csv
singlem appraise --metagenome_otu_tables metagenome.otu_table.csv --genome_otu_tables genomes.otu_table.csv
```
One may also accommodate some sequence differences, with `--imperfect`, or
output OTU tables of those sequences that match and those that do not (see
`singlem appraise -h`). Assessing assemblies is similar to assessing genomes,
except that when assessing bins duplicate markers from the same genome are
excluded as likely contamination.

An appraisal can also be represented visually, using `appraise --plot`:

![Image of appraise](https://raw.githubusercontent.com/wwood/singlem/master/appraise_plot.png)

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


### Installation

#### Installation via GNU Guix
The most straightforward way of installing SingleM is to use the GNU Guix package which is part of the ACE Guix package collection. This method installs not just the Python libraries required but the compiled bioinformatics tools needed as well. Once you have installed Guix, clone the ACE collection and install:
```
git clone https://github.com/Ecogenomics/ace-guix
GUIX_PACKAGE_PATH=ace-guix guix package --install singlem
```
Beyond installing GNU Guix, super-user privileges are not required for installation.

#### Installation via DockerHub
A docker image generated from the Guix package is available on DockerHub. After installing Docker, run the following, replacing `[RELEASE_TAG]` with a tag from https://hub.docker.com/r/wwood/singlem/tags:
```
docker pull wwood/singlem:[RELEASE_TAG]
```
If the sequence data to be analyzed is in the current working directory, SingleM can be used like so:
```
docker run -v `pwd`:`pwd` wwood/singlem:[RELEASE_TAG] pipe --sequences `pwd`/my.fastq.gz --otu_table `pwd`/my.otu_table.csv --threads 14
```

#### Installation via conda
Conda support at this time is partial because some dependencies are not packaged
for conda, and the following is not well tested, but it may aid your
installation. After installing conda and setting up the bioconda and conda-forge
channels,

```
conda create -n singlem nose python hmmer h5py matplotlib krona diamond orfm pplacer vsearch smafa tempdir biopython biom-format dendropy mfqe
conda activate singlem
pip install orator
pip install extern
pip install squarify
pip install graftm
pip install singlem
```

#### Installation via PyPI
SingleM has migrated to Python 3. To install the Python libraries required:
```
pip install graftm
pip install singlem
```
You may need super-user privileges.

SingleM also has the following non-Python dependencies:
* [smafa](https://github.com/wwood/smafa) >= 0.5.0
* [VSEARCH](https://github.com/torognes/vsearch)

Some dependencies of [GraftM](https://github.com/geronimp/graftM):
* [OrfM](https://github.com/wwood/OrfM) >= 0.2.0 
* [HMMER](http://hmmer.janelia.org/) >= 3.1b1 
* [mfqe](https://github.com/wwood/mfqe) >= 0.5.0
* [pplacer](http://matsen.fhcrc.org/pplacer/) >= 1.1.alpha17
* [KronaTools](http://sourceforge.net/p/krona/home/krona/) >= 2.4
* [diamond](https://github.com/bbuchfink/diamond) >= 0.9

## Help
If you have any questions or comments, send a message to the [SupportM mailing list](https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/supportm) or raise a [GitHib issue](https://github.com/wwood/singlem/issues).

### FAQ
#### Can you target the 16S rRNA gene instead of the default set of ribosomal proteins with SingleM?
Yes. By default, SingleM builds OTU tables from ribosomal protein genes rather than 16S because this in general gives more strain-level resolution due to redundancy in the genetic code. If you are really keen on using 16S, then you can use SingleM with a 16S SingleM package (spkg). There is a repository of auxiliary packages at https://github.com/wwood/singlem_extra_packages including a 16S package that is suitable for this purpose. The resolution won't be as high taxonomically, and there are issues around copy number variation, but it could be useful to use 16S for various reasons e.g. linking it to an amplicon study or using the GreenGenes taxonomy. For now there's no 16S spkg that gets installed by default, you have to use the `--singlem_packages` flag in `pipe` mode pointing to a separately downloaded package - see https://github.com/wwood/singlem_extra_packages/blob/master/README.md.

#### How should SingleM be run on multiple samples?
There are two ways. It is possible to specify multiple input files to the `singlem pipe` subcommand directly by space separating them. Alternatively `singlem pipe` can be run on each sample and OTU tables combined using `singlem summarise`. The results should be identical, though there are some performance trade-offs. For large numbers of samples (>100) it is probably preferable to run each sample individually or in smaller groups.

## License
SingleM is written by [Ben Woodcroft](http://ecogenomic.org/personnel/dr-ben-woodcroft) (@wwood) at the [Australian Centre for Ecogenomics (UQ)](http://ecogenomic.org/) and is licensed under [GPL3 or later](https://gnu.org/licenses/gpl.html).
