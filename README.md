Welcome to SingleM.

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

Currently SingleM concentrates on 15 single copy marker genes to provide fine-grained differentiation of species that is independent of the copy-number variation issues that hamper 16S analyses. SingleM is reasonably fast and is quite scalable, although there is much room for improvement. On average, each of the 15 genes better differentiates closely related lineages than a typical 16S amplicon-based study.

## Further processing of OTU tables
### Summarising OTU tables
Once an OTU table has been generated with the `pipe` command, it can be further processed in various ways using `summarise`:

Create a [Krona](https://sourceforge.net/p/krona/) plot of the community. The following command generates `my_krona*.html` files which can be viewed in a web browser:
```
singlem summarise --input_otu_table otu_table.csv --krona my_krona
```

Cluster sequences, collapsing them into OTUs with less resolution, but with more robustness against sequencing error:
```
singlem summarise --input_otu_table otu_table.csv --cluster --clustered_output_otu_table clustered.otu_table.csv
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
As SingleM generates OTUs that are independent of taxonomy, they can be used as input to beta diversity methods known to be appropriate for the analysis of 16S amplicon studies, of which there are many. We recommend [express beta diversity](https://github.com/dparks1134/ExpressBetaDiversity) (EBD) as it implements many different metrics with a unified interface. For instance to calculate Bray-Curtis beta diversity, first convert your OTU table to unifrac format using `singlem summarise`:
```
singlem summarise --input_otu_table otu_table.csv --unifrac otu_table.unifrac
```
The above commands generates 15 different unifrac format files, one for each marker gene used in SingleM. At this point, you need to choose one table to proceed with. Hopefully, the choice matters little, but it might pay to use multiple tables and ensure that the results are consistent.

To calculate beta diversity, use the EBD script `convertToEBD.py` to convert the unifrac format into ebd format, and calculate the diversity metric:
```
convertToEBD.py otu_table.unifrac.4.12.ribosomal_protein_L11_rplK.unifrac otu_table.ebd
ExpressBetaDiversity -s otu_table.ebd -c Bray-Curtis
```
Phylogenetic tree-based methods of calculating beta diversity can also be calculated, but `pipe` must be used to generate a new OTU table using the `diamond_example` taxonomy assignment method so that each OTU is assigned to a single leaf in the tree:
```
singlem pipe --sequences my_sequences.fastq.gz --otu_table otu_table.diamond_example.csv --threads 24 --assignment_method diamond_example
singlem summarise --otu_tables otu_table.diamond_example.csv --unifrac otu_table.diamond_example.csv
convertToEBD.py otu_table.diamond_example.unifrac otu_table.diamond_example.ebd
ExpressBetaDiversity -s otu_table.diamond_example.ebd -c Bray-Curtis -t `singlem get_tree --marker_name 4.21.ribosomal_protein_S19_rpsS`
```


### Creating and querying SingleM databases
It can be useful in some situations to search for sequences in OTU tables. For instance, you may ask "is the most abundant OTU or anything similar in samples B, C or D?" To answer this question make a SingleM database from sample B, C & D's OTU tables:
```
singlem makedb --otu_tables sample_B.csv sample_C.csv sample_D.csv --db_path sample_BCD.sdb
```
`.sdb` is the conventional file extension for SingleM databases. Then to query this database
```
singlem query --query_sequence TGGTCGCGGCGCTCAACCATTCTGCCCGAGTTCGTCGGCCACACCGTGGCCGTTCACAAC --db sample_BCD.sdb
```


### Appraising genome recovery efforts
To assess how well a set of genomes (or population genomes) represent those present in a metagenome, first run `pipe` on both the genomes and the raw reads, and then use `appraise`:
```
singlem appraise --metagenome_otu_tables metagenome.otu_table.csv --genome_otu_tables genomes.otu_table.csv
```
One may also accommodate some sequence differences, with `--imperfect`, or output OTU tables of OTUs in the genomes or not in the genomes with `--accounted_for_otu_table` and `--unaccounted_for_otu_table`.



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
docker run -v `pwd`:`pwd` wwood/singlem singlem pipe --sequences `pwd`/my.fastq.gz --otu_table `pwd`/my.otu_table.csv --threads 15
```

#### Installation via PyPI
To install the Python libraries required:
```
pip install graftm
pip install singlem
```
You may need super-user privileges.

SingleM also has the following non-Python dependencies:
* [BLAST+](http://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [VSEARCH](https://github.com/torognes/vsearch)

Some dependencies of [GraftM](https://github.com/geronimp/graftM):
* [OrfM](https://github.com/wwood/OrfM) >= 0.2.0 
* [HMMER](http://hmmer.janelia.org/) >= 3.1b1 
* [fxtract](https://github.com/ctSkennerton/fxtract)
* [pplacer](http://matsen.fhcrc.org/pplacer/) >= 1.1.alpha17
* [KronaTools](http://sourceforge.net/p/krona/home/krona/) >= 2.4
* [diamond](https://github.com/bbuchfink/diamond) >= 0.9

## Help
If you have any questions or comments, send a message to the [SupportM mailing list](https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/supportm) or raise a [GitHib issue](https://github.com/wwood/singlem/issues).

## License
SingleM is written by [Ben Woodcroft](http://ecogenomic.org/personnel/dr-ben-woodcroft) (@wwood) at the [Australian Centre for Ecogenomics (UQ)](http://ecogenomic.org/) and is licensed under [GPL3 or later](https://gnu.org/licenses/gpl.html).
