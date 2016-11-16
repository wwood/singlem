Welcome to SingleM.

SingleM is a tool to find the abundances of discrete operational taxonomic units (OTUs) directly from shotgun metagenome data, without heavy reliance of reference sequence databases. It is able to differentiate closely related species even if those species are from lineages new to science.

Where [GraftM](https://github.com/geronimp/graftM) can give a taxonomic overview of your community e.g. proportion of a community from a particular taxonomic family, SingleM gives you the ability to answer related but distinct questions such as:

* How many different kinds of TM6 do I have?
* What is the Chao 1 diversity of my sample?
* Are the Acidobacteria in sample 1 very closely related to the Acidobacteria in sample 2?
* Do I have population genomes for the main community members?
* How diverse are the Pelagibacteria relative to the Flavobacteria?
* Has my genome been observed in any samples submitted to the [SRA](http://www.ncbi.nlm.nih.gov/sra)?

##Generating an OTU table
An overview of your community can be obtained like so. Please use **raw** metagenome reads, not QC'd reads. QC'ing reads often makes them too short for SingleM to use.
```
singlem pipe --sequences my_sequences.fastq.gz --otu_table otu_table.csv --threads 24
```
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

##Further processing of OTU tables
###Summarising OTU tables
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

###Calculating beta diversity between samples
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


###Creating and querying SingleM databases
It can be useful in some situations to search for sequences in OTU tables. For instance, you may ask "is the most abundant OTU or anything similar in samples B, C or D?" To answer this question make a SingleM database from sample B, C & D's OTU tables:
```
singlem makedb --otu_tables sample_B.csv sample_C.csv sample_D.csv --db_path sample_BCD.sdb
```
`.sdb` is the conventional file extension for SingleM databases. Then to query this database
```
singlem query --query_sequence TGGTCGCGGCGCTCAACCATTCTGCCCGAGTTCGTCGGCCACACCGTGGCCGTTCACAAC --db sample_BCD.sdb
```



###Installation
SingleM is not currently available on pip, though we anticipate this in future. To install, clone from the GitHub repository and run from the checked out respository:
```
git clone https://github.com/wwood/singlem
./singlem/bin/singlem -h
```

SingleM also has the following dependencies:
* [GraftM](https://github.com/geronimp/graftM), which in itself has several dependencies :( and worse, SingleM currently requires changes in GraftM introduced after 0.9.5. This will be fixed soon.
* [BLAST+](http://blast.ncbi.nlm.nih.gov/Blast.cgi)
* [VSEARCH](https://github.com/torognes/vsearch)

##License
SingleM is written by [Ben Woodcroft](http://ecogenomic.org/personnel/dr-ben-woodcroft) (@wwood) at the Australian Centre for Ecogenomics (UQ) and is licensed under [GPL3 or later](https://gnu.org/licenses/gpl.html).
