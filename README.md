Welcome to the SingleM.

SingleM is a tool to find the abundances of discrete operational taxonomic units (OTUs) directly from shotgun metagenome data, without heavy reliance of reference sequence databases. It is able to differentiate closely related species even if those species are from lineages new to science.

Where [GraftM](https://github.com/geronimp/graftM) can give a taxonomic overview of your community e.g. proportion of a community from a particular taxonomic family, SingleM gives you the ability to answer related but distinct questions such as:

* How many different kinds of TM6 do I have?
* What is the Chao 1 diversity of my sample?
* Are the Acidobacteria in sample 1 very closely related to the Acidobacteria in sample 2?
* Do I have population genomes for the main community members?

An overview of your community can be obtained like so:
```
singlem pipe --sequences my_sequences.fastq.gz --otu_table output_otu_table.csv --threads 24
```
The output table consists of columns:
```
ribosomal_protein_S2_rpsB_gpkg  20111000_E2M    AGCCGCTGGAATCCGAAAATGCGGCCGTATATCTACGGCAAACGCAACATGATCCACATC    1       Root; k__Bacteria; p__Planctomycetes; c__Planctomycetia
ribosomal_protein_S2_rpsB_gpkg  20111000_E2M    GCCGCGTGGAACCCGAAGATGCAGCCGTACATCTACGGCAAGCGCAACGGGATCCACATC    2       Root; k__Bacteria; p__Planctomycetes
```
1. marker name
2. sample name
3. sequence of the OTU
4. number of reads detected from that OTU
5. "median" taxonomic classification of each of the reads in the OTU

Currently SingleM concentrates on 15 single copy marker genes to provide fine-grained differentiation of species that is independent of the copy-number variation issues that hamper 16S analyses. SingleM is reasonably fast and is quite scalable, although there is much room for improvement.

##License
SingleM is written by Ben Woodcroft (@wwood) and is licensed under [GPL3 or later](https://gnu.org/licenses/gpl.html).
