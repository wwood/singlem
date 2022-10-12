The SingleM `summarise` subcommand transforms OTU tables and taxonomic profiles into a variety of different formats. The `summarise` subcommand is useful for visualising the results of a SingleM analysis, and for performing downstream analyses.

# Summarising OTU tables by rarefying, clustering, etc.
Once an OTU table has been generated with the `pipe` command, it can be further processed in various ways using `summarise`:

Create a [Krona](https://sourceforge.net/p/krona/) plot of the community. The following command generates a Krona file `my_krona.html` which can be viewed in a web browser:
```
singlem summarise --input-otu-tables otu_table.csv --krona my_krona.html
```

Several OTU tables can be combined into one. Note that this is not necessary if the combined output is to be input again into summarise (or many other commands) - it is possible to just specify multiple input tables with `--input-otu-tables`. To combine:
```
singlem summarise --input-otu-tables otu_table1.csv otu_table2.csv --output-otu-table combined.otu_table.csv
```

Cluster sequences, collapsing them into OTUs with less resolution, but with more robustness against sequencing error:
```
singlem summarise --input-otu-tables otu_table.csv --cluster --clustered-output-otu-table clustered.otu_table.csv
```

Rarefy a set of OTU tables so that each sample contains the same number of OTU sequences:
```
singlem summarise --input-otu-tables otu_table.csv other_samples.otu_table.csv --rarefied-output-otu-table rarefied.otu_table.csv --number-to-choose 100
```

# Calculating beta diversity between samples
Ecological [beta-diversity](https://en.wikipedia.org/wiki/Beta_diversity) metrics, which measure differences between two microbial communities, can be calculated from SingleM profiles using OTU-based or phylogenetic tree-based approaches.

## OTU-based beta diversity measures
As SingleM generates OTUs that are independent of taxonomy, they can be used as input to beta diversity methods known to be appropriate for the analysis of 16S amplicon studies, of which there are many. We recommend [express beta diversity](https://github.com/dparks1134/ExpressBetaDiversity) (EBD) as it implements many different metrics with a unified interface. For instance to calculate Bray-Curtis beta diversity, first convert your OTU table to unifrac file format using `singlem summarise`. Note that this file format does not contain any phylogenetic information, even if the format is called 'unifrac'.
```
singlem summarise --input-otu-table otu_table.csv --unifrac-by-otu otu_table.unifrac
```
The above commands generates 14 different unifrac format files, one for each marker gene used in SingleM. At this point, you need to choose one table to proceed with. Hopefully, the choice matters little, but it might pay to use multiple tables and ensure that the results are consistent.

To calculate beta diversity that does not account for the phylogenetic relationships between the OTU sequences, use the EBD script `convertToEBD.py` to convert the unifrac format into ebd format, and calculate the diversity metric:
```
convertToEBD.py otu_table.unifrac.4.12.ribosomal_protein_L11_rplK.unifrac otu_table.ebd
ExpressBetaDiversity -s otu_table.ebd -c Bray-Curtis
```
## Phylogenetic tree-based beta diversity measures
Phylogenetic tree-based methods of calculating beta diversity can also be calculated, but `pipe` must be used to generate a new OTU table using the `diamond_example` taxonomy assignment method so that each OTU is assigned to a single leaf in the tree:
```
singlem pipe --sequences my_sequences.fastq.gz --otu-table otu_table.diamond_example.csv --threads 24 --assignment-method diamond_example
```
Then, use the `--unifrac-by-taxonomy` flag to create a unifrac format file indexed by taxonomy identifier:
```
singlem summarise --otu-tables otu_table.diamond_example.csv --unifrac-by-taxonomy otu_table.diamond_example.csv
convertToEBD.py otu_table.diamond_example.unifrac otu_table.diamond_example.ebd
```
Then, finally run `ExpressBetaDiversity` using the `-t` flag.
```
ExpressBetaDiversity -s otu_table.diamond_example.ebd -c Bray-Curtis -t <path_to_tree_in_singlem_package>
```
where `<path_to_tree_in_singlem_package>` is the newick format file in the SingleM package used to find the OTU sequences. This path can be found using `singlem get_tree`.
