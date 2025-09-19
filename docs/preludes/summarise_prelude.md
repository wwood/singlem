The SingleM `summarise` subcommand transforms taxonomic profiles and OTU tables into a variety of different formats. The `summarise` subcommand is useful for transforming the default output formats of `pipe`, visualising the results of a SingleM analysis, and for performing some downstream analyses.

# Converting a taxonomic profile to other formats

By default, the `pipe` command outputs a taxonomic profile in the format of a CSV file, when run with `-p`/`--taxonomic-profile`, which looks like this:
```
sample     coverage  taxonomy
marine0.1  3.64      Root; d__Archaea
marine0.1  0.02      Root; d__Bacteria
marine0.1  0.56      Root; d__Archaea; p__Thermoproteota
marine0.1  0.8       Root; d__Bacteria; p__Desulfobacterota
marine0.1  2.17      Root; d__Bacteria; p__Proteobacteria
...
```
The coverage is the estimated per-base read coverage of a species genome (or, for higher level taxons, the sum of their constituent species). Importantly, it is not inclusive of its descendents' coverage. For instance, the coverage of `Root; d__Bacteria` (`0.02`) does not include the coverage assigned to `p__Desulfobacterota` (`0.8`) or `p__Proteobacteria` (`2.17`). For more information on coverage, see the [Glossary](/Glossary#coverage-unfilled-coverage-and-filled-coverage).

For many applications, this format is inconvenient, so `summarise` mode provides some conversion options.

## Conversion to krona
```
singlem summarise --input-taxonomic-profile doco_example.profile \
    --output-taxonomic-profile-krona doco_example.profile.html
```
This will generate a new [krona](https://github.com/marbl/Krona) diagram, a hierarchical graphical representation of your community. Krona diagrams can also be output directly from pipe, using the `--taxonomic-profile-krona` option.

## Relative abundance and species-by-site
A commonly requested format for these data is the relative abundance, expressed as a species-by-site table. You can generate this from the above profile by running the below, which generates it on the phylum level. 
```
singlem summarise --input-taxonomic-profile doco_example.profile \
    --output-species-by-site-relative-abundance doco_example.species_by_site.csv \
    --output-species-by-site-level phylum
```
Where `doco_example.species_by_site.csv` is then generated as below, where each number is the relative abundance expressed as a percentage:
```
taxonomy                                marine0.1
unassigned                              50.9
Root; d__Archaea; p__Thermoproteota     7.79
Root; d__Bacteria; p__Desulfobacterota  11.13
Root; d__Bacteria; p__Proteobacteria    30.18
```
If there are multiple samples in the input files, then there will be multiple columns, one for each sample.

If you wish to generate a species-by-site table for another taxonomic level e.g. genus, you can specify the `--output-species-by-site-level` option:
```
singlem summarise --input-taxonomic-profile doco_example.profile \
    --output-species-by-site-relative-abundance doco_example.species_by_site.csv \
    --output-species-by-site-level genus
```

If you wish to generate a species-by-site table for each taxonomic level, you can use `--output-species-by-site-relative-abundance-prefix`:
```
singlem summarise --input-taxonomic-profile doco_example.profile \
    --output-species-by-site-relative-abundance-prefix myprefix
```
Usually, this will generate a different file for each taxonomic level (domain, phylum, .., species), with the prefix you specified. However, since the table here has only 2 taxonomic levels (domain and phylum), it will only generate those two files:
```
myprefix-domain.tsv
myprefix-phylum.tsv
```
More files will usually be generated (all the way down to species level), but the example only contains domain and phylum level taxonomic assignments.

## Long form with extras
In some cases, it is more convenient to keep the long form, but add some additional columns, including relative abundance and the [filled coverage](/Glossary#coverage-unfilled-coverage-and-filled-coverage):
```
singlem summarise --input-taxonomic-profile doco_example.profile \
    --output-taxonomic-profile-with-extras doco_example.with_extras.tsv
```
The generated `doco_example.with_extras.tsv` contains
```
sample     coverage  full_coverage  relative_abundance  level   taxonomy
marine0.1  0         7.19           100.0               root    Root
marine0.1  3.64      4.2            58.41               domain  Root; d__Archaea
marine0.1  0.02      2.99           41.59               domain  Root; d__Bacteria
marine0.1  0.56      0.56           7.79                phylum  Root; d__Archaea; p__Thermoproteota
marine0.1  0.8       0.8            11.13               phylum  Root; d__Bacteria; p__Desulfobacterota
marine0.1  2.17      2.17           30.18               phylum  Root; d__Bacteria; p__Proteobacteria
```

# Summarising OTU tables
The following describes how summarise can be used to transform [OTU tables](/Glossary#otu-table), rather than [taxonomomic profiles](/Glossary#taxonomic-profile). 

## Summarising OTU tables by rarefying, clustering, etc.
Once an OTU table has been generated with the `pipe` command, it can be further processed in various ways using `summarise`:

Create a [Krona](https://sourceforge.net/p/krona/) plot of the community. The following command generates a Krona file `my_krona.html` which can be viewed in a web browser:
```
singlem summarise --input-otu-tables otu_table.csv --krona my_krona.html
```

Several OTU tables can be combined into one. Note that this is not necessary if the combined output is to be input again into summarise (or many other commands) - it is possible to just specify multiple input tables with `--input-otu-tables`. To combine:
```
singlem summarise --input-otu-tables otu_table1.csv otu_table2.csv \
    --output-otu-table combined.otu_table.csv
```

Cluster sequences, collapsing them into OTUs with less resolution, but with more robustness against sequencing error:
```
singlem summarise --input-otu-tables otu_table.csv --cluster \
    --output-otu-table clustered.otu_table.csv
```

The `--clustered-output-otu-table` option can be used to output a clustered OTU table that includes information about which sequences have been clustered together:
```
singlem summarise --input-otu-tables otu_table.csv --cluster \
    --clustered-output-otu-table clustered_with_details.otu_table.csv
```

Rarefy a set of OTU tables so that each sample contains the same number of OTU sequences:
```
singlem summarise --input-otu-tables otu_table.csv \
    other_samples.otu_table.csv --rarefied-output-otu-table \
    rarefied.otu_table.csv --number-to-choose 100
```

## Calculating beta diversity between samples
Ecological [beta-diversity](https://en.wikipedia.org/wiki/Beta_diversity) metrics, which measure differences between two microbial communities, can be calculated from SingleM profiles using OTU-based approaches.

As SingleM generates OTUs that are independent of taxonomy, they can be used as input to beta diversity methods known to be appropriate for the analysis of 16S amplicon studies, of which there are many. We recommend [express beta diversity](https://github.com/dparks1134/ExpressBetaDiversity) (EBD) as it implements many different metrics with a unified interface. For instance to calculate Bray-Curtis beta diversity, first convert your OTU table to unifrac file format using `singlem summarise`. Note that this file format does not contain any phylogenetic information, even if the format is called 'unifrac'.
```
singlem summarise --input-otu-table otu_table.csv --unifrac-by-otu \
    otu_table-
```
The above commands generates a different unifrac format file for each marker. At this point, you need to choose one table to proceed with. It would probably be best to choose a marker that targets both Bacteria and Archaea e.g. `S3.5.ribosomal_protein_S2_rpsB`. Beyond choosing a marker that targets both domains, hopefully the choice matters little, but it might pay to use multiple tables and ensure that the results are consistent. 

To calculate beta diversity that does not account for the phylogenetic relationships between the OTU sequences, use the EBD script `convertToEBD.py` to convert the unifrac format into ebd format, and calculate the diversity metric:
```
convertToEBD.py otu_table-S3.5.ribosomal_protein_S2_rpsB.unifrac \
    otu_table.ebd
ExpressBetaDiversity -s otu_table.ebd -c Bray-Curtis
```
