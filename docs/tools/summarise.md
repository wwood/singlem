---
title: SingleM summarise
---
# singlem summarise
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

# TAXONOMIC PROFILE INPUT

**\--input-taxonomic-profiles** *INPUT_TAXONOMIC_PROFILES* [*INPUT_TAXONOMIC_PROFILES* \...]

  Input taxonomic profiles to be e.g. converted to krona HTML, or
    concatenated

# TAXONOMIC PROFILE OUTPUT

**\--output-taxonomic-profile** FILE

  Output a single output file containing taxonomic profiles of all
    input taxonomic profile files. Requires \--input-taxonomic-profiles

**\--output-taxonomic-profile-krona** FILE

  Output taxonomic profile to this file in Krona format.

**\--output-species-by-site-relative-abundance** FILE

  Output site by species relative abundance to this file

**\--output-species-by-site-level** {species,genus,family,order,class,phylum,domain}

  Output site by species level to this file. Requires
    \--output-species-by-site-relative-abundance.

**\--output-species-by-site-relative-abundance-prefix** PATH_PREFIX

  Output site by species relative abundance to this file prefix. One
    file will be written for each taxonomic level.

**\--output-filled-taxonomic-profile** FILE

  Output a taxonomic profile where the coverage of each taxon includes
    the coverage of each of its descendent taxons e.g. the d\_\_Bacteria
    entry includes the p\_\_Patescibacteria entry.

**\--output-taxonomic-profile-with-extras** FILE

  Output a taxonomic profile with extra information (coverage,
    \'filled\' coverage, relative abundance, taxonomy level).

**\--num-decimal-places** INT

  Number of decimal places to report in the coverage column of the
    \--output-taxonomic-profile-with-extras [default: 2].

**\--output-taxonomic-level-coverage** FILE

  Output summary of how much coverage has been assigned to each
    taxonomic level in a taxonomic profile to a TSV file.

# OTU TABLE INPUT

**\--input-otu-tables**, **\--input-otu-table** *INPUT_OTU_TABLES* [*INPUT_OTU_TABLES* \...]

  Summarise these tables

**\--input-otu-tables-list** *INPUT_OTU_TABLES_LIST*

  Summarise the OTU table files newline separated in this file

**\--input-archive-otu-tables**, **\--input-archive-otu-table** *INPUT_ARCHIVE_OTU_TABLES* [*INPUT_ARCHIVE_OTU_TABLES* \...]

  Summarise these tables

**\--input-archive-otu-table-list** *INPUT_ARCHIVE_OTU_TABLE_LIST*

  Summarise the archive tables newline separated in this file

**\--input-gzip-archive-otu-table-list** *INPUT_GZIP_ARCHIVE_OTU_TABLE_LIST*

  Summarise the list of newline-separated gzip-compressed archive OTU
    tables specified in this file

**\--stream-inputs**

  Stream input OTU tables, saving RAM. Only works with
    \--output-otu-table and transformation options do not work [expert
    option].

# OTU TABLE TRANSFORMATION

**\--cluster**

  Apply sequence clustering to the OTU table. Any dashes in OTU
    sequences will be replaced by N.

**\--cluster-id** *CLUSTER_ID*

  Sequence clustering identity cutoff if \--cluster is used [default:
    0.9666666666666667 i.e. 96.66666666666667%]

**\--taxonomy** *TAXONOMY*

  Restrict analysis to OTUs that have this taxonomy (exact taxonomy or
    more fully resolved)

**\--rarefied-output-otu-table** *RAREFIED_OUTPUT_OTU_TABLE*

  Output rarefied output OTU table, where each gene and sample
    combination is rarefied

**\--number-to-choose** *NUMBER_TO_CHOOSE*

  Rarefy using this many sequences. Sample/gene combinations with an
    insufficient number of sequences are ignored with a warning
    [default: maximal number such that all samples have sufficient
    counts]

**\--collapse-to-sample-name** *COLLAPSE_TO_SAMPLE_NAME*

  Merge all OTUs into a single OTU table, using the given sample name.
    Requires archive OTU table input and output.

**\--collapse-coupled**

  Merge forward and reverse read OTU tables into a unified table.
    Sample names of coupled reads must end in \'1\' and \'2\'
    respectively. Read names are ignored, so that if the forward and
    reverse from a pair contain the same OTU sequence, they will each
    count separately.

**\--collapse-paired-with-unpaired-archive-otu-table** *COLLAPSE_PAIRED_WITH_UNPAIRED_ARCHIVE_OTU_TABLE*

  For archive OTU tables that have both paired and unpaired
    components, merge these into a single output archive OTU table

# OTU TABLE OUTPUT

**\--output-otu-table** *OUTPUT_OTU_TABLE*

  Output combined OTU table to this file

**\--output-archive-otu-table** *OUTPUT_ARCHIVE_OTU_TABLE*

  Output combined OTU table to this file

**\--output-translated-otu-table** *OUTPUT_TRANSLATED_OTU_TABLE*

  Output combined OTU table to this file, with seqeunces translated
    into amino acids

**\--output-extras**

  Output extra information in the standard output OTU table

**\--krona** *KRONA*

  Name of krona file to generate. Note that this generates a krona
    file from the OTU table, not the taxonomic profile

**\--wide-format-otu-table** *WIDE_FORMAT_OTU_TABLE*

  Name of output species by site CSV file

**\--strain-overview-table** *STRAIN_OVERVIEW_TABLE*

  Name of output strains table to generate

**\--unifrac-by-otu** *UNIFRAC_BY_OTU*

  Output UniFrac format file where entries are OTU sequences

**\--unifrac-by-taxonomy** *UNIFRAC_BY_TAXONOMY*

  Output UniFrac format file where entries are taxonomies (generally
    used for phylogeny-driven beta diversity when pipe was run with
    \'\--assignment_method diamond_example\')

**\--clustered-output-otu-table** *CLUSTERED_OUTPUT_OTU_TABLE*

  Output an OTU table with extra information about the clusters. To
    simply cluster an OTU table, use \--cluster with \--output-otu-table
    instead.

**\--exclude-off-target-hits**

  Exclude hits that are not in the target domain of each SingleM
    package

**\--singlem-packages** *SINGLEM_PACKAGES* [*SINGLEM_PACKAGES* \...]

  Packages used in the creation of the OTU tables

**\--metapackage** *METAPACKAGE*

  Metapackage used in the creation of the OTU tables

**\--unaligned-sequences-dump-file** *UNALIGNED_SEQUENCES_DUMP_FILE*

  Output unaligned sequences from in put archive OTU table to this
    file. After each read name \'\~N\' is added which corresponds to the
    order of the read in the archive OTU table, so that no two sequences
    have the same read name. N\>1 can happen e.g. when the input file
    contains paired reads. \~0 does not necessarily correspond to the
    first read in the original input sequence set, but instead to the
    order in the input archive OTU table.

# OTHER GENERAL OPTIONS

**\--debug**

  output debug information

**\--version**

  output version information and quit

**\--quiet**

  only output errors

**\--full-help**

  print longer help message

**\--full-help-roff**

  print longer help message in ROFF (manpage) format

# AUTHORS

>     Ben J. Woodcroft, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
>     Samuel Aroney, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
>     Raphael Eisenhofer, Centre for Evolutionary Hologenomics, University of Copenhagen, Denmark
>     Rossen Zhao, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology

# EXAMPLES

Convert a taxonomic profile to a site-by-species table at the genus level:

  **\$ singlem summarise \--input-taxonomic-profiles \<profile1.tsv\>
    \<profile2.tsv\> \--output-species-by-site-relative-abundance
    \<output.tsv\> \--output-species-by-site-level genus**

Create a Krona diagram from a taxonomic profile:

  **\$ singlem summarise \--input-taxonomic-profiles \<profile1.tsv\>
    \--output-taxonomic-profile-krona \<output.html\>**

Add extra coverage and relative abundance information to a taxonomic profile:

  **\$ singlem summarise \--input-taxonomic-profiles \<profile1.tsv\>
    \--output-taxonomic-profile-with-extras \<output.tsv\>**
