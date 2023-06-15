---
title: SingleM summarise
---
# singlem summarise
The SingleM `summarise` subcommand transforms OTU tables and taxonomic profiles into a variety of different formats. The `summarise` subcommand is useful for visualising the results of a SingleM analysis, and for performing downstream analyses.

# Summarising OTU tables by rarefying, clustering, etc.
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
    --clustered-output-otu-table clustered.otu_table.csv
```

Rarefy a set of OTU tables so that each sample contains the same number of OTU sequences:
```
singlem summarise --input-otu-tables otu_table.csv \
    other_samples.otu_table.csv --rarefied-output-otu-table \
    rarefied.otu_table.csv --number-to-choose 100
```

# Calculating beta diversity between samples
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

INPUT
=====

**\--input-otu-tables**, **\--input-otu-table** *INPUT_OTU_TABLES* [*INPUT_OTU_TABLES* \...]

  Summarise these tables

**\--input-otu-tables-list** *INPUT_OTU_TABLES_LIST*

  Summarise the OTU table files newline separated in this file

**\--input-archive-otu-tables**, **\--input-archive-otu-table** *INPUT_ARCHIVE_OTU_TABLES* [*INPUT_ARCHIVE_OTU_TABLES* \...]

  Summarise these tables

**\--input-gzip-archive-otu-table-list** *INPUT_GZIP_ARCHIVE_OTU_TABLE_LIST*

  Summarise the list of newline-separated gzip-compressed archive OTU
    tables specified in this file

**\--input-archive-otu-table-list** *INPUT_ARCHIVE_OTU_TABLE_LIST*

  Summarise the archive tables newline separated in this file

**\--stream-inputs**

  Stream input OTU tables, saving RAM. Only works with
    \--output-otu-table and transformation options do not work [expert
    option].

**\--input-taxonomic-profiles** *INPUT_TAXONOMIC_PROFILES* [*INPUT_TAXONOMIC_PROFILES* \...]

  Convert these taxonomic profiles to krona HTML, output specified by
    \--output-taxonomic-profile-krona

TRANSFORMATION
==============

**\--cluster**

  Apply sequence clustering to the OTU table

**\--cluster-id** *CLUSTER_ID*

  Sequence clustering identity cutoff if \--cluster is used

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

OUTPUT
======

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

  Name of krona file to generate

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

  Output an OTU table with extra information about the clusters

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

**\--output-taxonomic-profile-krona** *OUTPUT_TAXONOMIC_PROFILE_KRONA*

  Output taxonomic profile to this file in Krona format. Requires
    \--input-taxonomic-profiles

OTHER GENERAL OPTIONS
=====================

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

AUTHORS
=======

>     Ben J. Woodcroft, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
>     Samuel Aroney, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
>     Rossen Zhao, Centre for Microbiome Research, School of Biomedical Sciences, Faculty of Health, Queensland University of Technology
