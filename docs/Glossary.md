# Glossary

## **Taxonomic profile**
A tab-separated table containing the estimated abundances of GTDB taxons in a metagenome. It is in TSV format with 3 columns, with each row corresponding to a taxon. A taxonomic profile may also be called a **condensed profile**, since it is the output of the `condense` algorithm within the main `pipe` workflow. Taxonomic profiles can be converted to other formats using `singlem summarise`. Columns:
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

## **OTU table**
A table containing window sequences per metagenome/contig and marker gene. It may be in default form (a TSV with 6 columns, like below), or an extended form with more detail in further columns. The default OTU table output from [pipe](/tools/pipe), [renew](/tools/renew) and [summarise](/tools/summarise) subcommands has 6 columns, with one sequence per row. The extended form OTU table and archive OTU tables have further information (see below). Columns of a default OTU table:
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

## **OTU table (extended form)**
The extended OTU table form generated with the `--output-extras` option to the [pipe](/tools/pipe), [renew](/tools/renew) and [summarise](/tools/summarise) subcommands, has all the columns of a regular OTU table, but with several additional columns which contain more information about each OTU:
  1. read_names - the names of the reads which encode the OTU sequence
  2. nucleotides_aligned - the number of nucleotides which aligned to the window (usually 60, but can be more or less if there are gaps or inserts)
  3. taxonomy_by_known? - whether the taxonomy of the OTU was determined by known genomes (TRUE) or by the reads themselves (FALSE). Currently this is a disused column and is always marked FALSE.
  4. read_unaligned_sequences - the raw sequences of the reads which encode the OTU sequence
  5. equal_best_hit_taxonomies - the taxonomies of the best hits to the OTU sequence, if there are multiple equally good hits. This is a JSON array of strings.


## **Archive OTU table**
Similar to an extended form OTU table, but in JSON form for machine readability and with formatting version recorded. The [renew](/tools/renew) subcommand which re-analyses a dataset requires this format of OTU table rather than the default tab-separated OTU table format. The canonical file extension for SingleM packages is `.json`.

## **SingleM package (spkg)** 
Reference data for one particular marker gene and its window position. The canonical file extension for SingleM packages is `.spkg`.

## **SingleM metapackage (smpkg)** 
A collection of SingleM packages, with additional indices. The canonical file extension for SingleM metapackages is `.smpkg`.

## **SingleM database** 
An OTU table which has been converted to SQLite3 format and sequence similarity search indexes. Canonically SingleM databases are named with the `.sdb` extension, but this is not enforced. SingleM databases are created with the [makedb](/advanced/makedb) subcommand, and queried with the [query](/advanced/query) subcommand.
