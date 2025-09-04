---
title: SingleM prokaryotic_fraction
---
# singlem prokaryotic_fraction
The SingleM `prokaryotic_fraction` subcommand (also available as `microbial_fraction`) estimates the fraction of reads in a metagenome that are microbial, compared to everything else e.g. eukaryote- or phage-derived. Here we define 'microbial' as either bacterial or archaeal, including their plasmids. It can help prioritise samples for deeper sequencing, forecast how much sequencing is required for a given set of samples, and identify problematic steps in genome-resolved metagenomic pipelines, for instance.

SingleM `prokaryotic_fraction` also estimates the average genome size (AGS) of microbial cells in the sample.

The main conceptual advantage of this method over other tools is that it does not require reference sequences of the non-microbial genomes that may be present (e.g. those of an animal host). Instead, it uses a SingleM taxonomic profile of the metagenome to "add up" the components of the community which are microbial. The remaining components are non-microbial e.g. host, diet, or phage. 

Roughly, the number of microbial bases is estimated by summing the genome sizes of each species/taxon after multiplying by their coverage in a taxonomic profile. The microbial fraction is then calculated as the ratio of microbial bases to the total metagenome size.

## Usage

To run `prokaryotic_fraction`, first run `pipe` on your metagenome.

```bash
$ singlem pipe --forward SRR9841429_1.fastq.gz --reverse SRR9841429_2.fastq.gz --threads 32 -p SRR9841429.profile
```
Then run `prokaryotic_fraction` on the profile.

```bash
singlem prokaryotic_fraction --forward SRR9841429_1.fastq.gz --reverse SRR9841429_2.fastq.gz -p SRR9841429.profile >SRR9841429.smf.tsv
```

The output you get is a tab-separated values file containing (with some columns omitted):

| sample        | bacterial_archaeal_bases  | metagenome_size | read_fraction | .. |
| ------------- | ------------------------- | --------------- | ------------- | -- |
|  SRR9841429_1 |                2059301599 |      2648177782 |         77.76 | .. |

So, for this sample, 77.76% of the reads are estimated to be microbial. 

A full table is shown below for this one sample, transposed for formatting:

| header | value |
| --- | --- |
| sample | SRR9841429_1 |
| bacterial_archaeal_bases | 2059301599 |
| metagenome_size | 2648177782 |
| read_fraction | 77.76 |
| average_bacterial_archaeal_genome_size | 3448186 |
| warning | |
| domain_relative_abundance | 1.10 |
| phylum_relative_abundance | 0.65 |
| class_relative_abundance | 1.72 |
| order_relative_abundance | 0.54 |
| family_relative_abundance | 1.03 |
| genus_relative_abundance | 8.01 |
| species_relative_abundance | 86.94 |

This table tells us that:

1. SingleM estimates 77.76% of the metagenome's reads are microbial (bacterial or archaeal).
2. The total number of bases assigned to microbial genomes is 2,059,301,599.
3. The total number of bases in the metagenome is 2,648,177,782.
4. The estimated average genome size of the microbial cells in the sample is 3,448,186 bp.
5. No warning about unreliability due to high prevalence of novel lineages was emitted (see below). The `warning` column is empty if no warning was emitted.
6. 1.10% of the profile is classified as domain level _and no further_. So 100-1.1=98.9% of the profile is classified to at least phylum level.
7. 0.65% of the profile is classified as phylum level _and no further_. So 100-1.1-0.65=98.25% of the profile is classified to at least class level.
8. and so on down to species level. 86.94% of the profile is classified to species level.

Note that the `read_fraction` and `relative_abundance` columns are percentages, i.e. 0.3 is 0.3%, not 30%. 

## Warnings and limitations

The method is least reliable in simple communities consisting of a small number of species that are missing from the reference database. The main challenge in these cases is that the genome sizes of novel species are hard to estimate accurately. Multiplying a coverage from the taxonomic profile against an uncertain genome length equals an uncertain number of bases assigned as microbial.

To detect such situations SingleM `prokaryotic_fraction` emits a warning when a sample's prokaryotic_fraction estimate could under- or overestimate the read fraction. Specifically, the warning is emitted when the 3 highest abundance lineages not classified to the species level would change the estimated read fraction of the sample by >2% if their genome size is halved or doubled. Prokaryotic fractions are also capped at 100% since values greater than this are impossible, but the original value can be recovered from the output if you calculate `bacterial_archaeal_bases / metagenome_size`.

One current limitation of the approach relates to multi-copy plasmids. In `prokaryotic_fraction`, the genome size of each microbial species is estimated as the sum of the chromosome and plasmid sizes, since these are the sequences available for each genome. However, in a metagenome, a species' plasmid may occur in multiple copies per cell (e.g. if the plasmid is 'high copy number'). SingleM `prokaryotic_fraction` does not account for plasmid copy number, leading to an underestimation of the prokaryotic fraction when plasmids are multi-copy. However, we consider this to be a minor issue, since plasmids are typically small compared to chromosomes. The average genome size estimate is unaffected by this limitation since by definition each plasmid counts only once regardless of its copy number.
# OPTIONS

# INPUT

**-p**, **\--input-profile** *INPUT_PROFILE*

  Input taxonomic profile file [required]

# READ INFORMATION [1+ ARGS REQUIRED]

**-1**, **\--forward**, **\--reads**, **\--sequences** sequence_file [sequence_file \...]

  nucleotide read sequence(s) (forward or unpaired) to be searched.
    Can be FASTA or FASTQ format, GZIP-compressed or not. These must be
    the same ones that were used to generate the input profile.

**-2**, **\--reverse** sequence_file [sequence_file \...]

  reverse reads to be searched. Can be FASTA or FASTQ format,
    GZIP-compressed or not. These must be the same reads that were used
    to generate the input profile.

**\--input-metagenome-sizes** *INPUT_METAGENOME_SIZES*

  TSV file with \'sample\' and \'num_bases\' as a header, where sample
    matches the input profile name, and num_reads is the total number
    (forward+reverse) of bases in the metagenome that was analysed with
    \'pipe\'. These must be the same reads that were used to generate
    the input profile.

# DATABASE

**\--taxon-genome-lengths-file** *TAXON_GENOME_LENGTHS_FILE*

  TSV file with \'rank\' and \'genome_size\' as headers [default: Use
    genome lengths from the default metapackage]

**\--metapackage** *METAPACKAGE*

  Metapackage containing genome lengths [default: Use genome lengths
    from the default metapackage]

# OTHER OPTIONS

**\--accept-missing-samples**

  If a sample is missing from the input-metagenome-sizes file, skip
    analysis of it without croaking.

**\--output-tsv** *OUTPUT_TSV*

  Output file [default: stdout]

**\--output-per-taxon-read-fractions** *OUTPUT_PER_TAXON_READ_FRACTIONS*

  Output a fraction for each taxon to this TSV [default: Do not
    output anything]

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
