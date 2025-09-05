---
title: Lyrebird pipe
---
# lyrebird pipe

# DESCRIPTION

Generate a taxonomic profile or OTU table for dsDNA phages from raw
sequences

# OPTIONS

# COMMON OPTIONS

**-1**, **\--forward**, **\--reads**, **\--sequences** sequence_file [sequence_file \...]

  nucleotide read sequence(s) (forward or unpaired) to be searched.
    Can be FASTA or FASTQ format, GZIP-compressed or not.

**-2**, **\--reverse** sequence_file [sequence_file \...]

  reverse reads to be searched. Can be FASTA or FASTQ format,
    GZIP-compressed or not.

**-f**, **\--genome-fasta-files** PATH [PATH \...]

  Path(s) to genome FASTA files. These are processed like input given
    with \--forward, but use higher default values for
    \--min-taxon-coverage and \--min-orf-length.

**-d**, **\--genome-fasta-directory** PATH

  Directory containing genome FASTA files. Treated identically to
    \--forward input with higher default values for
    \--min-taxon-coverage and \--min-orf-length.

**\--genome-fasta-list** PATH

  File containing genome FASTA paths, one per line. Behaviour matches
    \--forward with higher default values for \--min-taxon-coverage and
    \--min-orf-length.

**-x**, **\--genome-fasta-extension** EXT

  File extension of genomes in the directory specified with
    -d/\--genome-fasta-directory. [default: fna]

**-p**, **\--taxonomic-profile** FILE

  output a \'condensed\' taxonomic profile for each sample based on
    the OTU table. Taxonomic profiles output can be further converted to
    other formats using singlem summarise.

**\--taxonomic-profile-krona** FILE

  output a \'condensed\' taxonomic profile for each sample based on
    the OTU table

**\--otu-table** filename

  output OTU table

**\--threads** num_threads

  number of CPUS to use [default: 1]

**\--assignment-method** {smafa_naive_then_diamond,scann_naive_then_diamond,annoy_then_diamond,scann_then_diamond,diamond,diamond_example,annoy,pplacer}

  Method of assigning taxonomy to OTUs and taxonomic profiles
    [default: smafa_naive_then_diamond]

| Method                   | Description                                                                                                                                                                                                                                   |
|:-------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| smafa_naive_then_diamond | Search for the most similar window sequences \<= 3bp different using a brute force algorithm (using the smafa implementation) over all window sequences in the database, and if none are found use DIAMOND blastx of all reads from each OTU. |
| scann_naive_then_diamond | Search for the most similar window sequences \<= 3bp different using a brute force algorithm over all window sequences in the database, and if none are found use DIAMOND blastx of all reads from each OTU.                                  |
| annoy_then_diamond       | Same as scann_naive_then_diamond, except search using ANNOY rather than using brute force. Requires a non-standard metapackage.                                                                                                               |
| scann_then_diamond       | Same as scann_naive_then_diamond, except search using SCANN rather than using brute force. Requires a non-standard metapackage.                                                                                                               |
| diamond                  | DIAMOND blastx best hit(s) of all reads from each OTU.                                                                                                                                                                                        |
| diamond_example          | DIAMOND blastx best hit(s) of all reads from each OTU, but report the best hit as a sequence ID instead of a taxonomy.                                                                                                                        |
| annoy                    | Search for the most similar window sequences \<= 3bp different using ANNOY, otherwise no taxonomy is assigned. Requires a non-standard metapackage.                                                                                           |
| pplacer                  | Use pplacer to assign taxonomy of each read in each OTU. Requires a non-standard metapackage.                                                                                                                                                 |

**\--output-extras**

  give extra output for each sequence identified (e.g. the read(s)
    each OTU was generated from) in the output OTU table [default: not
    set]

# LESS COMMON OPTIONS

**\--archive-otu-table** filename

  output OTU table in archive format for making DBs etc. [default:
    unused]

**\--metapackage** *METAPACKAGE*

  Set of SingleM packages to use [default: use the default set]

**\--sra-files** sra_file [sra_file \...]

  \"sra\" format files (usually from NCBI SRA) to be searched

**\--read-chunk-size** num_reads

  Size chunk to process at a time (in number of reads). Requires
    \--sra-files.

**\--read-chunk-number** chunk_number

  Process only this specific chunk number (1-based index). Requires
    \--sra-files.

**\--output-jplace** filename

  Output a jplace format file for each singlem package to a file
    starting with this string, each with one entry per OTU. Requires
    \'pplacer\' as the \--assignment_method [default: unused]

**\--singlem-packages** *SINGLEM_PACKAGES* [*SINGLEM_PACKAGES* \...]

  SingleM packages to use [default: use the set from the default
    metapackage]

**\--assignment-singlem-db** *ASSIGNMENT_SINGLEM_DB*

  Use this SingleM DB when assigning taxonomy [default: not set, use
    the default]

**\--diamond-taxonomy-assignment-performance-parameters** *DIAMOND_TAXONOMY_ASSIGNMENT_PERFORMANCE_PARAMETERS*

  Performance-type arguments to use when calling \'diamond blastx\'
    during the taxonomy assignment step. [default: use setting defined
    in metapackage when set, otherwise use \'\--block-size 0.5
    \--target-indexed -c1\']

**\--evalue** *EVALUE*

  HMMSEARCH e-value cutoff to use for sequence gathering [default:
    1e-05]

**\--min-orf-length** length

  When predicting ORFs require this many base pairs uninterrupted by a
    stop codon [default: 72 for reads, 300 for genomes]

**\--restrict-read-length** length

  Only use this many base pairs at the start of each sequence searched
    [default: no restriction]

**\--translation-table** number

  Codon table for translation. By default, translation table 4 is
    used, which is the same as translation table 11 (the usual
    bacterial/archaeal one), except that the TGA codon is translated as
    tryptophan, not as a stop codon. Using table 4 means that the
    minority of organisms which use table 4 are not biased against,
    without a significant effect on the majority of bacteria and archaea
    that use table 11. See
    http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes
    for details on specific tables. [default: 4]

**\--filter-minimum-protein** length

  Ignore reads aligning in less than this many positions to each
    protein HMM when using \--no-diamond-prefilter [default: 24]

**\--max-species-divergence** INT

  Maximum number of different bases acids to allow between a sequence
    and the best hit in the database so that it is assigned to the
    species level. [default: 2]

**\--exclude-off-target-hits**

  Exclude hits that are not in the target domain of each SingleM
    package

**\--min-taxon-coverage** FLOAT

  Minimum coverage to report in a taxonomic profile. [default: 0.35
    for reads, 0.1 for genomes]

**\--working-directory** directory

  use intermediate working directory at a specified location, and do
    not delete it upon completion [default: not set, use a temporary
    directory]

**\--working-directory-dev-shm**

  use an intermediate results temporary working directory in /dev/shm
    rather than the default [default: the usual temporary working
    directory, currently /tmp]

**\--force**

  overwrite working directory if required [default: not set]

**\--filter-minimum-nucleotide** length

  Ignore reads aligning in less than this many positions to each
    nucleotide HMM [default: 72]

**\--include-inserts**

  print the entirety of the sequences in the OTU table, not just the
    aligned nucleotides [default: not set]

**\--known-otu-tables** *KNOWN_OTU_TABLES* [*KNOWN_OTU_TABLES* \...]

  OTU tables previously generated with trusted taxonomies for each
    sequence [default: unused]

**\--no-assign-taxonomy**

  Do not assign any taxonomy except for those already known [default:
    not set]

**\--known-sequence-taxonomy** FILE

  A 2-column \"sequence\<tab\>taxonomy\" file specifying some
    sequences that have known taxonomy [default: unused]

**\--no-diamond-prefilter**

  Do not parse sequence data through DIAMOND blastx using a database
    constructed from the set of singlem packages. Should be used with
    \--hmmsearch-package-assignment. NOTE: ignored for nucleotide
    packages [default: protein packages: use the prefilter, nucleotide
    packages: do not use the prefilter]

**\--diamond-prefilter-performance-parameters** *DIAMOND_PREFILTER_PERFORMANCE_PARAMETERS*

  Performance-type arguments to use when calling \'diamond blastx\'
    during the prefiltering. By default, SingleM should run in \<4GB of
    RAM except in very large (\>100Gbp) metagenomes. [default: use
    setting defined in metapackage when set, otherwise use
    \'\--block-size 0.5 \--target-indexed -c1\']

**\--hmmsearch-package-assignment**

  Assign each sequence to a SingleM package using HMMSEARCH, and a
    sequence may then be assigned to multiple packages. [default: not
    set]

**\--diamond-prefilter-db** *DIAMOND_PREFILTER_DB*

  Use this DB when running DIAMOND prefilter [default: use the one in
    the metapackage, or generate one from the SingleM packages]

**\--assignment-threads** *ASSIGNMENT_THREADS*

  Use this many processes in parallel while assigning taxonomy
    [default: 1]

**\--sleep-after-mkfifo** *SLEEP_AFTER_MKFIFO*

  Sleep for this many seconds after running os.mkfifo [default:
    None]

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

Get a taxonomic profile for dsDNA phages from paired read input:

  **\$ lyrebird pipe -1 \<fastq_or_fasta1\> -2 \<fastq_or_fasta2\> -p
    \<output.profile.tsv\>**

Get a taxonomic profile Krona diagram for dsDNA phages from single read input:

  **\$ lyrebird pipe -i \<fastq_or_fasta\> \--taxonomic-profile-krona
    \<output.profile.html\>**

Gather an OTU table (per marker sequence groupings) from paired reads:

  **\$ lyrebird pipe -1 \<fastq_or_fasta1\> -2 \<fastq_or_fasta2\>
    \--otu-table \<output.otu_table.tsv\>**
