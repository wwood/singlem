It performs three steps:

1. Finding discrete operational taxonomic units (OTUs) from shotgun metagenome data
2. Assign taxonomy to OTUs
3. Converting OTU information into a coherent taxonomic profile

In the first step, SingleM finds OTUs without heavy reliance on reference sequence databases. These OTUs are often able to differentiate species even if they are from lineages new to science. To do this, SingleM finds reads that encode conserved single copy marker genes. It specifically finds reads which cover short (20 amino acid / 60 base pair) highly conserved sections ("windows") within those genes. Sequences covering those small sections are OTU sequences that exist independent of taxonomy. By default, SingleM uses 35 bacterial and 37 archaeal single copy marker genes.

In the second step, taxonomy is assigned based on comparing the nucleotide sequence of the window to GTDB species representatives' window sequences. If none are similar enough (i.e. within 95% or 3bp of the 60), then diamond blastx is used.

A common analysis is to convert a set window sequences (i.e. a metagenome's OTU table) into a taxonomic community profile. This is achieved in the second step by considering the OTUs from the 59 different marker gene holistically, using trimmed means and expectation maximisation.
