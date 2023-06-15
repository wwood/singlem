The 'read_fraction' subcommand can be used to estimate the fraction of reads
from a metagenome that are assigned to Bacteria and Archaea compared to e.g.
phages or eukaryotes.

The advantage of this method is that it does not require a reference sequence  of
non-microbial genomes present in a metagenome (e.g. those of a animal host).
Instead, it uses a SingleM taxonomic profile of the metagenome to "add up" the
components of the community which are microbial, and assumes the remainder is
non-microbial (e.g. host or phage).
