The 'read_fraction' subcommand can be used to estimate the fraction of reads
from a metagenome that are assigned to Bacteria and Archaea compared to e.g.
phages or eukaryotes.

The advantage of this method is that it does not require a reference sequence of
non-microbial genomes present in a metagenome (e.g. those of an animal host).
Instead, it uses a SingleM taxonomic profile of the metagenome to "add up" the  
components of the community which are microbial, and assumes the remainder is
non-microbial (e.g. host, diet, or phage).

This tool can help prioritise samples for deeper seqeuncing, forecast how much sequencing
is required for a given set of samples, and identify problematic steps in genome-resolved
metagenomic pipelines. To learn more about this tool (including possible use cases), check
out this preprint {LINK TO BE INSERTED}.

A warning will be triggered when a sample's microbial_fraction estimate could under- or 
overestimate the read fraction by >10%. Additionally, if the 3 highest abundance lineages 
not classified to the species level would change the estimated read fraction of the sample 
by >2% if their genome size is halved or doubled. Microbial fractions are also capped at
100%, but the original value can be recovered from the output:
`bacterial_archaeal_bases / metagenome_size`

Here is an example output from 'microbial_fraction':

| sample      | bacterial_archaeal_bases | metagenome_size | read_fraction | warning |
| ----------- | ------------------------ | --------------- | ------------- | ------- |
| KRGM_94_05gb_M_1 | 473147661           | 512995675       | 92.23%
