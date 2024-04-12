The SingleM `microbial_fraction` ('SMF') subcommand estimates the fraction of reads in a metagenome that are microbial, compared to everything else e.g. eukaryote- or phage-derived. Here we define 'microbial' as either bacterial or archaeal.

The main conceptual advantage of this method over other tools is that it does not require reference sequences of the non-microbial genomes that may be present (e.g. those of an animal host). Instead, it uses a SingleM taxonomic profile of the metagenome to "add up" the components of the community which are microbial. The remaining components are non-microbial e.g. host, diet, or phage. 

Roughly, the number of microbial bases is estimated by summing the genome sizes of each species/taxon after multiplying by their coverage in the taxonomic profile. The microbial fraction is then calculated as the ratio of microbial bases to the total metagenome size.

This tool can help prioritise samples for deeper sequencing, forecast how much sequencing is required for a given set of samples, and identify problematic steps in genome-resolved metagenomic pipelines.

The method is least reliable in simple communities consisting of a small number of species that are missing from the reference database. The main challenge in these cases is that the genome sizes of novel species are hard to estimate accurately. Multiplying a coverage from the taxonomic profile against an uncertain genome length equals an uncertain number of bases assigned as microbial.

To detect such situations SingleM `microbial_fraction` emits a warning when a sample's microbial_fraction estimate could under- or overestimate the read fraction. Specifically, the warning is emitted when the 3 highest abundance lineages not classified to the species level would change the estimated read fraction of the sample by >2% if their genome size is halved or doubled. Microbial fractions are also capped at 100% since values greater than this are impossible, but the original value can be recovered from the output if you calculate `bacterial_archaeal_bases / metagenome_size`.

Here is an example output from 'microbial_fraction':

| sample      | bacterial_archaeal_bases | metagenome_size | read_fraction | warning |
| ----------- | ------------------------ | --------------- | ------------- | ------- |
| KRGM_94     | 473147661                | 512995675       | 92.23         |         |

The `read_fraction` column is a percentage, i.e. 0.3 is 0.3%, nt 30%. The `warning` column is empty if no warning was emitted.
