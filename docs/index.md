# SingleM

SingleM profiles bacterial and archaeal communities directly from short- or long-read shotgun metagenomes. It estimates relative abundance accurately while retaining sensitivity to lineages absent from reference genomes.

Install SingleM and Lyrebird with Bioconda:

```bash
conda create -c conda-forge -c bioconda --override-channels --name singlem singlem
conda activate singlem
singlem data
```

Generate a taxonomic profile from paired reads:

```bash
singlem pipe -1 sample.1.fastq.gz -2 sample.2.fastq.gz -p sample.profile.tsv --threads 4
```

Start with the [tutorial](tutorial.md), follow a [how-to guide](how-to/profile-metagenomes.md), read the [CLI reference](reference/singlem-pipe.md), or learn [how SingleM works](explanation/marker-window.md).

SingleM is developed by the [Woodcroft lab](https://research.qut.edu.au/cmr/team/ben-woodcroft/) and licensed under [GPL-3.0-or-later](https://www.gnu.org/licenses/gpl-3.0.html).
