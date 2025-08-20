The SingleM `supplement` subcommand adds genomes to a SingleM metapackage.

**TLDR**: Supplement a metapackage with new genomes like so:
```
singlem supplement --new-genome-fasta-files <genome1.fna> <genome2.fna> \
    --output-metapackage <supplemented.smpkg>
```

## GTDB-Tk
In order to add genomes to a metapackage, their taxonomy is required. SingleM `supplement` can generate this taxonomy for new genomes using GTDB-Tk. However, GTDB-Tk is not installed by default, and so it must be installed separately.

For instance, if you are using conda to manage dependencies, use:
```
conda install gtdbtk
```
Note that the version of GTDB-Tk installed must match the GTDB release of the metapackage. The [GTDB-Tk documentation](https://ecogenomics.github.io/GTDBTk/installing/index.html) provides guidance on installation and compatibility.
