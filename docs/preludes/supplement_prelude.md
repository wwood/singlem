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

## weebill database
When the input metapackage bundles a weebill database (as the default metapackage does), `supplement` adds the new genomes to it with `weebill db-add`, so that `singlem pipe` on reads can profile them with weebill too. This requires weebill >= 0.4.0. Only genomes which contribute at least one marker gene sequence to the new metapackage are added, so that every genome weebill reports can be mapped to a taxonomy. Sketching takes time roughly in proportion to the total size of the new genomes.

To create a metapackage without a weebill database instead, use `--no-weebill`. Then `singlem pipe` on reads with the output metapackage profiles from the marker genes alone.

When the input metapackage does not bundle a weebill database, one can be added with `--weebill-db-base <db.syl2db> --weebill-c <c>`, where the database contains the genomes of the input metapackage, and was built at the given `-c`. For instance, the `gtdb.syl2db` (built at `-c 100`) of the default metapackage from the same GTDB release can be used.

Supplementing a genome whose FASTA file name is already in the weebill database is an error. If a new genome is in the same species as a genome already in the metapackage, weebill splits that species' coverage between the two genomes. Since both have the same taxonomy, this does not affect the taxonomic profile.
