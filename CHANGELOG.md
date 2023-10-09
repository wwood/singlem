## v0.15.0
* Genomes that encode proteins with translation table 4 are now supported. This
  works by assuming all genomes have translation table 4, since regular sequence
  similarity search excludes inappropriately translated sequences from genomes
  which use table 11 (the standard bacterial table). Thanks to Dr. Andy Leu for
  useful test cases. NOTE: The `renew` mode is not sufficient for detecting
  these lineages, `pipe` must be run again from scratch.
* new_package_creation (beta): A snakemake pipeline included in the `extras`
  directory used to create new SingleM metapackages from scratch. In
  development. Thanks for @harmonydouwes @tvtv195 @JemmaSun for testing.
* Version S3.2.1 of the default metapackage released, which includes updated
  genome sizes for GTDB genomes (for use with `read_fraction`), now corrected
  based on CheckM v2 estimates of completeness and contamination. Thanks to
  @EisenRa for collaboration.
* `seqs`: Output the best window position to STDOUT.
* Other assorted bug fixes and documentation updates.

## v0.14.0
This release is a huge step forward for the SingleM software, comprising >750 git commits and several years work (particularly from @AroneyS and @EisenRa and @rzhao-2) since v0.13.2. 

There are so many changes that generating a CHANGELOG would take too long.

This release is equivalent to 1.0.0beta8, and is intended as a pre-release for version 1.0.0, but using a standard version number allows for a more streamlined release process.

## v0.10.0 to v0.13.2
(undocumented)

 ## v0.9.0
* Appraise can now generate 'appraisal plots'
* Use of smafa / SQLite rather than BLAST+ / VSEARCH for 'query' and clustering. SingleM
  databases (.sdb) will need to be regenerated.
* SingleM databases can now be queried via taxonomy or sample name
* Overhaul of command line help messages
* Appraise can now appraise assemblies as well as genomes
* Various bug fixes and enhancements

## v0.8.2
* Fix for installation through PyPI.

## v0.8.1
* Fix bug in singlem query where some results were omitted.
* Detect when max_target_seqs has been reached in singlem query.

## v0.8.0
* Overhauled makedb/query. Database creation is now faster and querying more accurate, especially for OTUs with gaps. Old databases should be re-generated.
* summarise: Added BIOM and wide format outputs - props to Steve Robbins and Louis Monaghan for the suggestions.
* appraise: Default to genus-level similarity cutoff.

## v0.7.1
* Now installable via pip / PyPI.

## v0.7.0
* Speed improvements for singlem query
* Memory improvements for singlem makedb
* db: Use and require the new diamond version / database format
* summarise: Only output to a single html for --krona where possible

