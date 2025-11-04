## v0.20.0

Major new function - Long-read input support (Nanopore >= R10.4.1 or PacBio HiFi recommended), thanks to @thepatientwait.

* Lyrebird database updated to v0.3.1, improving exclusion of off-target (non-phage) sequences
* `microbial_fraction` subcommand renamed to `prokaryotic_fraction` (old name retained as synonym)
* More flexible options for specifying genome input in `pipe` mode
* `appriase` mode: Add `--stream-inputs`
* [GlobDB R226 metapackage](https://fileshare.lisc.univie.ac.at/globdb/globdb_r226/taxonomic_profiling/globdb_r226_SingleM_metapackage.tar.gz) released

Thanks to @AroneyS, @rzhao-2, @EisenRa, @thepatientwait, @dspeth, @Anna-MarieSeelen, @luigallucci, @ilnamkang and others for contributions and testing.

## v0.19.0
Major new function - profiling of Caudoviricetes (aka "Caudovirales") phage communities (Lyrebird), thanks to @rzhao-2.

![Lyrebird](https://raw.githubusercontent.com/wwood/singlem/refs/tags/v0.19.0/docs/_include/lyrebird_resized.png?raw=1)

Other changes:
* Update default metapackage to GTDB R226
* admin: Use pixi instead of conda
* Use of diamond v2.1.10 specifically, to avoid segfault issues with diamond v2.1.11
* Clarify non-standard metapackage usage (#220)
* doc: Improve summarise --cluster (#210)

Thanks @rzhao-2 @AroneyS @ilnamkang Phil Hugenholtz @pchaumeil @zackhenny @thepatientwait

## v0.18.3
A small patch release

* `summarise`: Fix a regression
* docs: Minor fixes

## v0.18.1
A small patch release.

* Updates to singlem `supplement` and other modes for polars >1.0
* Pin dependencies to help future proof singlem

## v0.18.0
Combined changelog for v0.17.0 and 0.18.0

* Use of GTDB R220 reference metapackage by default
* `pipe`/`condense`: Improve algorithm by delaying some filtering steps, leading to more accurate taxonomic profiles
* `pipe`: update to [smafa](https://github.com/wwood/smafa) v0.8.0 for substantial speed improvement
* `microbial_fraction`: Remove `%` from column data and add average genome size estimation
* `supplement`: Change command line options in backwards incompatible way, clarifying their meaning
* `summarise`: Add `--output-taxonomic-profile-with-extras` output to add relative abundance etc. to taxonomic profiles
* `summarise`: Add `--output-species-by-site-relative-abundance-prefix` to create taxon-level specific relative abundances from taxonomic profiles
* `summarise`: Add `--output-taxonomic-level-coverage` to show how much coverage and number of taxa assigned to each level
* `pipe`: Faster processing when many genome fasta files are input
* `seqs`: Prioritise high-info HMM positions.
* dist: Fix singularity container
* assorted bug and documentation fixes

Thanks @AroneyS @EisenRa @jakobnissen @rzhao-2 @rrohwer @shaze @ellyyuyang @VadimDu @adityabandla @luispedro, and anonymous reviewers, among others.

The `microbial_fraction` mode now has its own citation - https://www.biorxiv.org/content/10.1101/2024.05.16.594470v1

## v0.16.0
This version tweaks the method which assign taxonomy to OTUs (increasing the species-level threshold) and the method which summarises the OTUs to create a final taxonomic profile (very low abundance lineages are given lower taxonomic resolution, rather than ignored completely). This improves the rate over "overclassification" i.e. when novel species are classified wrongly to the species level, and improves the `read_fraction` (now called `microbial_fraction`) estimates in complex / shallowly sequenced metagenomes.

We suggest recomputing community profiles using `renew` or `pipe` modes.

* pipe/renew: Change default species-level assignment from 3bp or closer, to 2bp or closer.
* pipe/renew/condense: Assign sub-min-taxon-coverage higher.
* read_fraction mode renamed to microbial_fraction

Thanks to Yu Yang, Caitlin Singleton, @MadsAlbertsen @EisenRa @BigDataBiology

## v0.15.1
Mostly minor bugfixes

* pipe: extract: Apply --evalue to hmmsearch thresholding.
* Fix for appraise --plot
* pipe: Dedup hmmsearch results during diamond package assignment.
* pipe/renew/condense: Prevent no_assign_taxonomy and taxonomic profile output.

Thanks @kalonji08 @AroneyS @harmonydouwes

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

