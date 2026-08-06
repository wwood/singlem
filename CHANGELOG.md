## Unreleased

* Archive OTU tables are now written at version 5, which stores the OTUs column-wise, stores fields that are constant across every OTU once rather than per OTU, and dereplicates read sequences into a shared list. Versions 1-4 are still read, and `ArchiveOtuTable` presents the same list of rows whichever version it read, but code that parses the JSON itself rather than going through `ArchiveOtuTable` needs updating - see the Archive OTU table entry in the glossary for the layout.
* Archive OTU tables version 5 also add a `reads_with_repaired_deletions` field, recording which of an OTU's reads have an ambiguous base that frameshift repair inserted.
* `pipe`: When processing one chunk of a dataset with `--read-chunk-size`/`--read-chunk-number`, the ambiguous bases left by frameshift repair are no longer resolved within the chunk. Which window a read ended up with otherwise depended on what else was in that chunk, so the chunks did not combine into the result the whole dataset would have given.
* `pipe`: Fix the `nucleotides_aligned` of an archive OTU table being ordered independently of `read_names` and `read_unaligned_sequences`, which were sorted by read name while it was not. The two orders only differ when the reads of one OTU do not all have the same aligned length, which takes an insert column inside the window's span that is a gap in some of them, so this was rarely visible - but `renew` reads the three as parallel.
* `summarise`: Add `--resolve-ambiguous-windows` (and `--max-frameshift-repair-divergence`), which resolves those bases while archive OTU tables are combined, where every chunk is present. Use it when combining the chunks of a dataset processed with `pipe --read-chunk-size`.
* `summarise`: `--collapse-to-sample-name` and `--collapse-paired-with-unpaired-archive-otu-table` now combine archive OTU tables of different versions, rather than rejecting them. The output takes the newest version's fields, and OTUs from an older archive have no value for the fields that version lacked. Previously an archive made before this release could not be combined with one made after it.
* `pipe`: Fix `--assignment-method annoy` failing outright on single-ended input with `UnboundLocalError: cannot access local variable 'read_name_to_fullseq'`. Pre-existing, and missed because the tests covering it are skipped unless `annoy` is installed, which it is not in CI.

* `pipe`: Frameshift repair, which uses the frameshifts DIAMOND reports during the prefilter to restore the reading frame of each read before it is aligned to the HMM, so reads carrying single base indels still yield a window, is now on by default. Substantially improves window recovery on Nanopore data, where indels rather than substitutions are the dominant error. Where a repaired deletion leaves a base of unknown identity, it is taken from the most abundant window within `--max-frameshift-repair-divergence` mismatches, and only for windows where repair actually inserted that base, not a pre-existing ambiguous base in the raw read. Disable with `--no-repair-frameshifts`. 

## v0.21.3

* Fix PyPI dependency list by generating `admin/requirements.txt` in CI before the wheel build; previous releases shipped without dependencies because the generated file was never committed to release tags.

## v0.21.2

* Fix PyPI dependency list structurally, and include bird_tool_utils.

## v0.21.1

* Fix PyPI dependency list

## v0.21.0

* Default metapackage updated to GTDB R232 (S6.5.0-R232)
* `pipe`: Substantially improved logging, including progress bars for prefilter, hmmsearch and DIAMOND assignment steps. Thanks @thepatientwait.
* `pipe`: Parallelised DIAMOND taxonomic assignment for faster runs. Thanks @thepatientwait.
* `pipe`: Add optional `--context-window` if the entire query sequence is too large for output size.
* `pipe`: Add zstd input support in places, including with `--read-chunk-size`/`--read-chunk-num` and `prokaryotic_fraction`
* `pipe`: Restore v0.19.0 orfm+mux pipeline for genome input, restoring performance lost in v0.20
* `pipe`: Use DIAMOND `--frameshift` when assigning taxonomy
* `pipe`: `--sra-files`: Numerous fixes and refinements to chunked SRA extraction, including detection of `kingfisher` failures and stripping of `.sra` suffix from sample names
* `pipe`/`renew`: Bugfix for missing sequences with Lyrebird
* `supplement`: Add `--output-matched-faa` to output matched protein sequences, with marker name included in sequence IDs
* `renew`: Accept zipped gzip archive OTU tables as input
* `lyrebird`: Add `extras/lyrebird_metapackage_creation` snakemake workflow for building Lyrebird metapackages from scratch
* `data`: Report DOI when database download or acquisition fails
* Switch internal data handling from pandas to polars
* Updated citation for `prokaryotic_fraction`
* docs: Improved GlobDB documentation (#267)
* Add `AGENTS.md` and `SKILL.md` to assist agentic coding tools
* Assorted bug fixes, test improvements, and documentation updates

Thanks @thepatientwait, @rzhao-2, @AroneyS, @EisenRa and others for contributions and testing, and thanks to @MCeciC (#269) and @magicprotoss (#267) for raising issues.

## v0.20.3

Minor bugfix release. Rare sequences tripped a bug in `pipe` mode.

Thanks to @asuq for reporting #265.

## v0.20.0 / 0.20.1 / 0.20.2

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

