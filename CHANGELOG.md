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

