# Architecture

> **Note:** Despite the request that generated this file mentioning "Rust", SingleM is
> a **Python** package. The only `Cargo.toml` files in the tree are vendored dependencies
> inside a Python virtualenv under `extras/`. This document describes the Python source in
> `singlem/`.

## Overview

SingleM profiles the taxonomic composition of shotgun metagenomes (short- and long-read)
without reference genomes for every organism, by focusing on short, conserved "windows"
within single-copy marker genes. Reads are searched against marker-gene HMMs/DIAMOND
databases, the ~20-60 bp window is extracted and translated as needed, identical windows
are collapsed into **OTUs**, and OTU windows are matched against a reference sequence
database to assign taxonomy. A separate "condense" step combines per-marker OTU signal
into a single coverage-weighted **taxonomic profile**. The design goal is robustness to
novel lineages: because it operates on short conserved windows rather than whole genomes,
it detects and quantifies organisms absent from reference databases.

## Module Map

Ordered roughly by where a reader should start.

- **`main.py`** — CLI argument parsing and dispatch (via the `bird_argparser` helper). Each
  subcommand (`pipe`, `renew`, `summarise`, `appraise`, `condense`, `makedb`, `query`,
  `data`, `supplement`, `metapackage`, `create`, `chainsaw`, `prokaryotic_fraction`, …)
  is validated here and routed to the corresponding module. Start here to trace any command.
- **`pipe.py`** — The core workflow. `SearchPipe.run()` / `run_to_otu_table()` orchestrate
  search → read extraction → window finding → taxonomy assignment → OTU table, and
  optionally chain into `Condenser` for a taxonomic profile. Also defines the intermediate
  result holders `SingleMPipeSearchResult`, `SingleMPipeSeparateSearchResult`,
  `SingleMPipeTaxonomicAssignmentResult` (thin wrappers over on-disk GraftM/DIAMOND output dirs).
- **`diamond_spkg_searcher.py`** — `DiamondSpkgSearcher` runs the DIAMOND "prefilter" that
  cheaply narrows reads to likely marker hits before per-package processing. `_search` in
  `pipe.py` alternatively drives HMM search via GraftM.
- **`pipe_sequence_extractor.py`** — Given search hits, extracts the reads/ORFs relevant to
  each marker package and sample. `PipeSequenceExtractor` produces `ExtractedReads` /
  `ExtractedReadSet`; helper functions align proteins to HMMs and filter via hmmsearch.
- **`metagenome_otu_finder.py`** — `MetagenomeOtuFinder.find_windowed_sequences()` /
  `find_best_window()` implement the central idea: choosing and extracting the fixed-width
  window columns from an alignment (protein or nucleotide).
- **`frameshift_repair.py`** — Pure helpers used by `pipe` frameshift repair (on by default, disable
  with `--no-repair-frameshifts`). `walk_btop()`
  turns DIAMOND's BTOP string from the prefilter into the positions of single-base indels in
  the read; `repair_frameshifts()` edits the read to restore its reading frame (inserting an
  `N` for a deleted base, dropping an inserted one); `resolve_ambiguous_windows()` then fills
  each `N` from the most abundant near-identical window of the same marker. Matters because an
  indel, unlike a substitution, breaks the translated alignment entirely — the dominant cause
  of lost windows on Nanopore reads.
- **`prefilter_pad.py`** — `PrefilterPadder` plus the pure helpers `window_alignment_positions()`
  and `pad_aligned_sequence()`. Emits a prefilter FASTA where every on-target sequence is
  padded to a fixed length (30aa + 20aa window + 30aa, `X`-padded) with the window in a
  constant position, so DIAMOND can reject alignments that miss the window. Reuses
  `MetagenomeOtuFinder`'s column logic so the window matches `pipe`.
- **`sequence_classes.py`** — Low-level sequence types: `Sequence`, `AlignedProteinSequence`,
  `UnalignedAlignedNucleotideSequence`, and `SeqReader` (a fast FASTA/FASTQ `readfq` generator).
- **`otu_table.py`, `otu_table_entry.py`, `archive_otu_table.py`, `otu_table_collection.py`** —
  The OTU data model. `OtuTableEntry` is one (marker, sample, window sequence, count,
  coverage, taxonomy) row; `OtuTable` is the plain TSV form; `ArchiveOtuTable` is the richer
  JSON form (retains per-read data, enabling `renew`); collection classes provide streaming
  iteration over many tables (`StreamingOtuTableCollection`). Archive version 5 writes the
  OTUs column-wise, hoists fields that are constant across every OTU into `constant_fields`,
  and dereplicates read sequences into a shared `reads` list; versions 1-4 are still read.
  In memory the table is a list of rows either way, so only `ArchiveOtuTable` itself cares.
- **`condense.py`** — `Condenser` turns per-marker OTUs into a `CondensedCommunityProfile`
  (a coverage tree of `WordNode`s), applying trimmed means and genus/species
  expectation-maximisation to resolve coverage across taxonomic levels. `...KronaWriter`
  renders Krona HTML. This is the "taxonomic profile" producer.
- **`metapackage.py`** — `Metapackage`: the bundle of marker packages plus the DIAMOND
  prefilter DB, taxonomy databases, and genome-size data. Handles `acquire`/`download`/
  `generate`/`verify` and exposes per-marker HMM/DIAMOND paths used by `pipe`.
- **`singlem_package.py`** — `SingleMPackage` (a.k.a. "spkg"): one marker gene's HMM, GraftM
  package, window position/size, target domains, and taxonomy hash. Versioned subclasses
  `SingleMPackageVersion1..4` add fields over time; always accessed via `SingleMPackage.acquire`.
- **`sequence_database.py`** — `SequenceDatabase` (the "sdb"): the reference window database
  used by `query`/`appraise`, with pluggable nearest-neighbour indexes (smafa-naive, annoy,
  nmslib, scann) over nucleotide/protein binary encodings. `makedb` builds it.
- **`querier.py`** — `Querier` searches query windows against an sdb by sequence similarity
  using the selected backend; defines `QueryInputSequence`, `QueryResult`.
- **`appraiser.py`** — `Appraiser` compares a metagenome OTU table against genome/assembly
  OTU tables to estimate what fraction of the community is represented (`AppraisalResult`,
  `AppraisalBuildingBlock`); supports a streaming mode.
- **`summariser.py`** — `Summariser`: format conversions and reports over OTU tables and
  taxonomic profiles (wide/long, Krona, UniFrac, rarefaction, translation, species-by-site).
- **`condense`/profile consumers:** `read_fraction.py` (`prokaryotic_fraction`/
  `microbial_fraction`), `strain_summariser.py`, `genome_size.py`.
- **`renew.py`** — Re-runs taxonomy assignment on an existing archive OTU table (reusing
  stored reads) against a new metapackage, skipping the expensive search step.
- **`supplement.py`** — Adds new genomes to a metapackage (builds new spkg entries / taxonomy).
- **`lyrebird.py`** — dsDNA-phage variant of the workflow (uses a Lyrebird metapackage).
- **`package_creator.py`, `regenerator.py`, `chainsaw.py`, `trim_package_hmms.py`** —
  spkg/metapackage construction and maintenance tooling.
- **`taxonomy.py`, `taxonomy_bihash.py`** — Taxonomy string parsing/splitting utilities.
- **`utils.py`, `run_via_os_system.py`, `ordered_set.py`, `biolib_lite/`** — Misc helpers and
  a vendored slice of biolib.

## Key Types & Data Flow

Central types:

- **`SingleMPackage` (spkg)** — one marker gene's search + placement resources and the window
  definition. Many spkgs form a **`Metapackage`**, which also holds the shared prefilter DB
  and taxonomy/genome-size data. Every stage consults these for HMMs, window positions, and taxonomy.
- **`OtuTableEntry` / `OtuTable` / `ArchiveOtuTable`** — the pivot data structure: a window
  sequence with its marker, sample, count, coverage and taxonomy. Everything downstream
  (query, appraise, summarise, condense) consumes OTU tables.
- **`CondensedCommunityProfile`** — the taxonomic profile: a coverage tree keyed by taxonomy,
  produced from OTU tables by `Condenser`.
- **`SequenceDatabase` (sdb)** — the searchable reference of window sequences behind
  `query`/`appraise`, distinct from the metapackage.

**Profiling (`pipe`), end to end:** input reads → `SearchPipe._search` (HMM via GraftM) or
`DiamondSpkgSearcher` (DIAMOND prefilter) narrow reads to marker hits → `PipeSequenceExtractor`
pulls the relevant reads/ORFs per spkg/sample → `MetagenomeOtuFinder` aligns them and extracts
the fixed window columns → identical windows are collapsed into `OtuTableEntry`s with counts
and coverage → `_assign_taxonomy` (smafa-naive then DIAMOND, or query/pplacer) labels each
window → an `OtuTable`/`ArchiveOtuTable` is written. If a taxonomic profile is requested, the
OTU table is streamed into `Condenser.condense`, which builds a `CondensedCommunityProfile`
and writes the profile and/or Krona.

**Database build/query:** `makedb` (`SequenceDatabase.create_from_otu_table`) encodes window
sequences into binary and builds nearest-neighbour indexes; `query` (`Querier`) encodes query
windows the same way and searches those indexes by divergence.

**Renew:** an `ArchiveOtuTable` retains per-read sequences, so `renew.py` re-assigns taxonomy
against a newer metapackage without re-searching the raw reads.

## Design Decisions & Invariants

- **Windows, not whole genes/genomes.** The fixed-width conserved window (defined per spkg in
  `singlem_position`/`window_size`) is the load-bearing idea — it makes identical-sequence
  collapse meaningful and gives sensitivity to novel lineages. `MetagenomeOtuFinder` is where
  this lives.
- **Versioned packages with subclassing.** `SingleMPackageVersion1..4` (and Metapackage
  `version` checks) exist so old databases keep working while new fields (window size, target
  domains, taxonomy hash, viral support) are added. Several features guard on
  `metapackage.version` (e.g. taxonomic profile requires v3+, viral requires v6). Always go
  through `acquire`, never construct directly.
- **On-disk intermediates as the interface.** `SingleMPipeSearchResult` and siblings are thin
  wrappers computing paths inside GraftM/DIAMOND output directories rather than holding data in
  memory — this keeps `pipe` memory-bounded on large samples but means the directory layout is
  an implicit contract between search, extraction and assignment stages.
- **Streaming OTU collections.** Collection/condense code is written to stream tables to bound
  memory (e.g. `condense` only materialises a list when Krona output forces a second pass).
- **Pluggable similarity backends.** `SequenceDatabase`/`Querier` support smafa-naive, annoy,
  nmslib and scann; the binary encoding of nucleotides/proteins is shared across them.
  (Rationale inferred from structure — the default is smafa-naive.)
- **Archive vs plain OTU table.** The archive (JSON) form must retain per-read data for `renew`
  to work; the plain TSV form is lossy. Don't assume they're interchangeable.
- OTU `data`/`fields` list mutation: `OtuTableEntry.add_found_data` deliberately copies `fields`
  before appending to avoid mutating a shared class-level default — a subtle invariant worth
  preserving.

## Entry Points

- **Add/modify a subcommand or its args:** `main.py` (argparse setup + the `args.subparser_name`
  dispatch chain near the bottom).
- **Change how metagenomes are profiled:** `pipe.py` (`SearchPipe.run_to_otu_table`), then
  `pipe_sequence_extractor.py` and `metagenome_otu_finder.py`.
- **Change taxonomic-profile logic (coverage, EM, trimming):** `condense.py`.
- **Add a marker gene / rebuild reference data:** `singlem_package.py`, `package_creator.py`,
  `metapackage.py` (and `supplement.py` to add genomes).
- **Build or query the reference sequence DB:** `sequence_database.py` (build) and `querier.py`
  (query); add a new NN backend by extending both.
- **New output format/report:** `summariser.py`.
- **Estimate community representation / read fraction:** `appraiser.py`, `read_fraction.py`.
- **Re-assign taxonomy on existing results:** `renew.py`.
- **Phage profiling:** `lyrebird.py`.

## Out of Scope

Update this file when modules are added, removed, or restructured — not for routine
function-level changes. It documents the shape of the codebase, not implementation detail.
