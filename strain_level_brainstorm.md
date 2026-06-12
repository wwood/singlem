# Strain-level information in the community profile, with reference-free compression

## The goal, restated

We want to (a) extend the condensed community profile *below species* to strain
(or sub-species) resolution, and (b) do so under a "compress once, assign many"
discipline:

1. **Compress** the `pipe` read input into a sample representation **without
   reference to any genome database** — the compression must depend only on the
   reads (and fixed, DB-independent definitions such as the marker HMMs or a hash
   function), not on which genomes we will later compare against.
2. **Assign** strain-level taxonomy *later*, from the compressed representation
   **alone** (no original reads), against **different reference genome databases**
   — swapping the database (a newer GTDB release, a clade-specific strain panel,
   a private MAG collection) re-profiles the *same* stored sketch.

This is the "sketch once, profile later" paradigm. SingleM and sylph already each
realise half of it; the proposals below push both to strain level and unify them,
and explain how the result plugs into the same clade-budget reconciliation that
the `condense` sylph-integration uses.

## What already separates compression from assignment

Two reference-free compressions of the reads are already produced today, and for
each the database only enters at assignment time:

- **Marker-window haplotypes (the archive OTU table).** `pipe` recovers short
  (~60 bp) protein-coding windows of single-copy marker genes from the reads. The
  windows are defined by the marker HMMs — *not* by any genome database — so the
  archive OTU table is a reference-free compression of the reads. It stores, per
  window, the exact nucleotide `sequence`, its count/coverage, and the
  contributing `read_unaligned_sequences`; taxonomy is a *separate* field.
  `renew` already re-derives taxonomy from an archive table against a *different*
  metapackage, with no access to the reads — exactly property (2), at
  species/genus level.

- **FracMinHash genome sketches (sylph `.sylsp`).** A bottom-hash subsample of
  every k-mer in the reads, defined by a hash function and a subsampling rate
  `-c` — again reference-free. `pipe --output-sylph-sketch` / `renew
  --input-sylph-sketch` already let a sketch be profiled later, and profiling the
  same sketch against a different `.syldb` re-assigns coverages without the reads.

The strain question is therefore not "how do we separate compression from
assignment" — that architecture exists — but "what do we put in the compressed
representation so that *strain* signal survives, and how do we assign and report
it."

## Where strain signal lives in each compression

- A marker window is ~60 bp of coding sequence; across ~14–60 single-copy markers
  that is ~1–3 kb of concatenated, strain-discriminating sequence per genome.
  Two strains of one species differing by SNPs in these windows appear as
  **distinct OTUs with distinct `sequence` values** — the signal is already
  captured losslessly; it is simply collapsed to species during condense.
- A FracMinHash sketch retains, at rate `1/c`, the k-mers that distinguish
  strains; sylph's containment-ANI can separate genomes down to a few percent
  divergence when the database contains the relevant strains and `c` is small
  enough.

## Proposals

### A. Marker-window haplotype strains (reference-free compression; swappable allele DB)

Keep the archive OTU table as the compression (already reference-free). Add a
strain-assignment step, analogous to `renew`'s taxonomy reassignment, that matches
each window's nucleotide `sequence` against a **strain-resolved reference of
marker alleles** — i.e. the per-genome window sequences of every genome in the
chosen database. Different databases (different genome sets, GTDB releases, or a
focused strain panel) yield different strain assignments from the *same* stored
windows.

- **Reported as:** sub-species nodes (an 8th `t__`/strain rank) in the condensed
  tree. `WordNode` already supports arbitrary depth, and the writers/Krona/
  `taxonomic_level_coverage_table` need only a "strain" level added.
- **De novo (fully reference-free) strains:** even without any allele database,
  windows of the *same* species that co-vary in abundance across markers and
  across samples can be clustered into strain haplotypes (a SingleM analogue of
  StrainPhlAn / DESMAN). This yields anonymous strains ("species X strain 1/2")
  with relative abundances, assignable to named strains later if a database is
  supplied.
- **Cost / limit:** tiny (kilobytes; already stored). Resolution is coarse —
  ~1–3 kb of conserved coding sequence distinguishes well-separated strains and
  tracks them across samples, but not fine microdiversity. Strains identical
  across all marker windows are indistinguishable by this route (hand them to B).

### B. Genome-sketch strains (reference-free `.sylsp`; swappable strain-panel DB)

Profile the stored FracMinHash sketch against a database that contains **multiple
genomes per species** (a strain panel). Sylph returns per-genome effective
coverage and ANI, i.e. strain-level coverages, and swapping the `.syldb`
re-profiles the same sketch — property (2) at genome resolution.

- **Improve resolution by sketching denser:** store the sketch at a smaller `-c`
  (more retained k-mers) so within-species k-mer differences survive; `c` is the
  one knob that trades sketch size against strain resolution, and it is fixed at
  compression time independent of the database.
- **Reported as:** the per-strain coverages become the species-subtree partition
  in condense (below).
- **Cost / limit:** sketch size grows as `1/c`; strain separation still requires
  the strains to be *in* the chosen database and to differ by more than sylph's
  ANI resolution. Blind to strains with no near reference (hand those to A's
  de novo route).

### C. A unified, reference-free "sample sketch" artifact

Define one container emitted by `pipe`, depending only on the reads:

- marker-window haplotypes + counts (the archive OTU table, sans the taxonomy
  field — or with it marked provisional),
- a FracMinHash sketch of the reads (`.sylsp`),
- read-count / base-count normalisation metadata (for coverage scaling and the
  sylph↔SingleM `alpha`),
- *optionally* per-marker SNV/allele-frequency tables (proposal D).

Assignment becomes a separate command — essentially today's `renew` + sylph
profile — taking `{artifact, one or more genome databases}` and producing a
strain-resolved community profile, with **no reads required and databases freely
swappable**. Most of the plumbing exists: `renew` already consumes the archive
table and (with `--input-sylph-sketch`) the sketch; this proposal is to package
them as one durable artifact and to formalise strain output. It also makes the
profile reproducible and GTDB-release-agnostic: the same artifact can be assigned
against r220 today and r226 next year.

### D. SNV / allele-frequency fingerprint (reference-free; for tracking + panels)

Beyond whole windows, store a compact per-sample vector of allele frequencies at
the variable positions of the marker windows (or at minimizer-anchored conserved
loci). This is a small, reference-free fingerprint that (i) tracks the *same*
strain across samples by fingerprint identity, and (ii) can be matched to
strain reference panels at assignment time. Smaller and more strain-discriminating
than raw windows because it keeps only the variable sites, at the cost of needing
a fixed locus definition.

### E. Denser reference-free sketch for fine strains (minimizers / marker-neighbourhood unitigs)

When markers are too coarse and a genome panel is unavailable, a denser
reference-free representation gives finer strains: a sourmash-style minimizer
sketch of all reads, or local assembly of unitigs in the neighbourhood of marker
loci. Either is computed from reads alone and aligned/contained against arbitrary
genome databases later. Larger than A–D but still far smaller than the reads, and
DB-swappable.

## How strain coverage enters the condensed profile

The clade-budget reconciliation already used for sylph integration extends one
rank deeper without new machinery. At the species→strain layer:

- **SingleM's species coverage is the budget.** The total coverage condense
  assigns to a species (leaf + any unresolved-within-species coverage) caps the
  strains beneath it, just as a genus capped its species.
- **Strain evidence sets the partition.** Window-haplotype strains (A/D) or
  genome-sketch strains (B) partition that budget; a residual "novel strain"
  node carries the within-species coverage no strain resolves (microdiversity /
  dark strains).
- **One-sided, never-reduce update.** The same signed draw-down rule implemented
  for `--sylph-reconcile` applies: a strain is lifted to its sketch/haplotype
  coverage by drawing from the species' unresolved budget, never reduced, and
  only coverage in excess of the species budget is added as new — so strain
  resolution is added without inflating or distorting the species totals SingleM
  is trusted for.

Because the budget/partition logic is rank-agnostic, the existing
`_apply_sylph_coverage` engine (and its order-independence proof) generalises to a
strain layer with little more than a deeper taxon array and a strain-rank-aware
writer.

## Recommended path

1. **B + C first.** Profiling the already-stored `.sylsp` against strain-panel
   databases is the lowest-effort strain layer and is the purest realisation of
   "sketch once, assign against different databases later." Package the archive
   table and sketch as the proposal-C artifact and add a strain rank to the
   condense tree and writers.
2. **A (de novo + allele DB) next**, for strains without a near reference and to
   exploit the strain signal already sitting unused in the window sequences.
3. **D/E** only where finer resolution is needed than A/B provide, accepting the
   larger sketch.

In all cases the database enters only at assignment, the compression is fixed by
the marker HMMs and the hash function, and the reconciliation that merges the
strain evidence into the community profile is the rank-generalised form of the
method already in `condense`.
