# Plan: weebill databases in supplemented metapackages

Background: `singlem pipe` on reads profiles with weebill when the metapackage
bundles a weebill two-stage database (`.syl2db`). It then condenses jointly with
`DEFAULT_JOINT_CONDENSE_ARGUMENTS` (`singlem/condense.py`), which pins every
species weebill reported at its own coverage and **fixes every other database
species at zero**. That makes it essential that (a) every genome weebill can
report can be mapped to a metapackage taxonomy, and (b) every species in the
metapackage is in the weebill database.

Change #1 (implemented separately) makes (a) hold for supplemented genomes. They
are matched by genome file name when there is no versioned GCA_/GCF_ accession.
This plan covers:

- **#2:** `singlem supplement` extends the weebill database, so (b) holds after
  supplementing.
- **#3:** making the genome → taxonomy lookup fast.

---

## #2 supplement extends the weebill database

### Behaviour

| Input metapackage | Flags | Output |
|---|---|---|
| bundles a `.syl2db` | (none) | the database grown with the new genomes via `weebill db-add` |
| bundles a `.syl2db` | `--no-weebill` | no weebill database, with a warning that `pipe` will then need `--no-weebill` |
| no `.syl2db` | (none) | no weebill database, as now |
| no `.syl2db` | `--weebill-db-base X.syl2db --weebill-c C` | `X` grown with the new genomes (e.g. a GlobDB supplement onto S7.5.0's `gtdb.syl2db`) |

Also:

- Erroring is replaced by extending. Remove the exception at
  `singlem/supplement.py:927`, and replace `test_weebill_db_in_input_metapackage_croaks`
  with tests of the behaviours above.
- If the input bundles more than one weebill database, error for now. Extending
  one would need a rule for which database the new genomes go into. Supporting
  several can come later, e.g. adding a separate small `.syl2db` instead of
  growing one.
- `--weebill-db-base` together with an input that already bundles a database is
  an error.

### Which genomes go in

Only genomes that end up with sequences in the new metapackage. Otherwise weebill
can report a genome that change #1 cannot map, and pinning then zeroes it.

1. Start from `new_genome_fasta_files` after quality filtering, dereplication and
   the "has a taxonomy" filter (`singlem/supplement.py:243`).
2. Intersect with the genomes that contributed at least one sequence to any spkg.
   `gather_hmmsearch_results` already knows these: it writes `genome‡transcript`
   names to `matched_transcripts_fna`. Return the genome set from there (or have
   `generate_new_metapackage` collect it), rather than re-deriving it from the
   finished read-name store.
3. Log how many genomes were excluded for having no marker hits. They are also
   absent from the SingleM side, so leaving them out of weebill is consistent.

### Commands

Run in `generate_new_metapackage`, before `Metapackage.generate`, in
`working_directory`:

```
weebill sketch -l <genome_list> -c <c> --min-spacing <spacing> -t <threads> -o <wd>/supplement_genomes
weebill db-add -D <old.syl2db> -o <wd>/<name> -t <threads> <wd>/supplement_genomes.syldb
```

- Sketch with `-l` rather than `db-add -g <genomes...>`: 10^5 genome paths would
  exceed the command-line length limit.
- `--min-spacing` must match the base database. `db-add -g` would inherit it, but
  the sketch-then-add route does not, and `db-add` only checks `-c` and `-k`.
  Default to 30, weebill's default, which `extras/update_metapackage` builds with.
  If weebill can report a database's spacing (`weebill inspect`), read it from
  there instead, and fail if they differ.
- Write the genome list with the exact paths supplement used. Change #1 keys on
  `remove_extension(basename(path))`, which is how supplement names genomes in
  read names.
- Name the output `<old db basename without .syl2db>` (e.g. `gtdb.syl2db`), or
  add `--weebill-db-name`. The `.syl2db` basename must be unique within the
  metapackage.
- `db-add` errors on a genome whose file name is already in the database, which
  catches re-supplementing the same genome. Keep that as an error (do not pass
  `--skip-existing`), and include the hint in our error message.

Then pass `weebill_db=[grown]`, `weebill_c=[c]` to `Metapackage.generate`
(`singlem/supplement.py:643`). It already copies the database in and writes
`weebill_dbs` into a version 7 `CONTENTS.json`.

### Checks before the expensive work

- If a weebill database will be built, check `shutil.which('weebill')` and that
  `weebill db-add --help` succeeds (weebill >= 0.4.0), before gene calling and
  hmmsearch. Fail early with "install weebill >= 0.4.0 or pass --no-weebill".
- Warn about each new genome whose genome key (change #1) clashes with an
  existing key in the metapackage. Such a clash only happens when a supplemented
  file name carries the versioned accession of a genome already present.

### Interaction with taxonomy

Unchanged from supplement today: if a new genome is in the same species as an
existing one, weebill now splits that species' coverage across two genomes with
the same taxonomy. condense sums per taxonomy, so this is harmless. Worth a
sentence in the docs.

### Packaging and docs

- Root `pixi.toml` already has `weebill >= 0.4.0`.
  `extras/update_metapackage/pixi.toml` is still `>= 0.3.0`; bump it when that
  pipeline is next run.
- Update the supplement docs (`docs/` and `--help`) for the new behaviour and
  flags. Update `extras/update_metapackage/README.md`, which says supplement
  cannot extend the database.
- Once #2 lands, `extras/add_weebill_to_supplemented/` is only needed for
  metapackages supplemented before it.

### Tests

The test data has no real genomes, which makes tests need a small weebill
database:

- Build a tiny `.syl2db` in `test/data` from 1 or 2 small genomes. `weebill
  sketch` + `db-convert` runs in seconds on test-sized FASTAs. Reuse the genomes
  in `test/data/supplement/`.
- `supplement` on a metapackage bundling it gives an output whose `CONTENTS.json`
  has `weebill_dbs`. `weebill inspect` on the output lists the new genome file.
- `--no-weebill` gives an output without `weebill_dbs`.
- `--weebill-db-base` on an input with no database gives one with it.
- A genome with no marker hits is not added.
- End to end, marked expensive if slow: supplement, then
  `WeebillProfiler.annotate` on a fake profile whose `Genome_file` is the
  supplemented FASTA path maps to the supplemented taxonomy.

### Runtime

Sketching is roughly linear in total genome size. For about 2×10^5 GlobDB
genomes, expect hours on 32 threads, plus `db-add`. This should be in the log
("Sketching N genomes for weebill"), and `--no-weebill` covers users who do not
want it.

---

## #3 fast genome key → taxonomy lookup

### Problem

`MetapackageReadNameStore.get_taxonomy_by_genome_accession` (renamed
`get_taxonomy_by_genome_key` by #1) runs `SELECT read_name, taxonomy` over the
whole table. It applies regexes in Python to every row, on every `pipe`/`renew`
that integrates weebill: 11M rows for the GlobDB metapackage, several minutes
on slow storage. Only the few hundred genomes weebill reports are needed.

### Design

Add a table, written when the store is generated:

```
genome_taxonomy(genome_key VARCHAR PRIMARY KEY, taxonomy VARCHAR)
```

- **Build:** in `MetapackageReadNameStore.generate`, while streaming
  `(read_name, taxonomy)` pairs, compute `genome_key_from_read_name(read_name)`
  (#1) and keep the first taxonomy seen per key, as the current lookup does.
  About 10^5–10^6 keys, so a dict in memory is fine. Write it with `COPY` in the
  same way as the other tables, then add the primary key or an index.
- **Lookup:** if the table exists, run `SELECT genome_key, taxonomy FROM
  genome_taxonomy WHERE genome_key IN (...)` in chunks of 900, like
  `get_taxonomy_of_reads`. Otherwise fall back to the current full scan, so old
  metapackages (including SQLite v3–v6 ones) keep working. Detect the table
  once, through `information_schema` (DuckDB) or `sqlite_master` (SQLite).
- **Consistency:** taxonomies that disagree for one genome key should not happen,
  since one genome has one taxonomy. Log a warning while building if they do.
- **No metapackage version bump is needed:** the table is optional and readers
  detect it. Optionally record `genome_taxonomy_table: true` in `CONTENTS.json`
  for `singlem metapackage --describe`.

### Retrofitting existing metapackages

`singlem metapackage` only builds stores from scratch. Add a small admin path,
either `singlem metapackage --add-genome-taxonomy-table <mpkg>` or a script in
`extras/`. It would scan the existing store once, add the table in place for
DuckDB, or write a new DuckDB store and repoint `CONTENTS.json` for SQLite. Use
it for the default S7.5.0 metapackage and the GlobDB one. That could also become
a final rule in `extras/add_weebill_to_supplemented/Snakefile`.

### Tests

- A store generated from test spkgs has the table, and the lookup returns the
  same mapping as the full-scan fallback for GTDB-style and supplemented read
  names.
- A store without the table (simulate by dropping it) still works through the
  fallback.
- The existing `test_genome_accession_to_taxonomy` (skipped when the r232
  metapackage is absent) is unchanged.

### Order

#3 is independent of #2 and can land first. Do it before #2 if the GlobDB
metapackage will be used in production soon, since the scan cost grows with
every supplement.
