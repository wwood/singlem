#!/usr/bin/env python3
"""Estimate the per-base sequencing error rate of a metagenome from a SingleM OTU table.

Self-contained: needs only the Python standard library, and no other file from
this project. It never takes a true error rate as input.

How it works
------------
SingleM windows are 60 bp of protein-coding marker gene, in frame. A window that
differs from its correct sequence differs for one of two reasons:

  * a sequencing error, which is unselected and therefore non-synonymous with
    probability f_e ~ 0.75 (set by the genetic code and the transition bias); or
  * real biological variation between the organisms sharing that window, which
    is under purifying selection and therefore non-synonymous with a much lower
    probability f_b ~ 0.19 (dN/dS ~ 0.08).

If a fraction f_obs of the observed mismatches are non-synonymous and the
overall mismatch rate is m per base, then

    error_rate = m * (f_obs - f_b) / (f_e - f_b)

i.e. the non-synonymous fraction is used to work out how much of the observed
divergence is error and how much is biology. f_e and f_b were fitted on
simulated communities (see result.md); the defaults are baked in per mode.

Modes
-----
Which sequence a window is compared against is the only thing that differs
between modes.

  denovo    Reference-free. Windows of one marker are clustered by Hamming
            distance, most-abundant-first, and each cluster's most abundant
            member is taken as the reference. Needs nothing but the OTU table,
            so it applies to a metagenome with no relative in any database.

  hub       Reference-based consensus. Each observed window is assigned to the
            nearest reference window, and then scored against the most abundant
            reference window of that ortholog cluster. Biological variation
            within a cluster shows up as mismatches, which is the regime the
            model is designed for. Requires --reference-windows.

  nearest   Reference-based, nearest match. Each observed window is scored
            against the single closest reference window. Biological variation is
            absorbed by the match, so almost all remaining mismatches are
            errors. This is the most accurate mode when the reference genuinely
            contains the organisms present; it degrades towards under-estimation
            when it does not. Requires --reference-windows.

Accuracy (all three fitted on the same 144 conditions, per-base absolute error)
------------------------------------------------------------------------------
  nearest   median 1.9e-4, 90th pct 8.1e-4, worst 1.2e-3
  denovo    median 4.6e-4, 90th pct 2.8e-3, worst 7.5e-3
  hub       median 8.1e-4, 90th pct 2.5e-3, worst 5.3e-3

De novo beats hub at the median despite using no reference at all, but has the
worse tail: its error grows with the true rate (median 3.0e-3 at 3.1% per base,
where hub manages 1.4e-3), because at high error rates its clusters start to
fragment. Below ~1% per base it is the better of the two.

The absolute error does not shrink as the true rate falls, so treat the output
as a rate +/- that absolute floor, not as a percentage. An estimate below ~0.002
in hub/denovo mode should be read as "at or below 0.2%", not as a number. A
negative estimate means the observed mismatches were more synonymous than the
biological baseline, i.e. the error contribution is indistinguishable from zero.

Examples
--------
  # reference-free, one estimate per sample in the table
  estimate_error_rate.py --otu-table pipe.otu_table.tsv --mode denovo

  # both reference-based schemes, using windows extracted from a metapackage
  estimate_error_rate.py --otu-table pipe.otu_table.tsv --mode nearest,hub \\
      --reference-windows truth_windows.tsv --reference-groups acc_to_genus.tsv

  # everything, pooled across samples, written to a file
  estimate_error_rate.py -i pipe.otu_table.tsv --mode all --pool-samples \\
      --reference-windows truth_windows.tsv -o error_estimate.tsv
"""
import argparse
import collections
import csv
import gzip
import math
import sys

# --------------------------------------------------------------------------
# Genetic code
# --------------------------------------------------------------------------
_BASES = "TCAG"
_AAS = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODON_TABLE = {
    a + b + c: _AAS[i]
    for i, (a, b, c) in enumerate(
        (x, y, z) for x in _BASES for y in _BASES for z in _BASES
    )
}
PURINES = frozenset("AG")
PYRIMIDINES = frozenset("CT")

WINDOW_LEN = 60

# Fitted calibration constants, from results/fit/error_model.json of the
# simulated sweep (4 communities x 3 strain levels x 6 error levels x 2 reps).
#   f_e  non-synonymous fraction of sequencing errors
#   f_b  non-synonymous fraction of biological variation
#   floor  typical absolute uncertainty of the estimate, per base
CALIBRATION = {
    "nearest": {"f_e": 0.7395, "f_b": 0.1424, "floor": 0.0004,
                "note": "fitted on 144 conditions; R2=0.998"},
    "hub": {"f_e": 0.7482, "f_b": 0.1895, "floor": 0.0011,
            "note": "fitted on 144 conditions; R2=0.978"},
    # Fitted on the same 144 conditions as the other two. f_e comes out below
    # the genetic-code prior (0.687 vs 0.760) and falls as the communities get
    # more heterogeneous, because abundance-ordered clustering splits divergent
    # organisms into separate clusters instead of pooling them: some biological
    # variation never becomes a mismatch, and the fitted f_e absorbs the
    # shortfall. So this f_e is an effective constant, not a code constant.
    "denovo": {"f_e": 0.6871, "f_b": 0.0104, "floor": 0.0010,
               "note": "fitted on 144 conditions; R2=0.976; f_e is effective, "
                       "not a genetic-code constant"},
}
MODES = ("denovo", "nearest", "hub")


# --------------------------------------------------------------------------
# Basic helpers
# --------------------------------------------------------------------------
def is_transition(a, b):
    return ({a, b} <= PURINES) or ({a, b} <= PYRIMIDINES)


def hamming(a, b):
    return sum(1 for x, y in zip(a, b) if x != y)


def hamming_le(a, b, limit):
    """Hamming distance, abandoning the count once it exceeds `limit`.

    Returns limit+1 as soon as that is certain. Matching an observed window
    against thousands of references is the hot loop of this script, and most
    comparisons are rejected within the first few bases.
    """
    d = 0
    for x, y in zip(a, b):
        if x != y:
            d += 1
            if d > limit:
                return d
    return d


def read_otu_table(path):
    """Tolerant reader for SingleM OTU tables (plain or gzipped).

    Accepts the usual column spellings; anything with a gene/marker column and a
    sequence column will do. Count falls back to 1.0 per row if absent.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    handle = sys.stdin if str(path) == "-" else opener(path, "rt")
    try:
        reader = csv.DictReader(handle, delimiter="\t")
        cols = {c.lower(): c for c in (reader.fieldnames or [])}

        def pick(*names):
            for n in names:
                if n in cols:
                    return cols[n]
            return None

        gene_c = pick("gene", "marker")
        seq_c = pick("sequence", "otu_sequence")
        n_c = pick("num_hits", "count", "coverage")
        samp_c = pick("sample", "sample_name")
        if not (gene_c and seq_c):
            sys.exit(f"{path}: no gene/sequence columns found in "
                     f"{reader.fieldnames}")

        for r in reader:
            yield {
                "gene": r[gene_c],
                "sequence": (r[seq_c] or "").upper(),
                "count": float(r[n_c]) if n_c and r.get(n_c) else 1.0,
                "sample": r[samp_c] if samp_c else "",
            }
    finally:
        if handle is not sys.stdin:
            handle.close()


# --------------------------------------------------------------------------
# Mismatch classification
# --------------------------------------------------------------------------
def classify_pair(ref, obs, weight, acc):
    """Accumulate synonymous / non-synonymous counts for one ref-obs pair.

    Each mismatching position is mutated into the reference codon on its own and
    the amino acid compared. Multiple mismatches in one window are treated
    independently against the reference codon background, which is adequate at
    the mismatch densities of interest (< ~5% per base).
    """
    for i, (r, o) in enumerate(zip(ref, obs)):
        if r == o or o not in "ACGT" or r not in "ACGT":
            continue
        ci, cp = divmod(i, 3)
        ref_codon = ref[3 * ci: 3 * ci + 3]
        if len(ref_codon) != 3 or set(ref_codon) - set("ACGT"):
            continue
        aa_ref = CODON_TABLE[ref_codon]
        aa_mut = CODON_TABLE[ref_codon[:cp] + o + ref_codon[cp + 1:]]

        if aa_mut == "*" and aa_ref != "*":
            kind = "nonsense"
        elif aa_ref == aa_mut:
            kind = "synonymous"
        else:
            kind = "nonsynonymous"

        acc[kind] += weight
        acc[f"pos{cp + 1}"] += weight
        acc["transition" if is_transition(r, o) else "transversion"] += weight
        acc["mismatches"] += weight


def prior_nonsyn_fraction(codon_weights, ti_tv):
    """Non-synonymous fraction expected if every mismatch were a random error.

    Enumerates all three single-base changes at all three positions of every
    reference codon actually used, weighted by the observed transition /
    transversion ratio. Nonsense changes are excluded from both numerator and
    denominator so this is directly comparable to the reported f_nonsyn.

    This is a prediction of f_e from codon composition alone, with no fitting.
    On the calibration data it agreed with the fitted f_e to within 0.02, so a
    large disagreement here is a warning that the windows may be out of frame or
    the reference inappropriate.
    """
    if ti_tv and ti_tv == ti_tv and math.isfinite(ti_tv):
        ti_w = ti_tv / (ti_tv + 2.0)
    else:
        ti_w = 1 / 3
    tv_w = (1.0 - ti_w) / 2.0

    num = den = nonsense = all_den = 0.0
    for codon, cw in codon_weights.items():
        aa = CODON_TABLE[codon]
        for cp in range(3):
            for alt in "ACGT":
                if alt == codon[cp]:
                    continue
                w = ti_w if is_transition(codon[cp], alt) else tv_w
                mut_aa = CODON_TABLE[codon[:cp] + alt + codon[cp + 1:]]
                all_den += cw * w
                if mut_aa == "*" and aa != "*":
                    nonsense += cw * w
                    continue
                den += cw * w
                if mut_aa != aa:
                    num += cw * w
    return (num / den if den else float("nan"),
            nonsense / all_den if all_den else float("nan"))


# --------------------------------------------------------------------------
# Reference construction: de novo
# --------------------------------------------------------------------------
def build_denovo_refs(rows, max_dist):
    """Cluster observed windows and pick a reference for each -> {(gene, seq): ref}.

    Greedy, most-abundant-first: the most abundant unassigned window of a marker
    seeds a cluster and absorbs every window within max_dist of it. Abundance
    ordering is what makes this work without a database — the error-free sequence
    is more abundant than any of its error variants, so it becomes the seed and
    the variants attach to it.
    """
    by_gene = collections.defaultdict(collections.Counter)
    for r in rows:
        by_gene[r["gene"]][r["sequence"]] += r["count"]

    ref_of = {}
    n_clusters = 0
    for gene, weight in by_gene.items():
        remaining = sorted(weight, key=lambda s: (-weight[s], s))
        while remaining:
            seed = remaining[0]
            members = [s for s in remaining if hamming_le(s, seed, max_dist) <= max_dist]
            for s in members:
                ref_of[(gene, s)] = seed
            n_clusters += 1
            drop = set(members)
            remaining = [s for s in remaining if s not in drop]
    return ref_of, n_clusters


# --------------------------------------------------------------------------
# Reference construction: from reference windows
# --------------------------------------------------------------------------
def cluster_orthologs(seqs, max_dist, owners):
    """Split the reference windows of one marker into orthologous clusters.

    Marker genes are frequently present as several divergent copies per genome
    (paralogs, or split HMM hits). Measured on GTDB r232, same-genome same-gene
    window pairs differ by ~28 of 60 bases at the median, whereas orthologous
    windows across species of one genus differ by ~10. Without this split, one
    hub reference per marker would be a paralog for a large share of reads and
    every position of the wrong copy would count as a mismatch, inflating the
    apparent mismatch rate several-fold.

    Complete linkage — every member within max_dist of every other — plus a hard
    rule that two windows carried by the same genome never join. Single linkage
    was tried first and chained paralogs together through intermediate windows,
    leaving a residual mismatch-rate floor even with no biological variation at
    all.

    `owners` maps sequence -> set of genomes carrying it.

    Returns a list of clusters, each a list of sequences.
    """
    seqs = sorted(seqs)
    dist = {}
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            d = hamming(seqs[i], seqs[j])
            dist[(seqs[i], seqs[j])] = dist[(seqs[j], seqs[i])] = d

    def compatible(a, b):
        if dist[(a, b)] > max_dist:
            return False
        return not (owners.get(a, frozenset()) & owners.get(b, frozenset()))

    remaining = list(seqs)
    clusters = []
    while remaining:
        # start from the sequence with the most compatible partners, so the
        # densest ortholog group forms first
        seed = max(remaining,
                   key=lambda s: (sum(1 for o in remaining
                                      if o != s and compatible(s, o)), s))
        members = [seed]
        for s in sorted(remaining,
                        key=lambda o: (dist[(seed, o)] if o != seed else -1)):
            if s != seed and all(compatible(s, m) for m in members):
                members.append(s)
        clusters.append(members)
        remaining = [s for s in remaining if s not in set(members)]
    return clusters


def build_reference(path, group_map, ortholog_max_dist, verbose=False):
    """Index reference windows for nearest and hub matching.

    -> candidates: {gene: [(sequence, cluster_key)]}
       hubs:       {cluster_key: hub_sequence}

    A cluster_key is (gene, group, cluster_index): one ortholog of one group of
    related genomes. Clustering runs per marker across all groups at once, so an
    ortholog cluster can span groups; the group then subdivides it. A window
    carried by more than one group is given no cluster_key at all, because it
    cannot say which group's consensus a read matching it should be scored
    against — such windows are still available as nearest matches but are not
    hub-eligible.

    The hub is the cluster's most abundant window, abundance here being the
    number of reference genomes carrying it. (The calibration weighted by the
    simulated relative abundances, which a real reference does not have; genome
    count is the closest available proxy and is the one place this script cannot
    reproduce the calibration exactly.)
    """
    per_gene = collections.defaultdict(dict)        # gene -> seq -> {group: n}
    owners = collections.defaultdict(               # gene -> seq -> {genomes}
        lambda: collections.defaultdict(set))
    n_rows = 0
    for row in read_otu_table(path):
        seq = row["sequence"]
        if len(seq) != WINDOW_LEN or set(seq) - set("ACGT"):
            continue
        genome = row["sample"] or "unknown"
        group = group_map.get(genome, group_map.get(_bare(genome), "all"))
        d = per_gene[row["gene"]].setdefault(seq, collections.Counter())
        d[group] += 1
        owners[row["gene"]][seq].add(genome)
        n_rows += 1

    candidates = collections.defaultdict(list)
    hubs = {}
    best = {}
    n_clusters = 0
    for gene, seqs in per_gene.items():
        if verbose and len(seqs) > 3000:
            print(f"# clustering {len(seqs)} reference windows for {gene}; "
                  f"this is O(n^2)", file=sys.stderr)
        own = {s: frozenset(g) for s, g in owners[gene].items()}
        clusters = cluster_orthologs(list(seqs), ortholog_max_dist, own)
        seq_cluster = {}
        for ci, members in enumerate(clusters):
            n_clusters += 1
            for s in members:
                seq_cluster[s] = ci
        for seq, groups in seqs.items():
            group = next(iter(groups)) if len(groups) == 1 else None
            key = (gene, group, seq_cluster[seq]) if group is not None else None
            candidates[gene].append((seq, key))
            if key is not None:
                w = sum(groups.values())
                if w > best.get(key, (None, -1))[1]:
                    best[key] = (seq, w)
    for k, (seq, _w) in best.items():
        hubs[k] = seq
    return dict(candidates), hubs, n_rows, n_clusters


def _bare(accession):
    """GTDB accessions appear as GB_GCA_000008085.1 and as GCA_000008085.1."""
    return accession.split("_", 1)[1] if accession[:3] in ("GB_", "RS_") else accession


def read_group_map(path):
    """accession -> group label, from a 2-column TSV (with or without a header).

    Groups bundle genomes whose windows should share one consensus reference —
    a GTDB genus is what the calibration used. Both the bare and the GB_/RS_
    prefixed forms of an accession are registered so either spelling matches.
    """
    out = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            acc, group = parts[0].strip(), parts[1].strip()
            if acc.lower() in ("accession", "genome", "sample", "sample_name"):
                continue
            for form in (acc, _bare(acc), f"GB_{acc}", f"RS_{acc}"):
                out.setdefault(form, group)
    return out


def match_nearest(seq, candidates, max_mm):
    """-> (reference_sequence, cluster_key), or (None, None) if absent/ambiguous.

    Ambiguous means two different reference sequences tie for closest; such a
    window carries no usable information about which reference it came from and
    is dropped.
    """
    best_d, best = max_mm + 1, []
    for cseq, ckey in candidates:
        d = hamming_le(seq, cseq, best_d)
        if d < best_d:
            best_d, best = d, [(cseq, ckey)]
        elif d == best_d and d <= max_mm:
            best.append((cseq, ckey))
    if best_d > max_mm or not best:
        return None, None
    if len({b[0] for b in best}) > 1:
        return None, None
    return best[0]


# --------------------------------------------------------------------------
# Estimation
# --------------------------------------------------------------------------
def estimate(rows, mode, f_e, f_b, args, reference=None):
    """Score every window of one sample and return a result dict."""
    acc = collections.Counter()
    codon_weights = collections.Counter()
    total_weight = unmatched_weight = 0.0
    n_windows = n_clusters = 0

    if mode == "denovo":
        ref_of, n_clusters = build_denovo_refs(rows, args.denovo_max_distance)
    else:
        candidates, hubs, _, n_clusters = reference

    for r in rows:
        seq, gene, w = r["sequence"], r["gene"], r["count"]
        n_windows += 1
        total_weight += w

        if mode == "denovo":
            ref = ref_of[(gene, seq)]
        else:
            cands = candidates.get(gene)
            if not cands:
                unmatched_weight += w
                continue
            nearest, key = match_nearest(seq, cands, args.max_mismatches)
            if nearest is None:
                unmatched_weight += w
                continue
            ref = nearest if mode == "nearest" else hubs.get(key)
            if ref is None:
                unmatched_weight += w
                continue

        acc["windows"] += w
        acc["bases"] += w * WINDOW_LEN
        for ci in range(WINDOW_LEN // 3):
            c = ref[3 * ci:3 * ci + 3]
            if not set(c) - set("ACGT"):
                codon_weights[c] += w
        classify_pair(ref, seq, w, acc)

    mm = acc["mismatches"]
    syn, non, nons = acc["synonymous"], acc["nonsynonymous"], acc["nonsense"]
    bases = acc["bases"]
    titv = (acc["transition"] / acc["transversion"]
            if acc["transversion"] else float("nan"))
    # nonsense is excluded from f_nonsyn: SingleM may drop stop-containing
    # windows, making that class censoring-prone
    f_obs = non / (syn + non) if (syn + non) else float("nan")
    m = mm / bases if bases else float("nan")
    est = m * (f_obs - f_b) / (f_e - f_b) if bases and (syn + non) else float("nan")
    prior_fe, prior_fns = prior_nonsyn_fraction(codon_weights, titv)

    return {
        "mode": mode,
        "estimated_error_rate": est,
        "estimated_error_rate_pct": est * 100 if est == est else float("nan"),
        "absolute_uncertainty": CALIBRATION[mode]["floor"],
        "n_window_rows": n_windows,
        "scored_weight": acc["windows"],
        "unmatched_weight": unmatched_weight,
        "total_weight": total_weight,
        "n_reference_clusters": n_clusters,
        "scored_bases": bases,
        "mismatches": mm,
        "mismatch_rate": m,
        "n_synonymous": syn,
        "n_nonsynonymous": non,
        "n_nonsense": nons,
        "f_nonsyn": f_obs,
        "f_nonsense": nons / mm if mm else float("nan"),
        "codon_pos1": acc["pos1"] / mm if mm else float("nan"),
        "codon_pos2": acc["pos2"] / mm if mm else float("nan"),
        "codon_pos3": acc["pos3"] / mm if mm else float("nan"),
        "ti_tv": titv,
        "f_e_used": f_e,
        "f_b_used": f_b,
        "f_e_prior_from_codons": prior_fe,
        "f_nonsense_prior_from_codons": prior_fns,
    }


FIELDS = ["sample", "mode", "estimated_error_rate", "estimated_error_rate_pct",
          "absolute_uncertainty", "mismatch_rate", "f_nonsyn", "f_e_used",
          "f_b_used", "f_e_prior_from_codons", "ti_tv", "codon_pos1",
          "codon_pos2", "codon_pos3", "f_nonsense",
          "f_nonsense_prior_from_codons", "mismatches", "scored_bases",
          "n_synonymous", "n_nonsynonymous", "n_nonsense", "n_window_rows",
          "scored_weight", "unmatched_weight", "total_weight",
          "n_reference_clusters"]


def fmt(v):
    if isinstance(v, float):
        return "nan" if v != v else f"{v:.6g}"
    return str(v)


# --------------------------------------------------------------------------
def build_parser():
    p = argparse.ArgumentParser(
        description=__doc__.split("Examples")[0],
        epilog="Examples:\n" + __doc__.split("Examples\n--------\n")[1],
        formatter_class=argparse.RawDescriptionHelpFormatter)

    io_g = p.add_argument_group("input / output")
    io_g.add_argument(
        "-i", "--otu-table", required=True, metavar="FILE",
        help="SingleM OTU table to estimate the error rate of, as written by "
             "`singlem pipe --otu-table`. Plain or .gz; '-' reads stdin. Needs "
             "a gene/marker column, a sequence column, and ideally num_hits and "
             "sample columns. Rows whose sequence is not exactly 60 bp are "
             "ignored. [required]")
    io_g.add_argument(
        "-o", "--output", metavar="FILE", default="-",
        help="Where to write the results TSV, one row per sample per mode. "
             "Default '-' for stdout.")
    io_g.add_argument(
        "--pool-samples", action="store_true",
        help="Treat every row of the OTU table as coming from one sequencing "
             "run and report a single estimate. By default each value of the "
             "sample column is estimated separately, which is usually what you "
             "want since the error rate is a property of a run.")
    io_g.add_argument(
        "--min-count", type=float, default=0.0, metavar="N",
        help="Discard window rows with num_hits below N before doing anything. "
             "USE WITH CARE: error variants are exactly the rare windows, so "
             "any positive value biases the estimate downwards. Only useful to "
             "make very large tables tractable. Default 0 (keep everything).")

    m_g = p.add_argument_group("what to estimate")
    m_g.add_argument(
        "--mode", default="denovo", metavar="LIST",
        help="Comma-separated list of reference schemes to use: 'denovo' "
             "(reference-free, needs nothing else), 'nearest' and/or 'hub' "
             "(both need --reference-windows), or 'all' for every applicable "
             "mode. See the description above for what each one measures and "
             "how accurate it is. Default 'denovo'.")

    r_g = p.add_argument_group("reference windows (for nearest / hub modes)")
    r_g.add_argument(
        "--reference-windows", metavar="FILE",
        help="OTU-table-formatted file of known-correct windows, one row per "
             "(genome, marker, sequence), with the genome accession in the "
             "sample column. For GTDB these can be read straight out of a "
             "SingleM metapackage's transcripts.sdb rather than recomputed. "
             "Required for --mode nearest and --mode hub.")
    r_g.add_argument(
        "--reference-groups", metavar="FILE",
        help="Two-column TSV mapping genome accession -> group label, which "
             "bundles genomes whose windows should share one consensus "
             "reference in hub mode. The calibration used GTDB genus. If "
             "omitted, all reference genomes go into a single group, which "
             "makes hub references coarser and biases hub estimates upwards. "
             "Ignored in nearest and denovo modes.")
    r_g.add_argument(
        "--max-mismatches", type=int, default=8, metavar="N",
        help="A window more than N bases from every reference window is "
             "unassignable and is dropped (counted in unmatched_weight). Too "
             "small and genuinely error-rich reads are lost, biasing the "
             "estimate down; too large and unrelated organisms get matched. "
             "Default 8, as calibrated.")
    r_g.add_argument(
        "--ortholog-max-distance", type=int, default=15, metavar="N",
        help="Complete-linkage Hamming cutoff that splits paralogous copies of "
             "a marker into separate reference clusters, so a hub reference is "
             "never a paralog of the read it scores. On GTDB r232 same-genome "
             "copies differ by ~28 bases and within-genus orthologs by ~10, so "
             "15 separates them cleanly. Default 15.")

    d_g = p.add_argument_group("de novo clustering")
    d_g.add_argument(
        "--denovo-max-distance", type=int, default=6, metavar="N",
        help="Hamming radius within which an observed window joins the cluster "
             "of a more abundant window. Must be wide enough to catch the error "
             "variants of a sequence but narrow enough not to merge distinct "
             "organisms. Default 6, which is what the spot-check used; raise it "
             "for very high error rates (>5%%), lower it for densely sampled "
             "communities of close relatives.")

    c_g = p.add_argument_group("calibration constants")
    c_g.add_argument(
        "--f-e", type=float, metavar="X",
        help="Override the non-synonymous fraction of sequencing errors. "
             "Defaults to the fitted value for the mode (nearest 0.7395, "
             "hub/denovo 0.7482). Only change this if you have recalibrated, or "
             "to substitute the f_e_prior_from_codons value the run reports.")
    c_g.add_argument(
        "--f-b", type=float, metavar="X",
        help="Override the non-synonymous fraction of biological variation. "
             "Defaults to the fitted value for the mode (nearest 0.1424, "
             "hub/denovo 0.1895); equivalently dN/dS = f_b/(3*(1-f_b)). Set it "
             "to 0 if you believe the windows contain no real variation at all, "
             "e.g. a single-isolate sequencing run.")
    c_g.add_argument(
        "--verbose", action="store_true",
        help="Report progress and warnings on stderr.")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    modes = ["denovo", "nearest", "hub"] if args.mode.strip() == "all" else [
        m.strip() for m in args.mode.split(",") if m.strip()]
    bad = [m for m in modes if m not in MODES]
    if bad:
        sys.exit(f"unknown mode(s) {bad}; choose from {', '.join(MODES)} or 'all'")
    needs_ref = [m for m in modes if m != "denovo"]
    if needs_ref and not args.reference_windows:
        if args.mode.strip() == "all":
            if args.verbose:
                print("# no --reference-windows given; running denovo only",
                      file=sys.stderr)
            modes = ["denovo"]
            needs_ref = []
        else:
            sys.exit(f"--mode {','.join(needs_ref)} requires --reference-windows")

    # Length is the only filter: windows containing ambiguity codes are kept,
    # because classify_pair ignores individual non-ACGT positions rather than
    # the whole window, and this is what the calibration did.
    rows = [r for r in read_otu_table(args.otu_table)
            if len(r["sequence"]) == WINDOW_LEN
            and r["count"] >= args.min_count]
    if not rows:
        sys.exit(f"{args.otu_table}: no usable 60 bp windows found")
    if args.verbose:
        print(f"# {len(rows)} window rows loaded", file=sys.stderr)

    reference = None
    if needs_ref:
        group_map = read_group_map(args.reference_groups) if args.reference_groups else {}
        reference = build_reference(args.reference_windows, group_map,
                                   args.ortholog_max_distance, args.verbose)
        if args.verbose:
            print(f"# reference: {reference[2]} windows -> {reference[3]} "
                  f"ortholog clusters", file=sys.stderr)
        if not args.reference_groups and "hub" in modes and args.verbose:
            print("# warning: no --reference-groups, so all reference genomes "
                  "share one group; hub estimates will read high",
                  file=sys.stderr)

    if args.pool_samples:
        by_sample = {"pooled": rows}
    else:
        by_sample = collections.defaultdict(list)
        for r in rows:
            by_sample[r["sample"] or "unknown"].append(r)

    out = sys.stdout if args.output == "-" else open(args.output, "w")
    try:
        out.write("\t".join(FIELDS) + "\n")
        for sample, srows in by_sample.items():
            for mode in modes:
                cal = CALIBRATION[mode]
                f_e = args.f_e if args.f_e is not None else cal["f_e"]
                f_b = args.f_b if args.f_b is not None else cal["f_b"]
                if f_e - f_b < 0.05:
                    sys.exit(f"f_e ({f_e}) and f_b ({f_b}) are too close: the "
                             "estimator divides by their difference, so errors "
                             "and biology must be distinguishable")
                res = estimate(srows, mode, f_e, f_b, args, reference)
                res["sample"] = sample
                out.write("\t".join(fmt(res[k]) for k in FIELDS) + "\n")
                if args.verbose:
                    print(f"# {sample}\t{mode}\t"
                          f"{res['estimated_error_rate']:.3g} "
                          f"+/- {res['absolute_uncertainty']:.3g}"
                          f"\t({cal['note']})", file=sys.stderr)
    finally:
        if out is not sys.stdout:
            out.close()


if __name__ == "__main__":
    main()
