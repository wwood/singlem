import logging

import numpy as np

from .condense import (
    WordNode,
    CondensedCommunityProfile,
    _canonical_species_key,
    _gtdb_string_to_wordnode_array,
)


def _is_species_string(taxonomy):
    '''True if the taxonomy string resolves to species level (last rank s__).'''
    key = _canonical_species_key(taxonomy)
    if key == '':
        return False
    return key.split(';')[-1].startswith('s__')


class JointDeconvolver:
    '''Joint SingleM + sylph taxonomic profiling by penalised, weighted,
    non-negative least-squares deconvolution.

    Replaces condense's heuristic EM + trimmed-mean tree-condense with a single
    optimisation per sample:

        minimise_{a>=0}  || W^(1/2) (b - A a) ||_2^2  +  l1_penalty * sum(a)

    where a is genome-equivalent coverage over an augmented GTDB tree (one column
    per candidate species plus one "novel-at-clade" column per clade SingleM
    signal only resolves to). The observation blocks are the SingleM marker
    coverages (b = M a), the sylph effective coverages (b = alpha * a), and
    high-weight zero rows for DB-candidate species sylph did not report. alpha is
    fit by variable projection; SingleM rows are robustified by IRLS.'''

    def solve(self, sample, sample_otus, sylph_hits, domain_marker_counts=None, alpha=None,
              l1_penalty=1.0, absence_weight=100.0, sylph_weight=50.0, min_markers=3,
              max_outer_iterations=25, tolerance=1e-4, prune_below=0.05,
              min_singlem_coverage=0.35, robust_coverage_floor=1.0, robust_min_weight=1e-3,
              alpha_min_coverage=0.1, coherence_weight=1000.0):
        try:
            from scipy.optimize import minimize
            from scipy.sparse import csr_matrix
        except ImportError:
            raise Exception("condense --joint requires scipy, which could not be imported")

        columns, singlem_rows, observed_marker_count, unique_marker_count, unique_marker_coverage = \
            self._build_columns_and_singlem_rows(sample_otus, sylph_hits)
        if len(columns) == 0:
            logging.warning("Sample {}: no SingleM or sylph taxa to deconvolve".format(sample))
            return CondensedCommunityProfile(sample, WordNode(None, 'Root'))

        num_columns = len(columns)

        # Marker-count normalisation: each column's coverage is averaged over its
        # taxon's full single-copy marker complement, so a column observed on only
        # a few markers is penalised toward zero by its unobserved markers (the
        # zero-padding the trimmed-mean condense applies). padding_weight[c] is the
        # number of unobserved markers for column c.
        domain_marker_counts = domain_marker_counts or {}
        padding_weight = np.zeros(num_columns)
        for idx, column in enumerate(columns):
            domain = column.key.split(';')[0].replace('d__', '')
            expected = domain_marker_counts.get(domain, observed_marker_count[idx])
            padding_weight[idx] = max(0.0, float(expected) - observed_marker_count[idx])

        # SingleM design matrix M (coefficients all 1) and target b_M.
        row_idx, col_idx = [], []
        b_singlem = np.empty(len(singlem_rows), dtype=float)
        for r, (cols, coverage, _) in enumerate(singlem_rows):
            b_singlem[r] = coverage
            for c in cols:
                row_idx.append(r)
                col_idx.append(c)
        num_singlem = len(singlem_rows)
        if num_singlem > 0:
            M = csr_matrix((np.ones(len(row_idx)), (row_idx, col_idx)),
                           shape=(num_singlem, num_columns))
            Mt = M.transpose().tocsr()
        else:
            M = Mt = None

        # Sylph presence rows (coefficient alpha) and absence rows (zero target).
        sylph_col, sylph_eff, absence_col = [], [], []
        for idx, column in enumerate(columns):
            if column.kind != 'species':
                continue
            hit = sylph_hits.get(column.key)
            if hit is not None:
                sylph_col.append(idx)
                sylph_eff.append(hit.eff_cov)
            else:
                absence_col.append(idx)
        sylph_col = np.array(sylph_col, dtype=int)
        sylph_eff = np.array(sylph_eff, dtype=float)
        sylph_w = np.full(len(sylph_col), float(sylph_weight))
        absence_col = np.array(absence_col, dtype=int)

        centres = self._column_marker_centres(columns, unique_marker_coverage, num_columns, min_markers)

        # Alpha, the scale between sylph's coverage units and SingleM's, is anchored to
        # the species' unambiguous marker coverages -- evidence that does not depend on
        # the fit. Estimating it from the fitted coverages instead is degenerate as soon
        # as the model defers to sylph: the fit then sits at a = e/alpha, which satisfies
        # the sylph rows for *any* alpha, so alpha is pinned only by whatever weight the
        # SingleM rows still carry, and it drifts. Alpha drifting is not a cosmetic
        # problem: the sylph-supported species are placed on sylph's scale while the
        # novel columns, which sylph cannot see, stay on SingleM's, and a profile whose
        # taxa are on two different scales has the wrong composition even when every
        # individual coverage looks plausible.
        if alpha is None:
            current_alpha = 1.0
            fit_alpha = len(sylph_col) > 0
            # Any unambiguous marker will do to anchor the scale, so this does not wait
            # for the min_markers the coherence ceiling insists on: a noisy anchor from
            # one marker is worth more than an alpha free to drift.
            ratios = []
            for i, c in enumerate(sylph_col):
                markers = unique_marker_coverage.get(int(c))
                if not markers:
                    continue
                singlem_coverage = float(np.median(list(markers.values())))
                if singlem_coverage > 0:
                    ratios.append(sylph_eff[i] / singlem_coverage)
            if len(ratios) > 0:
                current_alpha = float(np.median(ratios))
                fit_alpha = False
                logging.info("Anchored alpha={:.4f} to the unambiguous marker coverage of {} "
                             "sylph species".format(current_alpha, len(ratios)))
        else:
            current_alpha = float(alpha)
            fit_alpha = False
        logging.info("Joint deconvolution of sample {}: {} columns ({} sylph species, {} SingleM-only "
            "candidates), {} SingleM rows".format(sample, num_columns, len(sylph_col), len(absence_col), num_singlem))

        # Identifiability floor: a taxon with no sylph support must be uniquely
        # resolved by at least min_markers markers, otherwise its coverage is
        # fixed to zero. This suppresses false positives that ride on a single
        # uniquely-assigned marker amid windows shared with a sylph-supported
        # neighbour. Sylph-detected species (including sylph-only injections) are
        # exempt, since sylph itself provides the evidence.
        #
        # Novel columns are held to the same floor. Nothing else constrains them --
        # they have no sylph observation to answer to and no absence row -- so a novel
        # column is the cheapest place in the model for coverage that the species
        # columns cannot explain, and one that no read resolved to will still collect
        # whatever is left over. Requiring reads that resolve to the clade and no
        # deeper is what earns it a coverage.
        sylph_columns = set(int(c) for c in sylph_col)
        bounds = [(0.0, None)] * num_columns
        num_fixed = 0
        if min_markers > 0:
            for j, column in enumerate(columns):
                if j not in sylph_columns and unique_marker_count[j] < min_markers:
                    bounds[j] = (0.0, 0.0)
                    num_fixed += 1
        if num_fixed > 0:
            logging.info("Fixed {} sylph-unsupported columns to zero (< {} unique markers)".format(
                num_fixed, min_markers))

        a = np.zeros(num_columns)
        singlem_weights = np.ones(num_singlem)

        def optimise(current_a, current_alpha, singlem_weights):
            """Optimise once for the currently permitted candidate columns."""
            a = current_a
            def objective(x):
                f = 0.0
                g = np.zeros(num_columns)
                if num_singlem > 0:
                    residual = b_singlem - M.dot(x)
                    weighted_residual = singlem_weights * residual
                    f += float(np.dot(weighted_residual, residual))
                    g += -2.0 * Mt.dot(weighted_residual)
                if len(sylph_col) > 0:
                    a_sylph = x[sylph_col]
                    residual_sylph = sylph_eff - current_alpha * a_sylph
                    f += float(np.dot(sylph_w, residual_sylph * residual_sylph))
                    np.add.at(g, sylph_col, -2.0 * sylph_w * current_alpha * residual_sylph)
                if len(absence_col) > 0:
                    a_absent = x[absence_col]
                    f += absence_weight * float(np.dot(a_absent, a_absent))
                    np.add.at(g, absence_col, 2.0 * absence_weight * a_absent)
                # Marker-count normalisation (zero-padding of unobserved markers).
                f += float(np.dot(padding_weight, x * x))
                g += 2.0 * padding_weight * x
                # Marker coherence: a taxon's coverage may not exceed what its own
                # markers agree on. Coverage is near-uniform across a genome's
                # single-copy markers, so a column whose markers mostly read 3x is a
                # 3x organism even if two of them read 1100x -- those two are being
                # fed by something else that shares the window. Least squares is a
                # mean, so without this the two spikes carry the column away, and the
                # row weighting cannot see it: judged against the whole sample, 1100x
                # is merely abundant, and it is only against the column's *own* other
                # markers that it is absurd. One-sided, because a column sitting below
                # its markers' centre is the sylph and absence terms doing their job.
                excess = np.maximum(x - centres, 0.0)
                f += coherence_weight * float(np.dot(excess, excess))
                g += 2.0 * coherence_weight * excess
                f += l1_penalty * float(np.sum(x))
                g += l1_penalty
                return f, g

            for outer in range(max_outer_iterations):
                result = minimize(objective, a, jac=True, method='L-BFGS-B', bounds=bounds)
                new_a = np.maximum(result.x, 0.0)

                # Alpha update: the median of the per-species coverage ratios, not a
                # least-squares projection. Least squares weights each species by its
                # coverage squared, so the single most abundant species sets alpha
                # almost alone -- and if that species is over-fit (ambiguous coverage
                # in a clade has to land somewhere), alpha falls to accommodate it,
                # which loosens the sylph constraint on every other species and lets
                # the over-fit grow. The median of ratios is scale-invariant and
                # cannot be moved by one species, so alpha stays a calibration
                # between sylph's and SingleM's coverage units rather than becoming a
                # free parameter that absorbs the model's mistakes.
                new_alpha = current_alpha
                if fit_alpha and len(sylph_col) > 0:
                    a_sylph = new_a[sylph_col]
                    usable = a_sylph > alpha_min_coverage
                    if np.any(usable):
                        new_alpha = float(np.median(sylph_eff[usable] / a_sylph[usable]))
                    if not np.isfinite(new_alpha) or new_alpha <= 0:
                        new_alpha = current_alpha

                # IRLS reweighting of SingleM rows (bisquare on relative residuals).
                if num_singlem > 0:
                    singlem_weights = self._bisquare_weights(
                        b_singlem - M.dot(new_a), b_singlem,
                        coverage_floor=robust_coverage_floor, min_weight=robust_min_weight)

                converged = (np.max(np.abs(new_a - a)) < tolerance and abs(new_alpha - current_alpha) < tolerance)
                a, current_alpha = new_a, new_alpha
                if converged:
                    logging.debug("Joint deconvolution converged after {} iterations".format(outer + 1))
                    break
            return a, current_alpha, singlem_weights

        # First obtain provisional coverages for every eligible column. Then
        # remove low-coverage SingleM-only candidates and solve again. Repeat
        # until the active set is stable. Sylph-supported species are exempt:
        # genome-wide evidence can support them below SingleM's general taxon
        # coverage threshold.
        while True:
            a, current_alpha, singlem_weights = optimise(a, current_alpha, singlem_weights)
            newly_fixed = []
            if min_singlem_coverage is not None and min_singlem_coverage > 0:
                for j in range(num_columns):
                    if j not in sylph_columns and bounds[j] != (0.0, 0.0) and a[j] < min_singlem_coverage:
                        bounds[j] = (0.0, 0.0)
                        newly_fixed.append(j)
            if len(newly_fixed) == 0:
                break
            a[newly_fixed] = 0.0
            logging.info("Removed {} SingleM-only candidates below {:.2f} coverage; re-optimising".format(
                len(newly_fixed), min_singlem_coverage))

        logging.info("Joint deconvolution of sample {}: alpha={:.4f}, total coverage={:.2f}".format(
            sample, current_alpha, float(np.sum(a))))

        # Use the lower numerical/output floor only for sylph-supported species.
        # SingleM-only calls at any rank must satisfy the normal taxon threshold.
        for j in range(num_columns):
            threshold = prune_below if j in sylph_columns else min_singlem_coverage
            if threshold is not None and a[j] < threshold:
                a[j] = 0.0

        # Stash the solution for diagnostics / testing.
        self.fitted_alpha = current_alpha
        self.columns = columns
        self.solution = a
        self.coverage_by_key = {column.key: cov for column, cov in zip(columns, a)}

        return self._build_profile(sample, columns, a)

    def _build_columns_and_singlem_rows(self, sample_otus, sylph_hits):
        '''Return (columns, singlem_rows). columns is a list of _Column; each
        singlem row is (sorted_list_of_column_indices, coverage).

        A window resolving only to a clade (DIAMOND/genus hit, LCA-spanning query
        hit) loads onto that clade's novel leaf AND every candidate species within
        the clade, so the deconvolution can attribute a genus-resolved read to a
        known species that other evidence supports rather than double-counting it
        as novel.'''
        columns = []
        species_key_to_index = {}
        novel_key_to_index = {}

        def species_column(taxonomy):
            key = _canonical_species_key(taxonomy)
            if key not in species_key_to_index:
                species_key_to_index[key] = len(columns)
                columns.append(_Column('species', taxonomy, key))
            return species_key_to_index[key]

        def novel_column(taxonomy):
            key = _canonical_species_key(taxonomy)
            if key not in novel_key_to_index:
                novel_key_to_index[key] = len(columns)
                columns.append(_Column('novel', taxonomy, key))
            return novel_key_to_index[key]

        # markers observed per column: counted only from direct hits at the
        # column's own level (species-level for species columns, clade-level for
        # novel columns), so a species seen only via genus-level reads counts as
        # observed on zero of its own markers.
        observed_markers = {}  # column index -> set of markers
        unique_markers = {}    # column index -> set of markers resolving solely to it
        unique_marker_coverage = {}  # column index -> {marker: coverage} for those markers

        def observe(column_index, marker):
            observed_markers.setdefault(column_index, set()).add(marker)

        # Pass 1: create every column and record each OTU's resolved targets.
        otu_targets = []  # (marker, species_column_indices, clade_keys, coverage)
        for otu in sample_otus:
            best_hits = otu.equal_best_hit_taxonomies()
            if best_hits is None or len(best_hits) == 0:
                continue
            species_cols = set()
            clade_keys = set()
            for taxonomy in best_hits:
                if _is_species_string(taxonomy):
                    idx = species_column(taxonomy)
                    species_cols.add(idx)
                    observe(idx, otu.marker)
                else:
                    idx = novel_column(taxonomy)
                    observe(idx, otu.marker)
                    clade_keys.add(_canonical_species_key(taxonomy))
            # A window whose best hits tie across several species is not evidence
            # that one of those species is present: it resolves the read only to
            # their common ancestor. A novel organism in that clade produces
            # exactly this signature -- a tie spanning every DB species in the
            # clade -- and without a column for it the read must be attributed to
            # one of the known species (or, if none is admissible, dropped), so
            # novelty in a clade whose species are all in the DB is unrepresentable.
            # Load such windows onto the clade's novel column too, alongside the
            # tied species, and let the other evidence choose between them.
            # The novel column is added to this window's candidates directly rather
            # than as a clade key: a clade key expands to every species nested in
            # the clade, but a species-level tie already names the species the read
            # is compatible with, and the untied ones are positive evidence against.
            tie_column = None
            if len(species_cols) > 1 and len(clade_keys) == 0:
                lca = self._lowest_common_ancestor(best_hits)
                if lca is not None:
                    tie_column = novel_column(lca)
                    observe(tie_column, otu.marker)

            if len(species_cols) > 0 or len(clade_keys) > 0:
                otu_targets.append((otu.marker, species_cols, clade_keys, tie_column, otu.coverage))
            # A marker that resolves to a single taxon (one species, or one clade)
            # is unique evidence for it; ambiguous (shared) markers are not. Its
            # coverage is recorded too: these are the only observations that speak
            # about this taxon alone, so they are what the coherence ceiling is built
            # from.
            #
            # A species-level tie counts as unique evidence for the clade it ties
            # within. "Unique" must mean "nothing deeper than this taxon explains the
            # read", not "no other column appears in the row" -- and a window tied
            # across a genus's DB species is exactly the former: it resolves the read
            # to the genus and no further, which is what a clade-level hit is. The two
            # are the same observation reported differently, and which one SingleM emits
            # is an accident of how well the clade is sampled in the database. Withheld,
            # a novel organism whose windows always tie -- the common case in a genus
            # whose species are all in the DB, which is precisely where novelty is
            # hardest to see -- earns zero markers, is fixed to zero by the
            # identifiability floor before the optimiser runs, and its coverage vanishes
            # from the profile entirely rather than being misassigned.
            sole = None
            if len(species_cols) == 1 and len(clade_keys) == 0:
                sole = next(iter(species_cols))
            elif len(clade_keys) == 1 and len(species_cols) == 0:
                sole = novel_key_to_index[next(iter(clade_keys))]
            elif tie_column is not None:
                sole = tie_column
            if sole is not None:
                unique_markers.setdefault(sole, set()).add(otu.marker)
                coverages = unique_marker_coverage.setdefault(sole, {})
                coverages[otu.marker] = coverages.get(otu.marker, 0.0) + otu.coverage

        # Ensure every sylph-reported species has a column, even sylph-only ones. Sorted,
        # because column order must not depend on dict iteration order: the least-squares
        # problem has degenerate directions (a shared row cannot say which of its columns
        # owns the coverage), so a permutation of the columns moves where the optimiser
        # lands within them, and runs would not reproduce.
        for key in sorted(sylph_hits):
            hit = sylph_hits[key]
            if _is_species_string(hit.taxonomy):
                species_column(hit.taxonomy)

        # Map each clade key to the species columns nested within it.
        clade_to_species = {}
        for column_index, column in enumerate(columns):
            if column.kind != 'species':
                continue
            ranks = column.key.split(';')
            for depth in range(1, len(ranks)):
                clade_to_species.setdefault(';'.join(ranks[:depth]), []).append(column_index)

        # Pass 2: expand clade targets to novel leaf + nested species, and sum
        # the coverage of windows on the same marker resolving to the same
        # candidate set. Multiple windows of one marker (a genome's primary
        # window plus low-coverage error/fragment variants) are one measurement
        # of that marker's depth, not independent observations that would drag
        # the least-squares estimate down.
        aggregated = {}  # (marker, frozenset(cols)) -> summed coverage
        # Every window carrying a novel column is instead aggregated by (marker, clade),
        # unioning their candidate species. A clade may hold several novel organisms at
        # once, and each produces its own window on every marker -- a tie-set, a
        # clade-level hit, or a fragment. Left as separate rows they all assert the same
        # novel column, and least squares fits that column to the *mean* of the strains'
        # coverages when what the clade holds is their sum: a marker where twelve novel
        # strains each show 85x reads as twelve rows of 85x, so the column settles near
        # 85 instead of 1013, and the low fragment rows (1.6x) drag it down further,
        # since they too claim to measure the same thing. Merging restores additivity --
        # one row per marker, measuring the clade's whole novel coverage there, which is
        # exactly the per-marker total the trimmed-mean condense works from. A species
        # tied on this marker stays a candidate: if its window is one of those merged,
        # its coverage is part of the total the row observed.
        #
        # Markers on which each novel column's clade was actually resolved: some window
        # on that marker was tied across the clade, or hit it directly.
        clade_resolved_markers = {}
        for marker, species_cols, clade_keys, tie_column, coverage in otu_targets:
            if tie_column is not None:
                clade_resolved_markers.setdefault(tie_column, set()).add(marker)
            for clade_key in clade_keys:
                clade_resolved_markers.setdefault(novel_key_to_index[clade_key], set()).add(marker)

        clade_rows = {}  # (marker, novel column) -> [summed coverage, candidate columns]
        for marker, species_cols, clade_keys, tie_column, coverage in otu_targets:
            hidden = None
            novel_cols = set()
            if tie_column is not None:
                novel_cols.add(tie_column)
            else:
                hidden = self._hidden_novel_column(
                    species_cols, clade_keys, columns, novel_key_to_index,
                    marker, clade_resolved_markers)
                if hidden is not None:
                    novel_cols.add(hidden)
            for clade_key in clade_keys:
                novel_cols.add(novel_key_to_index[clade_key])

            if len(novel_cols) == 0:
                if len(species_cols) == 0:
                    continue
                key = (marker, frozenset(species_cols))
                aggregated[key] = aggregated.get(key, 0.0) + coverage
                continue

            cols = set(novel_cols)
            if hidden is not None:
                # This window is the species' own, and the coverage it observed is the
                # species' plus whatever novel organism shares its window here.
                cols |= species_cols
            else:
                # An ambiguous window: a tie across species, or a hit resolving only to
                # the clade. Its candidates are the clade's species -- except any whose
                # own window was directly observed on this marker. One organism has one
                # window per marker, so a species already seen in its own right on this
                # marker is not also in this marker's ambiguous coverage; its reads are
                # spent. Leaving it a candidate lets the model pay for the ambiguous
                # coverage twice over with the same species, which starves the novel
                # column of precisely the reads that make it novel.
                candidates = set(species_cols)
                for clade_key in clade_keys:
                    candidates.update(clade_to_species.get(clade_key, []))
                cols |= set(c for c in candidates if marker not in unique_markers.get(c, ()))
            # Merge on the deepest clade the window reached: that is the one whose
            # novel coverage this marker is measuring.
            deepest = max(novel_cols, key=lambda c: (len(columns[c].key.split(';')), c))
            entry = clade_rows.setdefault((marker, deepest), [0.0, set()])
            entry[0] += coverage
            entry[1].update(cols)

        # Each row carries its marker: a column's coverage is only credible if its
        # markers agree with each other, which cannot be checked without knowing
        # which marker each observation came from.
        singlem_rows = [(sorted(cols), coverage, marker) for (marker, cols), coverage in aggregated.items()]
        singlem_rows.extend(
            (sorted(cols), coverage, marker) for (marker, _), (coverage, cols) in clade_rows.items())

        observed_marker_count = [len(observed_markers.get(i, ())) for i in range(len(columns))]
        unique_marker_count = [len(unique_markers.get(i, ())) for i in range(len(columns))]
        return columns, singlem_rows, observed_marker_count, unique_marker_count, unique_marker_coverage

    def _column_marker_centres(self, columns, unique_marker_coverage, num_columns, min_markers):
        '''For each species column, the median coverage of the markers that resolve to
        it and nothing else -- the ceiling the coherence term holds it to. np.inf (no
        ceiling) for a column with fewer than min_markers such markers.

        Only unambiguous markers count. The obvious statistic -- apportion every row
        across its columns by their fitted shares and take each column's median -- is
        worthless here, because the apportionment is derived from the fit and so
        reproduces whatever the fit already believes: a column fitted at 408x is
        handed 37% of every shared row and duly "confirms" 408x. A shared row cannot
        testify about which of its columns the coverage belongs to; that is the whole
        reason it is shared. Only the markers where a taxon stands alone say what its
        coverage is, and being independent of the fit, they cannot be argued with.

        Novel columns get no ceiling, because the statistic does not mean anything for
        them. A species is one organism, so its markers agree and their median is its
        coverage. A novel column is the sum of however many novel organisms sit in the
        clade, and on any given marker only some of them resolve to the clade -- the
        rest are tied with species, or hidden inside a species row. Its unambiguous
        markers are therefore a mixture, not a measurement: for the Streptomyces of
        known50 they run 1094, 1003, 909, then 21, 19, 3.3, 1.6, and the median of 18
        is not a ceiling but a floor on how wrong a ceiling would be. Capping the novel
        column there would penalise it far harder than the species it competes with,
        and so drive the clade's coverage onto exactly the species we are trying to
        stop it landing on. The novel column is instead pinned by the rows, and by the
        marker floor that decides whether it may exist at all.'''
        centres = np.full(num_columns, np.inf)
        for c, markers in unique_marker_coverage.items():
            if columns[c].kind == 'species' and len(markers) >= min_markers:
                centres[c] = float(np.median(list(markers.values())))
        return centres

    def _hidden_novel_column(self, species_cols, clade_keys, columns, novel_key_to_index,
                             marker, clade_resolved_markers, min_ranks=6):
        '''The novel column a clade's novel organism must be hiding in on this row, or
        None.

        A window whose best hit is one species unambiguously is not proof that the
        species is what was sequenced. A novel member of the same genus has its own
        window at this marker, and if that window happens to coincide with the DB
        species' window -- these are short, and congeners often share them -- its reads
        land in the species' row with nothing to distinguish them. The species column is
        then the only thing in the row that can explain the coverage, so the novel
        organism's reads are booked to the species however much other evidence protests.

        The condition is per marker: on a marker where some window *did* resolve to the
        clade, we know where the novel organism's reads went (that row), and adding its
        column here as well would double-count it within the marker -- a novel organism
        contributes to exactly one row per marker, the one its window fell in. It is on
        the markers where the clade was *not* resolved that the organism is invisible
        and must be hiding inside a species row. Adding it there costs nothing when the
        clade holds no novel organism, because the novel column is then fixed to zero by
        the marker floor, which counts only clade-resolved markers.'''
        if len(clade_keys) > 0 or len(species_cols) == 0:
            return None  # already carries a novel column for its clade
        lineages = [columns[c].key.split(';') for c in species_cols]
        genus = []
        for ranks in zip(*lineages):
            if len(set(ranks)) != 1 or len(genus) == min_ranks:
                break
            genus.append(ranks[0])
        if len(genus) < min_ranks:
            return None  # the row's species do not share a genus
        column_index = novel_key_to_index.get(';'.join(genus))
        if column_index is None:
            return None  # no novel column for this clade at all
        if marker in clade_resolved_markers.get(column_index, ()):
            return None  # the clade resolved on this marker; its reads are elsewhere
        return column_index

    def _lowest_common_ancestor(self, taxonomies, min_ranks=6):
        '''The canonical key of the deepest clade containing every taxonomy given, or
        None if they share no rank, if the shared lineage is the whole thing, or if
        the clade is shallower than min_ranks (6 = genus, counting from domain).

        The rank floor is what keeps a novel column from becoming a sink for
        unresolvable coverage. Species tied within one genus look exactly like a
        novel member of that genus, but a tie spanning several genera is not evidence
        of a novel organism in the family -- it is a read that says little about
        anything, and a novel-at-family column would happily absorb whatever the
        deeper columns could not explain. Reads that genuinely resolve only to a high
        rank still reach these clades by the clade-hit path, which requires SingleM to
        have assigned them there.'''
        lineages = [_canonical_species_key(t).split(';') for t in taxonomies]
        shared = []
        for ranks in zip(*lineages):
            if len(set(ranks)) != 1:
                break
            shared.append(ranks[0])
        if len(shared) < min_ranks or len(shared) == min(len(lineage) for lineage in lineages):
            return None
        return ';'.join(shared)

    def _bisquare_weights(self, residuals, observed, c=4.685, coverage_floor=1.0, min_weight=1e-3):
        '''Bisquare IRLS weights, judged on each row's residual as a fraction of the
        coverage that row observed rather than on the raw residual.

        Marker coverage noise scales with coverage: a 10x residual on a 30x taxon is
        ordinary, while the same residual on a 2x taxon is not. A robust scale taken
        from raw residuals is set by the many low-coverage rows, so every row of the
        most abundant taxon in the sample exceeds it and is rejected as an outlier --
        and since bisquare rejection is absorbing (zero weight means zero gradient,
        so the fit never rises back toward the row), the dominant taxon collapses to
        near zero. Scaling the residual first makes the rejection rule size-blind: a
        row is an outlier only if the fit explains an unusual *fraction* of it.

        The scale divisor is the observed coverage, not the fitted one, so that it
        does not depend on the current iterate. Dividing by the fit re-creates the
        absorbing state by another route: a row whose column is currently fitted near
        zero has an enormous relative residual however ordinary it is, so it is
        rejected, so the column can never rise. Dividing by the observation instead
        means a row is judged the same way no matter where the fit currently sits,
        and since the first pass is unweighted, rejection is only ever decided
        against an already-fitted solution.

        The loss itself stays in absolute coverage units, matching the sylph and
        padding terms; only which rows are trusted is decided on the relative scale.
        min_weight keeps a rejected row's gradient alive so a fit can climb back.'''
        denominator = np.maximum(observed, coverage_floor)
        relative = residuals / denominator
        scale = 1.4826 * np.median(np.abs(relative))
        if not np.isfinite(scale) or scale <= 0:
            return np.ones_like(residuals)
        u = relative / (c * scale)
        weights = (1.0 - u * u) ** 2
        weights[np.abs(u) >= 1.0] = 0.0
        return np.maximum(weights, min_weight)

    def _build_profile(self, sample, columns, a, epsilon=1e-6):
        root = WordNode(None, 'Root')
        for column, coverage in zip(columns, a):
            if coverage <= epsilon:
                continue
            root.add_words(_gtdb_string_to_wordnode_array(column.taxonomy), float(coverage))
        return CondensedCommunityProfile(sample, root)


class _Column:
    __slots__ = ('kind', 'taxonomy', 'key')

    def __init__(self, kind, taxonomy, key):
        self.kind = kind          # 'species' or 'novel'
        self.taxonomy = taxonomy  # original taxonomy string (for tree building)
        self.key = key            # canonical species/clade key
