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
              l1_penalty=1.0, absence_weight=100.0, sylph_weight=1.0, min_markers=3,
              max_outer_iterations=25, tolerance=1e-4, prune_below=0.05):
        try:
            from scipy.optimize import minimize
            from scipy.sparse import csr_matrix
        except ImportError:
            raise Exception("condense --joint requires scipy, which could not be imported")

        columns, singlem_rows, observed_marker_count, unique_marker_count = \
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
        for r, (cols, coverage) in enumerate(singlem_rows):
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

        if alpha is None:
            current_alpha = 1.0
            fit_alpha = len(sylph_col) > 0
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
        sylph_columns = set(int(c) for c in sylph_col)
        bounds = [(0.0, None)] * num_columns
        num_fixed = 0
        if min_markers > 0:
            for j, column in enumerate(columns):
                if column.kind == 'species' and j not in sylph_columns and unique_marker_count[j] < min_markers:
                    bounds[j] = (0.0, 0.0)
                    num_fixed += 1
        if num_fixed > 0:
            logging.info("Fixed {} sylph-unsupported columns to zero (< {} unique markers)".format(
                num_fixed, min_markers))

        a = np.zeros(num_columns)
        singlem_weights = np.ones(num_singlem)

        for outer in range(max_outer_iterations):
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
                f += l1_penalty * float(np.sum(x))
                g += l1_penalty
                return f, g

            result = minimize(objective, a, jac=True, method='L-BFGS-B', bounds=bounds)
            new_a = np.maximum(result.x, 0.0)

            # Variable-projection alpha update.
            new_alpha = current_alpha
            if fit_alpha and len(sylph_col) > 0:
                a_sylph = new_a[sylph_col]
                denominator = float(np.dot(sylph_w, a_sylph * a_sylph))
                if denominator > 0:
                    new_alpha = float(np.dot(sylph_w, sylph_eff * a_sylph) / denominator)

            # IRLS reweighting of SingleM rows (bisquare on residuals).
            if num_singlem > 0:
                singlem_weights = self._bisquare_weights(b_singlem - M.dot(new_a))

            converged = (np.max(np.abs(new_a - a)) < tolerance and abs(new_alpha - current_alpha) < tolerance)
            a, current_alpha = new_a, new_alpha
            if converged:
                logging.debug("Joint deconvolution converged after {} iterations".format(outer + 1))
                break

        logging.info("Joint deconvolution of sample {}: alpha={:.4f}, total coverage={:.2f}".format(
            sample, current_alpha, float(np.sum(a))))

        # L-BFGS-B only approximates the L1 kink at zero, leaving tiny non-zero
        # values; clear them so the profile is genuinely sparse.
        a[a < prune_below] = 0.0

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
            if len(species_cols) > 0 or len(clade_keys) > 0:
                otu_targets.append((otu.marker, species_cols, clade_keys, otu.coverage))
            # A marker that resolves to a single taxon (one species, or one clade)
            # is unique evidence for it; ambiguous (shared) markers are not.
            if len(species_cols) == 1 and len(clade_keys) == 0:
                unique_markers.setdefault(next(iter(species_cols)), set()).add(otu.marker)
            elif len(clade_keys) == 1 and len(species_cols) == 0:
                unique_markers.setdefault(novel_key_to_index[next(iter(clade_keys))], set()).add(otu.marker)

        # Ensure every sylph-reported species has a column, even sylph-only ones.
        for hit in sylph_hits.values():
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
        for marker, species_cols, clade_keys, coverage in otu_targets:
            cols = set(species_cols)
            for clade_key in clade_keys:
                cols.add(novel_key_to_index[clade_key])
                cols.update(clade_to_species.get(clade_key, []))
            if len(cols) == 0:
                continue
            key = (marker, frozenset(cols))
            aggregated[key] = aggregated.get(key, 0.0) + coverage

        singlem_rows = [(sorted(cols), coverage) for (marker, cols), coverage in aggregated.items()]

        observed_marker_count = [len(observed_markers.get(i, ())) for i in range(len(columns))]
        unique_marker_count = [len(unique_markers.get(i, ())) for i in range(len(columns))]
        return columns, singlem_rows, observed_marker_count, unique_marker_count

    def _bisquare_weights(self, residuals, c=4.685):
        absolute = np.abs(residuals)
        scale = 1.4826 * np.median(absolute)
        if not np.isfinite(scale) or scale <= 0:
            return np.ones_like(residuals)
        u = residuals / (c * scale)
        weights = (1.0 - u * u) ** 2
        weights[np.abs(u) >= 1.0] = 0.0
        return weights

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
