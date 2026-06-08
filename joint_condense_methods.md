## Methods

### Integration of genome-resolved abundances into marker-gene taxonomic profiles

SingleM derives a taxonomic profile of a metagenome from single-copy marker genes: short "windows" of each marker are recovered from reads, assigned to taxa, and condensed into per-taxon coverages (in genome-equivalent units) by a sequence of expectation–maximisation (EM) steps and a trimmed-mean summarisation across markers. This procedure quantifies the total community, including lineages with no database representative ("dark matter"), but is comparatively insensitive at low abundance and cannot resolve ambiguity between database species that share marker-window sequences. Sylph, by contrast, estimates genome-resolved coverages by k-mer containment against a reference database and is both more sensitive and more precise for species that have a reference genome, but is blind to novel lineages.

We integrate the two by supplying SingleM's `condense` step with a sylph profile of the same sample. Both tools are run against databases built from the same Genome Taxonomy Database (GTDB) release, so their taxon strings are identical and every species SingleM can detect is one sylph can also quantify; this removes any cross-database taxon mapping and reduces the integration to a single question — how much marker support SingleM has for each taxon. We implemented two integration strategies, a lightweight additive injection and a joint deconvolution, selected at run time.

### Sylph profile ingestion and taxonomic harmonisation

Sylph effective coverages are provided as a tab-separated table annotated with GTDB taxonomy (genome accessions mapped to GTDB species via the release's taxonomy files). Each taxon string is reduced to a canonical key by removing the synthetic root rank and any empty ranks and concatenating the remainder, so that sylph and SingleM taxa can be matched exactly. As a guard against database-release mismatch, the fraction of sylph species absent from the metapackage's species set is computed and a warning is emitted if it exceeds one half.

### Coverage-scale calibration (alpha)

SingleM coverage (mean per-base depth over single-copy marker regions) and sylph effective coverage (`Eff_cov`) are commensurate only up to a per-sample scale factor $\alpha$, reflecting differences in how each statistic responds to read length, sequencing error and k-mer subsampling. Where a fixed $\alpha$ is required, it is estimated by regressing SingleM coverage on sylph effective coverage through the origin,

$$
\alpha \;=\; \frac{\sum_{s\in\mathcal{A}} e_s\, c_s}{\sum_{s\in\mathcal{A}} e_s^{2}},
$$

over the anchor set $\mathcal{A}$ of species detected by both tools at SingleM coverage $c_s \ge 10\times$ (with $e_s$ the sylph effective coverage). When fewer than three such anchor species are available the regression is considered under-powered and $\alpha$ defaults to 1. A user-supplied value overrides the estimate.

### Additive injection of sylph-only species

The first strategy retains SingleM's profile unchanged and augments it with species that sylph detects but SingleM misses entirely. After the standard condense tree has been built, every sylph species absent from the SingleM profile is added as a new species leaf with coverage $\alpha\,e_s$. Because these injections occur after the EM, they bypass SingleM's proximity-based pruning. To keep SingleM authoritative over each clade's total, the injected coverage is reconciled against the existing profile: it is drawn down from the nearest ancestral internal node's unresolved ("novel") coverage, descending toward the root, and only coverage in excess of that budget is added to the community total. Injection therefore converts unresolved clade-level coverage into concrete species where sylph evidence permits, and adds genuinely new coverage only where sylph exceeds SingleM's clade budget, avoiding double counting.

### Joint deconvolution of SingleM and sylph

The second strategy replaces SingleM's EM and trimmed-mean summarisation with a single penalised, weighted, non-negative least-squares (NNLS) deconvolution per sample that explains the marker and sylph observations jointly.

**Variables.** We solve for a non-negative coverage vector $a$ over an augmented taxonomy: one entry per *candidate species* (the union of all species appearing in SingleM equal-best-hit sets and all sylph-reported species) and one *novel-at-clade* entry for each internal clade to which some marker signal resolves only partially. Novel entries make the unresolved fraction of each clade an explicit variable rather than a residual.

**SingleM observation block.** Each marker contributes one observation per distinct candidate set. For a marker $m$ and a set of candidate columns $\mathcal{C}$, the coverages of all windows of $m$ whose equal-best set maps to $\mathcal{C}$ are summed into a single target $b_{(m,\mathcal{C})}$, with unit design coefficients on the columns in $\mathcal{C}$. Aggregating windows in this way treats a genome's coverage on a marker as one measurement of that marker's depth, rather than as several independent measurements (a marker's true window plus low-coverage error- or fragment-derived window variants), which would otherwise bias the least-squares estimate downward. Windows that resolve only to a clade (genus-level DIAMOND assignments, or query assignments whose equal-best set spans a clade) load onto that clade's novel-at-clade column *and* every candidate species nested within the clade; this allows clade-level reads of a known species to be attributed to that species when other evidence supports it, rather than being counted separately as novelty.

**Marker-count normalisation.** SingleM coverage is an average over a taxon's full single-copy marker complement, so taxa observed on few markers must be down-weighted accordingly. For each column we add a penalty proportional to the number of *unobserved* markers, $w^{\text{pad}}_j = \max(0,\, N_{d(j)} - k_j)$, where $k_j$ is the number of distinct markers on which column $j$ was directly observed and $N_{d(j)}$ is the number of markers targeting its domain $d(j)$. This reproduces the zero-padding of the trimmed-mean summarisation and suppresses spurious taxa supported by only a handful of mis-assigned markers.

**Sylph observation block and absence constraint.** Each sylph-reported species $s$ contributes a row $e_s = \alpha\, a_s$. Because sylph is the more sensitive instrument, a database species that sylph does *not* report is treated as strong evidence of absence: each such candidate species receives a high-weight zero-target row, with weight $w^{\text{abs}}$ (default 100), that drives its coverage toward zero. Unexplained marker signal that could belong to an absent species is thereby reinterpreted as clade-level novelty rather than attributed to that species.

**Identifiability floor.** A species that sylph does not detect must, in addition, be uniquely resolved by at least a minimum number of markers (default three) to be retained; otherwise its coverage is constrained to zero. A marker counts toward this total only when it resolves to that species alone, mirroring the unique-marker criterion SingleM uses internally to whitelist confidently-detected species. This removes false positives that would otherwise persist on a single uniquely-assigned marker amid windows shared with a sylph-supported near neighbour: the shared coverage re-routes to the supported species and the lone discordant marker is treated as misassignment. Sylph-detected species (including sylph-only species, for which sylph itself is the evidence) and novel-at-clade entries (which represent dark matter and are governed only by the sparsity penalty) are exempt.

**Objective.** Collecting the SingleM rows into a sparse matrix $M$ with targets $b$ and diagonal weights $W$, the per-sample estimate is

$$
\hat a \;=\; \arg\min_{a \ge 0}\;
\big\| W^{1/2}(b - M a)\big\|_2^{2}
\;+\; \sum_{s} w^{\text{syl}}_s\,(e_s - \alpha\,a_s)^2
\;+\; \sum_{s\notin\text{sylph}} w^{\text{abs}} a_s^{2}
\;+\; \sum_j w^{\text{pad}}_j a_j^{2}
\;+\; \lambda \sum_j a_j .
$$

The $\ell_1$ term (penalty $\lambda$, default 1) induces sparsity, subsuming SingleM's proximity pruning and minimum-coverage cutoff; non-negativity is biologically required and further sparsifies the solution. Identifiability of collinear (window-sharing) species columns, which motivates SingleM's EM, is restored by the sylph rows, which add a near-orthogonal constraint for every candidate species under the shared-database assumption.

**Robustness.** To guard against per-marker coverage outliers (paralogy, horizontal transfer, contamination) that the trimmed mean otherwise removes, the SingleM rows are reweighted by iteratively reweighted least squares: after each solve, weights are recomputed from the residuals $r$ with Tukey's bisquare function, $w_i = (1 - (r_i / c\hat\sigma)^2)^2$ for $|r_i| < c\hat\sigma$ and $0$ otherwise, with tuning constant $c = 4.685$ and scale $\hat\sigma = 1.4826\,\operatorname{median}(|r|)$.

**Scale estimation by variable projection.** Rather than pre-fitting $\alpha$, it is estimated jointly with $a$ by variable projection: each outer iteration solves the NNLS for $a$ at the current $\alpha$, then updates $\alpha$ in closed form, $\alpha = \sum_s w^{\text{syl}}_s e_s a_s \big/ \sum_s w^{\text{syl}}_s a_s^2$, over sylph-reported species. The shared database makes the jointly-detected species set large, so $\alpha$ is well determined. A user-supplied $\alpha$ fixes the value and disables this step.

**Optimisation and output.** Each weighted, $\ell_1$-penalised NNLS subproblem is solved with the bounded limited-memory BFGS algorithm (L-BFGS-B) using analytic gradients; the outer reweighting/variable-projection loop iterates to convergence (default tolerance $10^{-4}$, maximum 25 iterations). Because L-BFGS-B only approximates the $\ell_1$ kink at the origin, coverages below a small floor (default 0.05) are set to zero. The solution is written into the SingleM coverage tree, with species columns placed at their leaves and novel-at-clade columns placed as the corresponding internal node's own (unresolved) coverage, yielding a profile of the same form as standard `condense`.

Viral profiling, in which markers are not single-copy and design coefficients are not binary, is not supported by the joint mode and falls back to the standard algorithm.

### Implementation

Both strategies are implemented in SingleM's `condense` module (Python); the joint deconvolution uses NumPy and SciPy (sparse linear algebra and `scipy.optimize`). The additive-injection mode is enabled by supplying a sylph profile; the joint mode is enabled additionally with a command-line flag, with the $\ell_1$ penalty, absence weight and $\alpha$ exposed as parameters.

### In silico validation

A mock metagenome was simulated to evaluate the integration. Paired-end reads (2 × 150 bp) were generated with wgsim at a reduced base-error rate (0.1%, to approximate modern Illumina chemistry) from two *Methanobacteriaceae* genomes drawn from GTDB r232 — one at 10× and one at 0.5× coverage — and combined. Sylph (`-c 200`) was run against a GTDB r232 database and its output annotated with GTDB taxonomy; SingleM `pipe` was run against the corresponding GTDB r232 metapackage to produce an archive OTU table; both were passed to `condense`. The high-coverage genome was recovered by both tools, whereas the low-coverage genome fell below SingleM's marker-detection sensitivity and was recovered only via sylph, providing a direct test of both the additive-injection and joint paths. The complete workflow (read simulation, sylph and SingleM execution, and taxonomy annotation) was implemented in Snakemake to support reproduction.
