# Scientific Methodology Roadmap

> Status: proposed methodology for technical discussion. The items below are
> not claimed as implemented or experimentally validated by the current code.

## What should be improved

1. **Separate representation from biological confidence.** `bin` should remain
   the presence/absence baseline. For `cnt`, the primary candidate should be
   sublinear term frequency, `tf = log(1 + count)`, followed by an optional
   inverse-frequency factor,
   `idf = log((N + 1) / (df + 1)) + 1`, and L1 or L2 normalization per complex.
   This prevents repeated contacts and larger ligands from dominating the
   representation. Raw count, `log1p`, TF-IDF, and binary fingerprints should
   be compared under the same external validation protocol. None should be
   called a biological confidence score until it has been calibrated against
   independent evidence. Count-aware similarity must use weighted Tanimoto or
   cosine distance rather than silently coercing counts to binary values.

2. **Replace the blind level matrix with an information-gain stopping rule.**
   At every expansion level, the pipeline should record normalized Shannon
   entropy of feature frequencies, Jensen-Shannon divergence from the previous
   level, new-node and new-edge fractions, graph density, largest-component
   fraction, hub dominance, and community stability. A candidate stopping rule
   may require two consecutive levels with marginal divergence below a
   preregistered threshold and increasing density or hub dominance. Thresholds
   must be fitted on development data and frozen before final evaluation; they
   must not be selected on the test set. This is a project-specific decision
   rule, not an established biological law.

3. **Control ascertainment and study bias.** Every edge should retain database,
   assay, publication, species, interaction type, and evidence count. Results
   should be repeated with curated-only, high-throughput-only, physical-only,
   and combined evidence strata. Node rankings should be compared with and
   without publication/evidence covariates, while topology statistics should be
   standardized against degree-preserving nulls. Published analyses have shown
   substantial literature bias in curated PPI layers, and topology-driven
   prediction can be biased by hubs unless degree is normalized
   ([Luck et al.-associated multilayer analysis, 2021](https://doi.org/10.1038/s41467-021-26674-1),
   [Kovács et al., 2019](https://doi.org/10.1038/s41467-019-09177-y)).

4. **Quantify uncertainty and stability.** Results should include confidence
   intervals obtained by bootstrap over complexes/ligands, sensitivity to
   fingerprint size, and rank stability across seeds, level/growth settings,
   evidence filters, and database releases. Splits must be made by biological
   entity or scaffold when ordinary random splitting would leak highly related
   examples across training and test sets.

## What should be added

1. **Degree-preserving null models.** For each observed network, at least 1,000
   accepted double-edge-swap randomizations should be generated after excluding
   self-loops and duplicate edges. Degree, centrality, modularity, path-length,
   and motif statistics should be reported as
   `z = (observed - mean(null)) / sd(null)`, together with the empirical
   `p = (1 + number(null >= observed)) / (B + 1)` and false-discovery-rate
   correction. Maslov and Sneppen established the degree-preserving comparison
   for molecular networks ([Science, 2002](https://pubmed.ncbi.nlm.nih.gov/11988575/)).
   Weighted networks additionally require a strength-aware null or a separate
   permutation of weights; binary rewiring alone does not test weighted claims.

2. **A compact, nonredundant topology panel.** Export degree/strength,
   betweenness, closeness or harmonic centrality for disconnected graphs,
   eigenvector centrality, PageRank, k-core/coreness, clustering coefficient,
   component membership, and null-standardized scores. Personalized PageRank
   should be used only when biologically justified seed nodes are declared.
   Computationally expensive metrics should expose exact and approximate modes
   for HPC-scale graphs.

3. **Robust community analysis.** Leiden should be the default modularity/CPM
   optimizer because it addresses poorly connected communities produced by
   Louvain ([Traag et al., 2019](https://doi.org/10.1038/s41598-019-41695-z)).
   Infomap supplies a complementary flow-based partition
   ([Rosvall and Bergstrom, 2008](https://pubmed.ncbi.nlm.nih.gov/18216267/)).
   Multiple random seeds and resolution values should be evaluated, with
   variation of information or adjusted mutual information used to quantify
   partition stability. The original Louvain result remains a useful benchmark
   ([Blondel et al., 2008](https://doi.org/10.1088/1742-5468/2008/10/P10008)).

4. **Versioned external biological validation.** Identifiers must first be
   normalized to a stable namespace and species. STRING physical and functional
   networks must not be mixed; production requests should use a versioned API
   endpoint as recommended by the
   [official STRING API](https://string-db.org/help/api/). BioGRID validation
   should preserve experimental-system metadata and record the monthly database
   version through its [official REST service](https://wiki.thebiogrid.org/doku.php/biogridrest).
   IntAct evidence can be obtained in PSI-MITAB through the
   [official PSICQUIC service](https://www.ebi.ac.uk/intact/documentation/user-guide).
   Source-held-out or time-held-out evaluation is preferable to validating
   against a database that contributed to model construction.

5. **Evaluation without false negatives.** Unobserved protein pairs are not
   confirmed noninteractions. Sensitivity and precision can be reported against
   positives, but specificity and ROC require a declared negative reference
   set. Negative sampling by subcellular separation can itself inflate accuracy
   ([Ben-Hur and Noble, 2006](https://pmc.ncbi.nlm.nih.gov/articles/PMC1810313/)).
   When reliable negatives are unavailable, positive-unlabelled evaluation and
   several sensitivity analyses should replace a single arbitrary negative set.
   For strongly imbalanced test sets, PR-AUC should accompany ROC-AUC because
   precision-recall curves are more informative in that regime
   ([Davis and Goadrich, 2006](https://doi.org/10.1145/1143844.1143874)).

6. **Automated functional enrichment.** Both an over-representation analysis
   of a declared top set and a rank-based analysis of all nodes should be
   supported. The statistical universe must be the proteins that could have
   entered the experiment, not automatically the whole proteome. GO's official
   guidance explicitly defines enrichment relative to a background
   ([Gene Ontology](https://www.geneontology.org/docs/go-enrichment-analysis/)).
   Reactome's Analysis Service returns over-representation p-values and FDR
   values and provides a programmatic endpoint
   ([Reactome](https://reactome.org/dev/analysis)). Outputs should record
   database release, identifier mapping coverage, background, test, effect size,
   raw p-value, Benjamini-Hochberg FDR, and a redundancy-reduced term table.

## Impact versus effort

| Priority | Proposal | Scientific impact | Coding effort | Recommended milestone |
|---|---|---:|---:|---|
| 1 | Preserve raw `bin`/`cnt`; add `log1p`, TF-IDF, normalization, and weighted similarity | High | Low-medium | Immediate |
| 2 | Add topology panel plus degree-preserving null distributions and empirical p-values | Very high | Medium | First scientific release |
| 3 | Add versioned provenance, identifier mapping, evidence strata, and database snapshots | Very high | Medium | First scientific release |
| 4 | Add Leiden/Infomap stability analysis | High | Medium | Second milestone |
| 5 | Add entropy/divergence diagnostics and preregistered automatic stopping | High | Medium | Second milestone |
| 6 | Add STRING/BioGRID/IntAct held-out validation with careful negative policy | Very high | High | Validation release |
| 7 | Add GO/Reactome over-representation and ranked enrichment | High | Medium-high | Validation release |

The highest impact-to-effort starting point is a reproducible representation
benchmark: retain `bin` and raw `cnt`, add `log1p` and TF-IDF variants, freeze
train/validation/test splits, and compare them using the same external evidence.
The next increment should couple the topology panel to degree-preserving nulls;
centrality without a null distribution or sensitivity analysis remains mainly
descriptive.
