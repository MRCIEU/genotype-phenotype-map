# Clustering analysis

Scripts for coloc clustering validation and target–indication (T–I) pair support /
truth-set construction used to evaluate clustering methods.

## Paths

From `pipeline_steps/constants.R`:

- `ti_pairs_data_dir` → `$DATA_DIR/ti_pairs/` (Minikel inputs, matching tables, caches)
- Versioned pipeline artefacts → `$RESULTS_DIR/<results_version>/`
- Analysis outputs → `$RESULTS_DIR/<results_version>/analysis/clustering/`

### Expected layout under `$DATA_DIR/ti_pairs/`

All inputs and caches are flat under `ti_pairs/`:

```text
ti_pairs/
  # Minikel et al. (https://zenodo.org/records/10783210)
  Supplement_ST01.csv
  assoc.tsv.gz
  sim.tsv.gz
  # Trait–MeSH matching (inputs + output of 01)
  gpmap_matched_traits.csv
  gpmap_matched_traits_filtered_corrected.csv
  target-indicationpairs_gpmapevidence.tsv   # written by 01
  # Scored T-I tables written by 02; read by 03–04
  GPMAP_T-Ipairs_allmatchedstudies.tsv
  GPMAP_T-Ipairs_unique.tsv
  # Caches written / read by 02–03 (shared across include_trans modes)
  gpmap_indications.rda
  gpmap_indications_rarevariants.rda
  gpmap_tipairwisecolocs.rda
  gpmap_indications_pairwisecolocs.rda
  coloc_pairwise_indications.txt
  gpmap_indications_ids.txt
  exclude_studies.txt                        # optional
  # Per-mode caches / outputs; include_trans = FALSE writes the *_notrans variants
  gpmap_ticolocs.rda
  gpmap_tisharedrare.rda
  tipairs_preinfomap.rda                     # pre-infomap T-I support (from 03)
  tipairs_launched_truthset.tsv              # truth set (from 04)
```

## Workflow

### T–I pair scoring / truth set

| Script | Cadence | Purpose |
|--------|---------|---------|
| `01_extracting_ti_pairs.R` | Once (or when trait–MeSH matching changes) | Build `target-indicationpairs_gpmapevidence.tsv` |
| `02_gpmap_support_for_ti_pairs.Rmd` | Re-run when clusters / results version change | Score GC support; write `GPMAP_T-Ipairs_*.tsv`. `recompute_cache` controls derived caches; `recompute_api_cache` controls the slow `gpmapr::trait()` pulls (defaults to `recompute_cache`); `include_trans = FALSE` drops trans QTL markers post-pull and writes `*_notrans` outputs. |
| `03_ti_pair_validation_of_coloc.Rmd` | Re-run per clustering comparison | Compare pairwise / pre-infomap / final clusters. `recompute_cache` controls `tipairs_preinfomap[<suffix>].rda`; `include_trans` must match 02. |
| `04_pull_truth_links.R` | After 02+03 | Write `tipairs_launched_truthset[<suffix>].tsv`. Pass `--include_trans FALSE` for the no-trans set. |

### Clustering method comparison

| Script | Purpose |
|--------|---------|
| `coloc_clustering_validation.Rmd` | Compare clustering parameterisations across LD blocks |
| `clustering_post_analysis.Rmd` | Mode comparison + truth-set recovery diagnostics |
| `rank_truthset_ld_blocks.R` | Rank LD blocks by truth-set gaps |
| `clustering_comparison_helpers.R` | Shared helpers for the Rmds above |

### Examples

```bash
cd scripts/analysis/clustering

Rscript 01_extracting_ti_pairs.R --results_version 1.0.0

# Mode A: include trans QTLs. Refreshes the shared API cache (recompute_api_cache = TRUE).
Rscript -e 'rmarkdown::render("02_gpmap_support_for_ti_pairs.Rmd",
  params = list(results_version = "1.0.0", recompute_cache = TRUE,
                recompute_api_cache = TRUE, include_trans = TRUE),
  output_file = "02_gpmap_support_for_ti_pairs.html",
  output_dir = "/local-scratch/projects/genotype-phenotype-map/results/1.0.0/analysis/clustering")'

# Mode B: exclude trans QTLs. Re-uses the shared API cache (recompute_api_cache = FALSE),
# recomputes downstream caches and writes *_notrans outputs.
Rscript -e 'rmarkdown::render("02_gpmap_support_for_ti_pairs.Rmd",
  params = list(results_version = "1.0.0", recompute_cache = TRUE,
                recompute_api_cache = FALSE, include_trans = FALSE),
  output_file = "02_gpmap_support_for_ti_pairs_notrans.html",
  output_dir = "/local-scratch/projects/genotype-phenotype-map/results/1.0.0/analysis/clustering")'

Rscript -e 'rmarkdown::render("03_ti_pair_validation_of_coloc.Rmd",
  params = list(results_version = "1.0.0", recompute_cache = TRUE, include_trans = TRUE),
  output_file = "03_ti_pair_validation_of_coloc.html")'

Rscript -e 'rmarkdown::render("03_ti_pair_validation_of_coloc.Rmd",
  params = list(results_version = "1.0.0", recompute_cache = TRUE, include_trans = FALSE),
  output_file = "03_ti_pair_validation_of_coloc_notrans.html")'

Rscript 04_pull_truth_links.R --results_version 1.0.0 --include_trans TRUE
Rscript 04_pull_truth_links.R --results_version 1.0.0 --include_trans FALSE

# Load cached intermediates (default)
Rscript -e 'rmarkdown::render("02_gpmap_support_for_ti_pairs.Rmd",
  params = list(results_version = "1.0.0", recompute_cache = FALSE, include_trans = TRUE),
  output_dir = "/local-scratch/projects/genotype-phenotype-map/results/1.0.0/analysis/clustering")'

Rscript -e 'rmarkdown::render("coloc_clustering_validation.Rmd", params = list(results_version = "1.0.0"))'
Rscript -e 'rmarkdown::render("clustering_post_analysis.Rmd", params = list(results_version = "1.0.0"))'

Rscript rank_truthset_ld_blocks.R --results_version 1.0.0 --top_n 20
```
