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
  # Caches written / read by 02–03
  gpmap_indications.rda
  gpmap_indications_rarevariants.rda
  gpmap_ticolocs.rda
  gpmap_tisharedrare.rda
  gpmap_tipairwisecolocs.rda
  tipairs_preinfomap.rda                     # pre-infomap T-I support (from 03)
  exclude_studies.txt                        # optional
```

## Workflow

### T–I pair scoring / truth set

| Script | Cadence | Purpose |
|--------|---------|---------|
| `01_extracting_ti_pairs.R` | Once (or when trait–MeSH matching changes) | Build `target-indicationpairs_gpmapevidence.tsv` |
| `02_gpmap_support_for_ti_pairs.Rmd` | Re-run when clusters / results version change | Score GC support; write `GPMAP_T-Ipairs_*.tsv`. Param `recompute_cache` controls cache load vs overwrite. |
| `03_ti_pair_validation_of_coloc.Rmd` | Re-run per clustering comparison | Compare pairwise / pre-infomap / final clusters. Param `recompute_cache` controls `tipairs_preinfomap.rda`. |
| `04_pull_truth_links.R` | After 02+03 | Write `tipairs_launched_truthset.tsv` |

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

# Load cached intermediates (default)
Rscript -e 'rmarkdown::render("02_gpmap_support_for_ti_pairs.Rmd",
  params = list(results_version = "1.0.0", recompute_cache = FALSE),
  output_dir = "/local-scratch/projects/genotype-phenotype-map/results/1.0.0/analysis/clustering")'

# Recompute and overwrite caches under $DATA_DIR/ti_pairs/
Rscript -e 'rmarkdown::render("02_gpmap_support_for_ti_pairs.Rmd",
  params = list(results_version = "1.0.0", recompute_cache = TRUE),
  output_dir = "/local-scratch/projects/genotype-phenotype-map/results/1.0.0/analysis/clustering")'

Rscript -e 'rmarkdown::render("03_ti_pair_validation_of_coloc.Rmd",
  params = list(results_version = "1.0.0", recompute_cache = FALSE))'

Rscript 04_pull_truth_links.R --results_version 1.0.0

Rscript -e 'rmarkdown::render("coloc_clustering_validation.Rmd", params = list(results_version = "1.0.0"))'
Rscript -e 'rmarkdown::render("clustering_post_analysis.Rmd", params = list(results_version = "1.0.0"))'

Rscript rank_truthset_ld_blocks.R --results_version 1.0.0 --top_n 20
```
