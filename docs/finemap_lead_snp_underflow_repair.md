# Finemap lead-SNP underflow tie-break bug: fix and rollout runbook

Date: 2026-09-30

## Summary

The lead SNP chosen per finemapped credible set was arbitrary whenever multiple
SNPs' LBF-based p-values underflowed to the same double (typically 0), because
the choice was made with `which.min(LBF_P)` and `LBF_P = 2 * pnorm(-z)`
underflows to exactly 0 for strong signals. With a 79-way tie at 0, `which.min`
returns the first row in file order (lowest BP), not the strongest SNP.

Confirmed on GWAS upload `c6d2e20f-6ece-d093-78ee-37534547c466`
(EUR/10/109285213-113568394, credible set 1): the TCF7L2 signal's lead was
recorded as `10:112974337_A_G` (arbitrary lowest-BP underflowed SNP) instead of
`10:112995025_C_T` (LBF 2438). The wrong `snp`/`bp`/`min_p` propagate into
`compiled_extracted_studies.tsv` / `study_extractions` / DBs / static web.

`min_p` itself is functionally unchanged by the fix: it was already the true
global minimum p-value; only the labelled lead SNP (and `cis_trans` when the
corrected lead crosses +/-1 Mb of the original signal bp) changes.

## Code fix (already applied, verified, lint-clean)

- `pipeline_steps/gwas_calculations.R`: new `convert_lbf_to_log_p_value()`
  (`log(2) + pnorm(-abs(z), log.p = TRUE)`) - log-scale p-values cannot
  underflow, so ranking never ties at 0. Note: `2 * pnorm(..., log.p = TRUE)`
  would compute `log(p^2)`, not `log(2p)` - do not "simplify" back to it.
- `pipeline_steps/finemap_studies_in_ld_block.R`
  `split_susie_result_into_conditional_gwases()`: `LBF_P` is now log-scale
  (`LBF_LOG_P`); the bad-SNP filter compares against
  `log(lowest_p_value_threshold)`; `which.min` is exact;
  `new_min_p <- exp(min(LBF_LOG_P, na.rm = T))` keeps `min_p` a valid p-value.
  `LBF_P` is dropped before files are written, so on-disk finemapped formats
  are unchanged.
- `pipeline_steps/finemap_studies_in_ld_block.R` (two sites, the
  `process_unfinemapped_gwas` / failed-susie paths): `which.min(gwas$P)` is
  replaced by `order(gwas$P, -abs(gwas$Z))[1]` so P-value ties/underflow are
  broken by the most extreme Z.
- `pipeline_steps/cluster_studies_in_ld_block.R` duplicate-removal:
  `which.min(min_p)` ties broken deterministically by
  `order(min_p, bp)[1]`.

## Scan results (2026-09-30, main-pipeline machine)

Produced by `scripts/find_finemap_lead_snp_inconsistencies.R` (run inside the
apptainer image; only rows with `min_p < 1e-300` are re-verified since only
subnormal/zero-range p-values can tie in double precision):

- `scripts/finemap_lead_snp_inconsistencies.tsv` - all 4,974 re-verified rows
  with old/new lead `snp`, `bp`, `min_p`, `cis_trans`, `lead_snp_changed`
- `scripts/finemap_lead_snp_affected_blocks.tsv` - affected study/blocks
  (study_name, ld_block, variant_type, n_changed_rows)
- `scripts/finemap_lead_snp_inconsistencies_c6d2e20f_upload.tsv` - the
  affected row for that upload's local copy (lead `_1`:
  `10:112974337_A_G` -> `10:112995025_C_T`)

Headline: 2,029 of 4,974 at-risk rows (41%) have an arbitrary lead SNP, across
1,275 studies / 1,341 study-blocks / 549 distinct LD blocks (out of ~731k
finemapped success rows; only ultra-strong signals were ever at risk). 66 rows
flip cis/trans (63 trans->cis, 3 cis->trans); mean lead-bp shift ~60 kb, max
1.34 Mb. Zero unreadable files. Rows with a single underflowed signal (2,945)
were already correct.

## Phase 3 - repair + downstream re-run (in this order)

### 3.0 Scan each machine that owns data

The main-pipeline tree scan here is done. On the upload server re-run the scan
so its `ld_blocks/gwas_upload/<guid>/` trees (incl. c6d2e20f) and its
pipeline-tree copy are covered:

```bash
cd <repo>/scripts   # inside the apptainer image, as run_pipeline.sh does
DATA_DIR=<DATA_DIR> apptainer exec --bind <DATA_DIR>:<DATA_DIR> \
  docker://mrcieu/genotype-phenotype-map:1.0.0 \
  Rscript find_finemap_lead_snp_inconsistencies.R \
  --diff_file finemap_lead_snp_inconsistencies.tsv \
  --affected_blocks_file finemap_lead_snp_affected_blocks.tsv \
  --mc_cores 16
```

Standalone upload copies (a directory containing `compiled_extracted_studies.tsv`)
can be scanned with `--study_dir <path>`.

### 3.1 Repair per-block metadata

Only rewrites `snp`, `bp`, `min_p`, `cis_trans` rows in per-block
`finemapped_studies.tsv`; the finemapped `_N.tsv.gz` / `_with_lbf.tsv.gz` files
are never modified, and susie is never re-run.

```bash
# dry run first, review the report, then:
DATA_DIR=<DATA_DIR> apptainer exec --bind <DATA_DIR>:<DATA_DIR> \
  docker://mrcieu/genotype-phenotype-map:1.0.0 \
  Rscript repair_finemap_lead_snps.R \
  --diff_file finemap_lead_snp_inconsistencies.tsv --apply TRUE
```

It skips rows whose metadata drifted since the scan (re-scan if so) and skips
`study_dir:` rows (standalone copies - repair the canonical tree instead).
Whole-file rewrites may reformat `min_p` floats at ulp level for untouched rows;
this is semantically irrelevant (all downstream comparisons are
`min_p <= threshold` with thresholds of 1e-6 / 1.5e-4).

Order matters: repair the main-pipeline tree BEFORE re-running upload
coloc/clustering (upload coloc partner selection at
coloc_studies_in_ld_block.R:405 and clustering's merge of main-pipeline coloc
pairs depend on main-pipeline lead bps). If both machines hold synced copies of
the main tree, repair on the owner and let the normal sync flow propagate it,
so a later sync does not overwrite the repair.

### 3.2 Main pipeline: invalidate coloc/clustering for affected blocks, then re-run

The coloc script only computes pairs missing from existing results (anti-join
at coloc_studies_in_ld_block.R:323), so stale per-block results must be deleted
to recompute pair sets from the repaired lead bps:

```bash
cd <repo root>
export DATA_DIR=<DATA_DIR>
tail -n +2 scripts/finemap_lead_snp_affected_blocks.tsv | cut -f2 | sort -u | while read -r block; do
  rm -f "$DATA_DIR/ld_blocks/$block/coloc_pairwise_results.tsv.gz"
  rm -f "$DATA_DIR/ld_blocks/$block"/clustered_results*.tsv.gz
  rm -f "$DATA_DIR/ld_blocks/$block"/igraph_clustered_results*.rds
  rm -f "$DATA_DIR/ld_blocks/$block/coloc_complete" "$DATA_DIR/ld_blocks/$block/clustering_complete"
done
./run_pipeline.sh
```

Snakemake then re-runs coloc -> clustering for exactly those blocks, then
`compile_results` -> `create_results_db` -> `create_static_web_files` (these
aggregate everything, so they re-run fully). The global-BFDR rule is commented
out of the Snakefile, so there is no cross-block BFDR ripple. `run_pipeline.sh`
runs `identify_studies_to_process.R` first, which the Snakefile requires.

### 3.3 Uploads (upload server), per affected guid

After repairing the guid's per-block trees in 3.1:

```bash
cd <repo>/pipeline_steps   # inside the apptainer image
export DATA_DIR=<DATA_DIR>
for block in <affected blocks of this guid from that machine's scan>; do
  rm -f "$DATA_DIR/ld_blocks/gwas_upload/<guid>/$block/coloc_pairwise_results.tsv.gz"
  rm -f "$DATA_DIR/ld_blocks/gwas_upload/<guid>/$block"/clustered_results*.tsv.gz
  rm -f "$DATA_DIR/ld_blocks/gwas_upload/<guid>/$block"/igraph_clustered_results*.rds
  rm -f "$DATA_DIR/ld_blocks/gwas_upload/<guid>/$block/coloc_complete" \
        "$DATA_DIR/ld_blocks/gwas_upload/<guid>/$block/clustering_complete"
  Rscript coloc_studies_in_ld_block.R --ld_block "$block" \
    --completed_output_file /tmp/coloc_complete_<guid> \
    --worker_guid <guid> \
    --worker_p_value_threshold <p>
    # add: --gwas_upload_ids_to_compare "guid1,guid2"   if the upload was processed with compare-guids
  Rscript cluster_studies_in_ld_block.R --ld_block "$block" \
    --completed_output_file /tmp/clustering_complete_<guid> \
    --worker_guid <guid>
done
```

`<p>` comes from the upload's `study_metadata.json` (`p_value_threshold`; for
c6d2e20f it is 5e-06), as does `compare_with_upload_guids`.

Then regenerate the upload's compiled outputs (mirrors
worker/pipeline_worker.R:323; this also re-uploads them to the Oracle bucket):

```bash
cd <repo>/pipeline_steps
Rscript -e '
source("constants.R")
update_directories_for_worker("<guid>")
source("../worker/compile_pipeline_results.R")
gwas_info <- list(metadata = jsonlite::fromJSON("<DATA_DIR>/gwas_upload/<guid>/study_metadata.json"))
compile_results(gwas_info)
'
```

`gwas_upload.db` is created empty by the main pipeline's `create_results_db`
step (3.2) and populated by the existing upload-server/web flow from these
compiled outputs - follow that flow after recompiling.

### 3.4 Verify

Re-run the scan (3.0) on each machine - expect 0 changed rows. Spot-check that
the c6d2e20f block's `_1` lead is now `10:112995025_C_T` in
`compiled_extracted_studies.tsv`, and that affected blocks'
`coloc_pairwise_results.tsv.gz` were rebuilt (fresh `bp_distance` values).

Known benign notes:

- 281 rows store `min_p = 0` where the true value is a subnormal (~1e-315) -
  identical under every downstream `min_p <= threshold` filter, no action
  needed.
- Unfinemapped rows (`finemap_message != "success"`) cannot be re-verified from
  disk (no P column in their files); the fix for their `which.min(P)` paths is
  forward-looking only.

## Phase 4 - tests

- Regression test "Finemapping lead SNP choice is not broken by p-value
  underflow" is in `tests/testthat/test_pipeline.R`. Its assertions are
  R-version-safe (verified on host R 4.6.1 and the image's R 4.2.2).
- Run: `make test` (~15 min; needs the test dataset at
  `/local-scratch/projects/genotype-phenotype-map/test/data/` and all packages
  from `docker/requirements.R`), or subsets via
  `Rscript tests/testthat/run_tests.R --only worker|pipeline`. On success,
  commit `tests/testing_complete.txt` (CI requires it with the branch name).
- `make lint` (lintr 3.3.0.1 via the container) is clean apart from the
  pre-existing `scripts/fix_se_sentinal.R` trailing-newline issue.
