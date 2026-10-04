# AGENTS.md

You are an expert, highly autonomous developer. Your goal is to write clean, maintainable, and pragmatic code, while strictly adhering to the user's architectural guidelines.

# Output and Communication Style

- **No Preamble:** Start answering immediately. Do not say "Sure, I can help with that."
- **Keep explanations brief:** Answer the question first, then offer to go deeper. Bullet points are preferred over long paragraphs.
- **Never summarize user prompts:** Do not restate or paraphrase what the user asked you to do.
- **Stop writing when finished:** Do not pad responses or give unsolicited advice unless it relates to a critical bug/security issue.

## Repo overview

GWAS/QTL summary-statistics pipeline: ingest studies (opengwas/besd/tsv/summary_stats) → extract regions → per-LD-block standardisation → imputation → finemapping → coloc → clustering → compile → SQLite DBs → static web files. Snakemake (`Snakefile`) orchestrates R scripts in `pipeline_steps/`. A separate Redis-based worker (`worker/pipeline_worker.R`) processes GWAS uploads and shares the same extraction/standardisation/imputation/finemap code.

## Environment

- All R/Python run inside apptainer image `docker://mrcieu/genotype-phenotype-map:1.0.0` (built from `docker/Dockerfile`). `run_pipeline.sh` builds `APPTAINER_VARS` with required bind mounts.
- `run_pipeline.sh` loads `.env` (use `.env.pipeline_local` / `.env.pipeline_worker` as templates; `.env` is gitignored) and accepts `--version/-v` and `--block-list` args. It must run `identify_studies_to_process.R` **before** `snakemake` — the `Snakefile` reads `$DATA_DIR/pipeline_metadata/studies_to_process.tsv` at import time and exits if it is missing/empty. Snakemake's profile is the repo root (`config.yaml`).
- The `Snakefile` does `os.chdir('pipeline_steps')`; R scripts resolve paths via `constants.R` from `DATA_DIR`/`RESULTS_DIR` env vars.
- Setting `TEST_RUN` makes code read LD blocks from `tests/data/ld_blocks.tsv` instead of `pipeline_steps/data/ld_blocks.tsv`.

## Commands

- `make lint` — lintr on all R files. lintr is pinned to **3.3.0.1** (in `.github/workflows/main.yml` and `docker/requirements.R`); a different version gives different results.
- `make format` — styler; `make lint-summary` — lint summary.
- `make test` — `Rscript tests/testthat/run_tests.R`, ~15 min.

## Testing (critical)

- Tests are **not** run in CI. Run locally before merging. They require the large test datasets at the hard-coded paths `/local-scratch/projects/genotype-phenotype-map/test/e2e/` (pipeline/worker) and `.../test/unit/` (step tests; fixtures built once by `tests/generate_unit_fixtures.R`) (see `tests/testthat/run_tests.R`) plus all R packages from `docker/requirements.R`.
- On success, tests overwrite `tests/testing_complete.txt` with `SUCCESS: All tests passed on branch: <branch>`. CI fails if this file is missing or lacks the PR branch name — so after running tests, **commit `tests/testing_complete.txt`**.
- Subset runs: `Rscript tests/testthat/run_tests.R --only worker|pipeline|unit`; add `--dont_delete` to keep test output for debugging.

## Architecture notes

- Per-LD-block work lives under `$DATA_DIR/ld_blocks/{ancestry}/{chr}/{start}-{stop}/`, gated by `*_complete` sentinel files (e.g. `standardisation_complete`, `imputation_complete`, `finemapping_complete`, `coloc_complete`, `clustering_complete`, `compare_rare_complete`) used as Snakemake outputs. Rules skip a block by touching the sentinel when no studies are in it.
- `pipeline_steps/constants.R` is the shared config/source of truth for every R script; `database_definitions.R` holds DB schemas; `gwas_calculations.R` holds shared stats code. `SVG_*`/plot helpers are in `svg_helpers.R`.
- Results are written under `$RESULTS_DIR/<PIPELINE_VERSION>/` (`current` by default); `latest/` mirrors it.

## Implementation Guide

- **Apply requested behavior only:** Do not refactor, extract, or reorganize files unless explicitly asked. If a larger change is suggested, ask first before implementing.
- **Minimal Changes Unless Otherwise Requested:** I want you to always make more minimal changes unless I ask you more explicitly.
- **Keep it simple:** Before adding a private helper method, count the call sites. If there is only one, inline it unless the user asks otherwise.
- You don't need to always verify your changes by writing bespoke scripts and running them.  Sometimes they take a very long time to run.  Ask before doing any of that, I do have tests for this reason.


## Code Style

- **Follow Convension:** Always follow the convension, or library choice of what already exists in the file or repository.  Only use ASCII characters in code
- **Simplicity:** Favor readable, explicit code over clever, condensed one-liners. 
- Default to using dplyr, and don't use `.data$` syntax if you don't have to
- Follow `.lintr`: 2-space indent, 120-col line limit, `<-` (and `<<-`) assignment, explicit `return()`. Object-name linting is disabled.
