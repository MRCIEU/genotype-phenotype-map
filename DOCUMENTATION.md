# Documentation

## Repository structure

| Directory | Contents |
|-----------|----------|
| `pipeline_steps/` | Core pipeline R scripts (extraction, standardisation, imputation, finemapping, coloc, compile); `impute_studies_in_ld_block.py`; `data/` with `study_list.csv`, `study_sources.csv`, `ld_blocks.tsv`, etc.; shared helpers (`constants.R`, `database_definitions.R`, `gwas_calculations.R`, `svg_helpers.R`) |
| `worker/` | `pipeline_worker.R` (Redis-based job processor), `redis_client.R` |
| `scripts/` | One-off utilities: `sync_to_servers.R`, `remove_studies_from_pipeline.R`, `update_studies_in_pipeline.R`, `create_test_db_from_results.R`; `data_extraction/` (BESD formatting, rare-variant pulls); `data_formatting/` (traits, VEP); `analysis/` (RMarkdown reports); `prepare_ld_data/` |
| `docker/` | `Dockerfile` (main pipeline image), `uploadDockerfile`, `requirements.R`, `requirements.txt` |
| `tests/` | `testthat/` tests, `data/` fixtures, `lint.R`, `lint_format.R`, `lint_summary.R` |
| `docs/` | Method docs and illustrations (`basic_simulation/`, `cistrans_pqtl/`, `region_illustration/`) |
| `.github/` | CI workflows |

**Root:** `Snakefile`, `config.yaml`, `run_pipeline.sh`, `Makefile`

---

## Adding new data to the pipeline

Studies are configured in `pipeline_steps/data/study_list.csv`. Each row defines a set of studies to ingest. Required columns:

| Column | Description |
|--------|-------------|
| `id_pattern` | Glob pattern to match study IDs (e.g. `ebi-a-*`, `*` for all in directory) |
| `data_type` | One of: `phenotype`, `gene_expression`, `splice_variant`, `transcript`, `protein`, `methylation`, `metabolome`, `cell_trait`, `plasma_protein` |
| `source` | Must exist in `pipeline_steps/data/study_sources.csv` |
| `data_format` | `opengwas`, `besd`, `tsv`, or `summary_stats` |
| `bespoke_parsing` | Parsing variant (e.g. `none`, `gtex_sqtl`, `godmc`) |
| `data_location` | Path to the data on disk |
| `ancestry` | e.g. `EUR` |
| `reference_build` | `GRCh37` or `GRCh38` |
| `p_value_threshold` | e.g. `1.5e-4` |
| `variant_type` | `common` or `rare_exome` |
| `metadata_type` | `json`, `tsv`, or `csv` (`csv` required for `summary_stats`) |
| `coverage` | `dense` or `sparse` |

The pipeline discovers studies by matching `data_location` + `id_pattern` against the filesystem. Study IDs are derived from folder and file names.

### summary_stats format

For delimited GWAS uploads (csv/tsv), use `data_format=summary_stats` and `metadata_type=csv`. Point `data_location` at a directory containing `metadata.csv`, where each row describes one summary-stat file (same fields as the worker upload JSON):

| Column | Description |
|--------|-------------|
| `study_name` | Unique study ID (underscores become dashes) |
| `file_location` | Path to the summary-stat file (absolute, or relative to `data_location`) |
| `trait` | Trait display name |
| `sample_size` | Sample size |
| `category` | `continuous` or `categorical` |
| `file_type` | `csv` (delimited; VCF not yet supported) |
| `source` | Unique source slug for this file/paper (used in the DB `study_sources` table) |
| `name` | Display name for the source (e.g. paper or consortium name) |
| `url` | Source URL (may be empty) |
| `doi` | DOI URL or identifier (may be empty) |
| `SNP`, `RSID`, `CHR`, `BP`, `EA`, `OA`, `P`, `BETA`, `OR`, `SE`, `EAF`, `N` | Optional column map: value is the source column name in the file |

Ancestry, reference build, and p-value threshold come from the `study_list.csv` row. Extraction reuses `extract_regions_from_summary_stats.R` (same script as the worker path).

For `study_list.csv`, use `source=manual` (or another catch-all already in `study_sources.csv`). Each study’s `source` slug is taken from `metadata.csv` and stored on the study as usual. Citation fields (`name` / `url` / `doi`) stay in `metadata.csv` only and are merged into the DB `study_sources` table when databases are created — they are not written to `studies_processed`. Use a unique `source` slug per paper/dataset.

To exclude studies, add them to `pipeline_steps/data/ignore_studies.tsv` or `pipeline_steps/data/ignore_studies_rare.tsv`.

### study_sources.csv

`pipeline_steps/data/study_sources.csv` defines the data sources referenced in `study_list.csv`. Each row has:

| Column | Description |
|--------|-------------|
| `source` | Unique identifier (must match `source` in study_list.csv) |
| `name` | Display name |
| `url` | Source URL |
| `doi` | DOI for citation |

When adding a new source to `study_list.csv`, add a corresponding row to `study_sources.csv`. For one-off summary-stat ingestions, prefer putting citation fields in `metadata.csv` instead (see above).

---

## Ancillary data needed to run pipeline

These must exist under `$DATA_DIR` before running:

1. **LD reference panel** (`$DATA_DIR/ld_reference_panel_hg38/<ANCESTRY>/<CHR>/<START>-<STOP>.*`)
   - Precomputed LD matrices per LD block
   - Format: `.unphased.vcor1.gz`, `.ldeig.rds`, `.tsv` (info file)
   - LD blocks are defined in `pipeline_steps/data/ld_blocks.tsv` (from [ldetect](https://github.com/jmacdon/LDblocks_GRCh38))

2. **Variant annotation** (`$DATA_DIR/variant_annotation/`)
   - `vep_annotations_hg38.tsv.gz` – VEP annotations for SNP consequences
   - `gene_info.tsv` – Gene metadata

3. **Liftover data** (`$DATA_DIR/liftover/`)
   - Chain files for GRCh37 → GRCh38 (for studies on hg19/GRCh37)

4. **Study data**
   - Paths specified in `study_list.csv`; typically under `/local-scratch/data/`
   - Must match the expected layout for the chosen `data_format` (opengwas, besd, tsv, summary_stats)

---

## Data organisation

### Data directory (`$DATA_DIR`)

Populated by the pipeline:

| Path | Contents |
|------|----------|
| `study/<study_id>/` | Per-study outputs: `extracted_snps.tsv`, and subdirs `extracted/`, `standardised/`, `imputed/`, `finemapped/` |
| `ld_blocks/<ancestry>/<chr>/<start>-<stop>/` | LD-block outputs: extracted regions, standardisation, imputation, finemapping, coloc results |
| `pipeline_metadata/` | `studies_to_process.tsv`, `updated_ld_blocks_to_colocalise.tsv`, `logs/snakemake.log` |
| `logs/` | Pipeline worker logs (e.g. `pipeline_worker_YYYYMM.log`) |
| `rsync_to_server/` | Staging directory for `scripts/sync_to_servers.R` before rsync to Oracle |
| `study/*/svgs/` | Per-study SVG visualisations (extraction plots) |
| `svgs/` | svg files that have been re-organised after db creation |

Pre-existing (you must provide):

| Path | Contents |
|------|----------|
| `ld_reference_panel_hg38/<ancestry>/<chr>/` | LD matrices per block |
| `variant_annotation/` | VEP annotations, gene info |
| `liftover/` | Chain files for coordinate liftover |

### Results directory (`$RESULTS_DIR`)

| Path | Contents |
|------|----------|
| `latest/` | Latest run outputs: `studies_processed.tsv.gz`, `traits_processed.tsv.gz`, `study_extractions.tsv.gz`, `analysis/` |
| `current/` | Current pipeline run outputs (also symlinked/copied for “latest”) |
| `<timestamp>/` or `<version>/` | Archived results for a specific run (e.g. `2025_02_24-12_30/` or `1.0.0/`), copied from `current/` on success |

### Results files and creating scripts

| File(s) | Script |
|---------|--------|
| `studies_processed.tsv.gz`, `traits_processed.tsv.gz`, `study_extractions.tsv.gz` | `compile_results.R` |
| `coloc_pairwise_results.tsv.gz` | `compile_results.R` (aggregates per-LD-block results from `coloc_studies_in_ld_block.R`, then annotated by `cluster_studies_in_ld_block.R`) |
| `coloc_clustered_results.tsv.gz` | `compile_results.R` (aggregates per-LD-block clustering from `cluster_studies_in_ld_block.R`, using the global threshold from `calculate_global_bfdr.R`) |
| `pipeline_metadata/global_bfdr.tsv` | `calculate_global_bfdr.R` (synced to worker servers for GWAS upload clustering) |
| `rare_results.tsv.gz` | `compile_results.R` (aggregates from `compare_rare_studies_in_ld_block.R` per LD block) |
| `studies.db`, `associations.db`, `associations_full.db`, `coloc_pairs.db`, `coloc_pairs_full.db`, `ld.db`, `gwas_upload.db` | `create_db_from_results.R` |
| `static_web/opengwas_ids.json`, `phenotype_id_map.json`, `robots.txt`, `sitemap.xml`, `pipeline_summary.html` | `create_static_web_files.R` |
| `static_web/svgs/` | `create_static_web_files.R` (via `prepare_svg_files_for_use`) |

---

## Concepts

- **Study:** A single summary-statistics dataset (e.g. a GWAS like BMI, or a gene–tissue eQTL like GTEx Whole Blood WASH7P).
- **LD block:** Genome regions with minimal LD overlap, from [ldetect](https://github.com/jmacdon/LDblocks_GRCh38).
- **Docker / Apptainer:** All CLI tools, R and Python packages run in [mrcieu/genotype-phenotype-map](https://hub.docker.com/repository/docker/mrcieu/genotype-phenotype-map).
- **Pipeline:** Snakemake orchestrates extraction, standardisation, imputation, finemapping, and colocalisation. See `Snakefile` and `config.yaml`.
