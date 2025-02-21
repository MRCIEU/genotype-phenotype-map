source('constants.R')

library(dplyr)
library(duckdb)
library(data.table)
library(validate)
library(tidyr)
library(R.utils)
library(purrr)
library(furrr)
library(parallel)

parser <- argparser::arg_parser('Create DuckDB from pipeline results')
#INPUT
parser <- argparser::add_argument(parser, '--results_dir', help = 'Results directory of pipeline', type = 'character')
#OUTPUT
parser <- argparser::add_argument(parser, '--studies_db_file', help = 'DuckDB studies file', type = 'character')
parser <- argparser::add_argument(parser, '--associations_db_file', help = 'DuckDB associations file', type = 'character')

args <- argparser::parse_args(parser)

main <- function() {
  all_study_blocks <- fread(file.path(args$results_dir, "all_study_blocks.tsv"))
  str(all_study_blocks)
  raw_coloc_results <- fread(file.path(args$results_dir, "raw_coloc_results.tsv"))
  str(raw_coloc_results)
  # rare_results <- fread(file.path(results_dir, "rare_results.tsv"))
  # str(rare_results)
  results_metadata <- fread(file.path(args$results_dir, "results_metadata.tsv"))
  str(results_metadata)
  studies_processed <- fread(file.path(args$results_dir, "studies_processed.tsv"))
  str(studies_processed)
  variant_annotations <- fread(file.path(args$results_dir, "variant_annotations.tsv"))
  str(variant_annotations)
  variant_annotations_full <- fread(file.path(variant_annotation_dir, "vep_annotations_hg38.tsv.gz"))
  names(variant_annotations)
  names(variant_annotations_full)

  results_metadata <- results_metadata %>%
    tidyr::separate(ld_block, into=c("pop", "chr", "start", "end"), sep="[/-]", remove=FALSE) %>%
    mutate(chr=as.integer(chr), start=as.integer(start), end=as.integer(end))


  vl <- mclapply(1:nrow(results_metadata), \(i) {
    message(i)
    chr <- results_metadata$chr[i]
    start <- results_metadata$start[i]
    end <- results_metadata$end[i]
    ind <- variant_annotations_full$CHR == chr & variant_annotations_full$BP >= start & variant_annotations_full$BP < end
    tibble(SNP=variant_annotations$SNP[ind], ld_block=results_metadata$ld_block[i])
  }, mc.cores=50) %>% bind_rows()

  variant_annotations_full <- left_join(variant_annotations_full, vl, by="SNP")
  str(variant_annotations_full)


  coloc <- raw_coloc_results %>%
    filter(posterior_prob > 0.5) %>%
    mutate(id=1:n()) %>%
    separate_longer_delim(cols=traits, delim=", ")

  s <- all_study_blocks %>%
    select(study, unique_study_id, chr, bp, min_p, cis_trans, ld_block, known_gene)

  coloc <- coloc %>%
    left_join(s, by=c("traits"="unique_study_id"))

  str(coloc)


  # all_study_blocks includes variants that don't colocalise with anything. Need to get those, but they don't have variant IDs
  head(all_study_blocks)
  dim(all_study_blocks)
  temp <- select(variant_annotations_full, SNP, CHR, BP)
  all_study_blocks <- left_join(all_study_blocks, temp, by=c("chr"="CHR", "bp"="BP")) %>% filter(!duplicated(unique_study_id))
  dim(all_study_blocks)

  # Find finemapped results that colocalise with nothing
  no_coloc <- subset(all_study_blocks, !unique_study_id %in% coloc$traits) %>% rename(candidate_snp=SNP, traits=unique_study_id) %>% mutate(id=(1:nrow(.))+nrow(coloc))
  coloc <- bind_rows(coloc, no_coloc)

  # Add ld_block to variant_annotations
  temp <- coloc %>% filter(!duplicated(candidate_snp)) %>% select(candidate_snp, ld_block)
  variant_annotations <- left_join(variant_annotations, temp, by=c("SNP"="candidate_snp"))
  dim(variant_annotations)
  head(variant_annotations)


  # Generate LD
  ld_dir <- file.path(data_dir, "/ld_reference_panel_hg38")

  generate_ld_obj <- function(ld_dir, ld_block, all_study_blocks) {
    message(ld_block)
    ld <- suppressMessages(fread(file.path(ld_dir, paste0(ld_block, ".unphased.vcor1"))))
    # ld <- readr::read_tsv(file.path(ld_dir, paste0(ld_block, ".unphased.vcor1")), col_names=FALSE)
    ldvars <- suppressMessages(scan(file.path(ld_dir, paste0(ld_block, ".unphased.vcor1.vars")), character()))
    names(ld) <- ldvars
    ld$lead <- ldvars
    ind <- which(ldvars %in% all_study_blocks$SNP)
    ld <- ld[ind,]
    ldl <- tidyr::pivot_longer(ld, cols=-lead, names_to="variant", values_to="r") %>% filter(r^2 > 0.8 | variant %in% ld$lead) %>% filter(lead != variant) %>% mutate(ld_block=ld_block)
    return(ldl)
  }

  ld_blocks <- unique(all_study_blocks$ld_block)
  generate_ld_obj(ld_dir, ld_blocks[1], variant_annotations)
  ldl <- mclapply(ld_blocks, \(x) generate_ld_obj(ld_dir, x, variant_annotations), mc.cores=30) %>% bind_rows()


  processed_con <- dbConnect(duckdb::duckdb(), args$studies_db_file)
  dbWriteTable(processed_con, "all_study_blocks", all_study_blocks)
  dbWriteTable(processed_con, "results_metadata", results_metadata)
  dbWriteTable(processed_con, "studies_processed", studies_processed)
  dbWriteTable(processed_con, "variant_annotations", variant_annotations_full)
  dbWriteTable(processed_con, "coloc", coloc)
  dbWriteTable(processed_con, "ld", ldl)

  dbDisconnect(processed_con, shutdown=TRUE)

  # Extract summary statistics
  varids <- unique(all_study_blocks$SNP)
  length(varids)
  varids_info <- tibble(varids) %>% separate(varids, into=c("chr", "other"), sep=":", remove=FALSE)
  varids_list <- lapply(1:22, \(x) {
    subset(varids_info, chr==x)$varids
  })

  sources <- unique(studies_processed$source)

  assocs <- lapply(sources, \(x) assocs_source(studies_processed, varids_list, x, 30))
  unlink(args$associations_db_file)
  assocs_con <- dbConnect(duckdb::duckdb(), args$associations_db_file)
  dbWriteTable(assocs_con, "assocs", assocs[[1]])
  for(i in 2:length(assocs)) {
    if(nrow(assocs[[i]]) > 0) {
      dbAppendTable(assocs_con, "assocs", assocs[[i]])
    }
  }
  dbDisconnect(assocs_con, shutdown=TRUE)

  ensure_dbs_are_valid()
}


extract_variants <- function(varids_list, path, study="study") {
  file_list <- list.files(path) %>% file.path(path, .) %>% grep("pre_filter", ., value=TRUE, invert=TRUE) %>% grep("dentist", ., value=TRUE, invert=TRUE)
  if(length(file_list) == 0) {
    return(NULL)
  }
  ext <- lapply(file_list, \(y) {
    chr <- strsplit(basename(y), "_")[[1]][2] %>% as.numeric()
    tryCatch({
      fread(y) %>% filter(SNP %in% varids_list[[chr]]) %>% select(SNP, BETA, SE, IMPUTED, P, EAF) %>% mutate(study=study)
    }, error=function(e) {
      message(e)
      return(NULL)
    })
  }) 
  ext <- ext[!sapply(ext, is.null)] %>% bind_rows()
  return(ext)
}

assocs_source <- function(studies_processed, varids_list, source, mc.cores=10) {
  assocs <- mclapply(studies_processed$study_name[studies_processed$source == source], \(x) {
    message(x)
    path <- file.path(data_dir, "/study/", x, "imputed")
    tryCatch({
      extract_variants(varids_list, path, x)
    }, error=function(e) {
      message(e)
      return(NULL)
    })
  }, mc.cores=mc.cores)
  assocs <- assocs[!sapply(assocs, \(x) inherits(x, "try-error"))]
  return(assocs %>% bind_rows())
}

ensure_dbs_are_valid <- function() {
  #TODO: Check that the databases are valid
}

main()