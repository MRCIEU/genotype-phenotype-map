library(data.table)
library(jsonlite)

# Raw DICE_DB_1 eQTL VCFs (GRCh37, one file per immune cell type).
vcf_dir <- "/local-scratch/data/raw/dice"

# Intermediate esd/flist output (fed to make_besd.sh).
work_dir <- "/local-scratch/data/tmp_processes/dice"
esd_dir <- file.path(work_dir, "esd")
flist_dir <- file.path(work_dir, "flists")
liftover_work_dir <- file.path(work_dir, "liftover")

# GRCh38 gene coordinates (already lifted over) and the hg19 -> hg38 liftOver chain.
gene_map_file <- "/local-scratch/projects/genotype-phenotype-map/data/liftover/gene_name_map_full.tsv"
lift_over_binary <- "/local-scratch/projects/genotype-phenotype-map/data/liftover/liftOver"
lift_over_chain <- "/local-scratch/projects/genotype-phenotype-map/data/liftover/hg19ToHg38.over.chain.gz"

# DICE cell type labels and eQTL sample sizes (donors with expression + genotype,
# derived from the DICE TPM donor columns). New labels are added to cell_types in
# pipeline_steps/constants.R.
dice_metadata <- data.frame(
  file_prefix = c(
    "B_CELL_NAIVE", "MONOCYTES", "M2", "NK", "TREG_MEM", "CD4_NAIVE", "CD4_STIM",
    "TREG_NAIVE", "TFH", "TH1", "THSTAR", "TH17", "TH2", "CD8_NAIVE", "CD8_STIM"
  ),
  cell_type = c(
    "B IN", "Mono C", "Mono NC", "NK", "Treg Mem", "CD4 NC", "CD4 STIM",
    "Treg Naive", "TFH", "TH1", "THSTAR", "TH17", "TH2", "CD8 NC", "CD8 STIM"
  ),
  sample_size = c(91, 91, 90, 90, 89, 88, 88, 89, 89, 81, 88, 89, 89, 89, 88),
  stringsAsFactors = FALSE
)

awk_program <- c(
  "BEGIN { FS=\"\\t\"; OFS=\"\\t\" }",
  "/^#/ { next }",
  "{",
  "  gene=\"\"; gs=\"\"; pv=\"\"; beta=\"\"; stat=\"\";",
  "  n=split($8, a, \";\");",
  "  for (i=1; i<=n; i++) {",
  "    if (a[i] ~ /^Gene=/) { sub(/^Gene=/, \"\", a[i]); gene=a[i] }",
  "    else if (a[i] ~ /^GeneSymbol=/) { sub(/^GeneSymbol=/, \"\", a[i]); gs=a[i] }",
  "    else if (a[i] ~ /^Pvalue=/) { sub(/^Pvalue=/, \"\", a[i]); pv=a[i] }",
  "    else if (a[i] ~ /^Beta=/) { sub(/^Beta=/, \"\", a[i]); beta=a[i] }",
  "    else if (a[i] ~ /^Statistic=/) { sub(/^Statistic=/, \"\", a[i]); stat=a[i] }",
  "  }",
  "  print $1, $2, $4, $5, gene, gs, pv, beta, stat",
  "}"
)

parse_vcf <- function(vcf_file, awk_file) {
  awk_command <- paste0("zcat ", vcf_file, " | awk -f ", awk_file)
  return(data.table::fread(
    cmd = awk_command,
    header = FALSE,
    sep = "\t",
    col.names = c("chr", "bp", "ref", "alt", "gene", "gene_symbol", "p", "beta", "statistic"),
    showProgress = TRUE
  ))
}

prepare_variants <- function(dat) {
  dat <- dat[grepl("^chr", chr)]
  dat[, chr := sub("^chr", "", chr)]
  dat <- dat[chr %in% as.character(1:22)]
  dat[, gene := sub("\\..*", "", gene)]
  dat <- dat[gene != ""]
  dat[, `:=`(
    bp = as.integer(bp),
    beta = as.numeric(beta),
    statistic = as.numeric(statistic),
    p = as.numeric(p)
  )]
  dat <- dat[!is.na(bp) & !is.na(beta) & !is.na(statistic) & !is.na(p)]
  # DICE VCFs carry the ANOVA statistic (beta / se) rather than se.
  dat[, se := abs(beta / statistic)]
  return(dat)
}

liftover_variants <- function(dat, prefix) {
  unique_variants <- unique(dat[, .(chr, bp)])
  bed <- data.table(
    chrom = paste0("chr", unique_variants$chr),
    start = unique_variants$bp - 1,
    end = unique_variants$bp,
    name = paste0(unique_variants$chr, ":", unique_variants$bp)
  )
  bed_in <- file.path(liftover_work_dir, paste0(prefix, ".bed"))
  bed_out <- file.path(liftover_work_dir, paste0(prefix, ".hg38.bed"))
  unmapped <- file.path(liftover_work_dir, paste0(prefix, ".unmapped.bed"))

  data.table::fwrite(bed, bed_in, sep = "\t", col.names = FALSE)
  system2(lift_over_binary, args = c(bed_in, lift_over_chain, bed_out, unmapped))

  lifted <- data.table::fread(
    bed_out,
    header = FALSE,
    sep = "\t",
    col.names = c("chrom", "start", "end", "name"),
    showProgress = FALSE
  )
  lifted <- lifted[!duplicated(name)]
  lifted[, `:=`(
    bp38 = start + 1,
    chr37 = sub(":.*", "", name),
    bp37 = as.numeric(sub(".*:", "", name))
  )]

  dat <- merge(
    dat,
    lifted[, .(chr = chr37, bp = bp37, bp38)],
    by = c("chr", "bp"),
    all.x = FALSE
  )
  dat[, bp := as.integer(bp38)]
  dat[, bp38 := NULL]
  return(dat)
}

standardise_dice_alleles <- function(dat) {
  dat[, A1 := toupper(alt)]
  dat[, A2 := toupper(ref)]
  to_flip <- dat$A1 > dat$A2
  if (any(to_flip)) {
    temp <- dat$A2[to_flip]
    dat$A2[to_flip] <- dat$A1[to_flip]
    dat$A1[to_flip] <- temp
    dat$beta[to_flip] <- -1 * dat$beta[to_flip]
  }
  return(dat)
}

join_gene_annotations <- function(dat, gene_map) {
  dat <- merge(dat, gene_map, by.x = "gene", by.y = "ENSEMBL_ID", all.x = TRUE)
  dat[, ProbeBp := BP_START]
  # Fall back to the most significant variant for genes absent from the map.
  dat[is.na(ProbeBp), ProbeBp := bp[which.min(p)], by = gene]
  dat[, Gene := ifelse(is.na(GENE_NAME) | GENE_NAME == "", gene_symbol, GENE_NAME)]
  dat[, Orientation := "+"]
  dat[, `:=`(BP_START = NULL, GENE_NAME = NULL)]
  return(dat)
}

write_esd_files <- function(dat, esd_out_dir) {
  dir.create(esd_out_dir, recursive = TRUE, showWarnings = FALSE)
  dat[, {
    data.table::fwrite(
      .SD[, .(Chr = chr, SNP, Bp = bp, A1, A2, Freq = NA_real_, Beta = beta, se, p)],
      file.path(esd_out_dir, paste0(gene, ".esd")),
      sep = "\t",
      na = "NA",
      quote = FALSE
    )
    NULL
  }, by = gene]
  return(invisible(NULL))
}

write_flist <- function(dat, flist_file, esd_out_dir) {
  probes <- unique(dat[, .(Chr = chr, ProbeID = gene, ProbeBp = ProbeBp, Gene = Gene, Orientation = Orientation)])
  probes[, Chr := as.integer(Chr)]
  probes[, GeneticDistance := 0]
  probes[, PathOfEsd := file.path(esd_out_dir, paste0(ProbeID, ".esd"))]
  data.table::setcolorder(
    probes,
    c("Chr", "ProbeID", "GeneticDistance", "ProbeBp", "Gene", "Orientation", "PathOfEsd")
  )
  data.table::fwrite(probes, flist_file, sep = "\t")
  return(invisible(NULL))
}

write_metadata_json <- function(metadata, json_file) {
  jsonlite::write_json(
    list(
      cell_type = metadata$cell_type,
      tissue = "Whole Blood",
      sample_size = metadata$sample_size,
      ancestry = "EUR",
      cis_trans = "cis",
      category = "continuous"
    ),
    json_file,
    pretty = TRUE,
    auto_unbox = TRUE
  )
  return(invisible(NULL))
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  vcf_files <- list.files(vcf_dir, pattern = "\\.vcf\\.gz$", full.names = TRUE)
  if (length(args) > 0) {
    prefixes <- sub("\\.vcf\\.gz$", "", basename(vcf_files))
    vcf_files <- vcf_files[prefixes %in% args]
  }
  if (length(vcf_files) == 0) {
    stop("No VCF files found to process")
  }

  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(esd_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(flist_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(liftover_work_dir, recursive = TRUE, showWarnings = FALSE)

  awk_file <- file.path(work_dir, "parse_dice_vcf.awk")
  writeLines(awk_program, awk_file)

  gene_map <- data.table::fread(
    gene_map_file,
    select = c("ENSEMBL_ID", "GENE_NAME", "BP_START"),
    showProgress = FALSE
  )
  gene_map <- unique(gene_map, by = "ENSEMBL_ID")

  for (vcf_file in vcf_files) {
    prefix <- sub("\\.vcf\\.gz$", "", basename(vcf_file))
    metadata <- dice_metadata[dice_metadata$file_prefix == prefix, ]
    if (nrow(metadata) != 1) {
      stop(paste0("No cell type metadata defined for VCF prefix: ", prefix))
    }

    message(paste0("Processing ", prefix))
    dat <- parse_vcf(vcf_file, awk_file)
    dat <- prepare_variants(dat)
    dat <- liftover_variants(dat, prefix)
    dat <- standardise_dice_alleles(dat)
    dat <- join_gene_annotations(dat, gene_map)
    dat[, SNP := paste0(chr, ":", bp, "_", A1, "_", A2)]

    esd_out_dir <- file.path(esd_dir, prefix)
    write_esd_files(dat, esd_out_dir)
    write_flist(dat, file.path(flist_dir, paste0(prefix, ".flist")), esd_out_dir)
    write_metadata_json(metadata, file.path(flist_dir, paste0(prefix, ".json")))

    rm(dat)
    gc()
  }

  message("Finished formatting DICE eQTLs. Run make_besd.sh to build BESD files.")
  return(invisible(NULL))
}

main()
