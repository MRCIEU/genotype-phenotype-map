library(data.table)
library(jsonlite)

# Locate the repo so we can reuse the shared extraction helpers. Override with
source(file.path("pipeline_steps/common_extraction_functions.R"))

# Raw DICE_DB_1 eQTL VCFs (GRCh37, one file per immune cell type).
vcf_dir <- "/local-scratch/data/raw/dice"

# Intermediate esd/flist output (fed to make_besd.sh).
work_dir <- "/local-scratch/data/tmp_processes/dice"
esd_dir <- file.path(work_dir, "esd")
flist_dir <- file.path(work_dir, "flists")
liftover_work_dir <- file.path(work_dir, "liftover")

# GRCh38 gene coordinates and the hg19 -> hg38 liftOver resources used by the
# pipeline (see pipeline_steps/constants.R).
liftover_dir <- "/local-scratch/projects/genotype-phenotype-map/data/liftover"
gene_map_file <- file.path(liftover_dir, "gene_name_map_full.tsv")
bcftools <- "/home/bcftools/bcftools"
hg19_fasta <- file.path(liftover_dir, "hg19.fa")
hg38_fasta <- file.path(liftover_dir, "hg38.fa")
lift_over_chain <- file.path(liftover_dir, "hg19ToHg38.over.chain.gz")

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

# Raw DICE lung eQTL TSVs (GRCh38, one file per immune cell subset). The files are
# named after the download page; non-ASCII characters are transliterated before
# matching so this source stays ASCII. Sample sizes are the MatrixEQTL t-statistic
# degrees of freedom plus the 7 fitted parameters (intercept, genotype, 1 genotype
# PC, 2 PEER factors, sex, age), matching the 120 study participants.
dice_lung_metadata <- data.frame(
  file_prefix = c(
    "DICE_LUNG_eQTL_Atypical B cells", "DICE_LUNG_eQTL_CD14+ monocytes",
    "DICE_LUNG_eQTL_CD16+ monocytes", "DICE_LUNG_eQTL_CD4+ CTL (1)",
    "DICE_LUNG_eQTL_CD4+ CTL (2)", "DICE_LUNG_eQTL_CD4+ TCM cells",
    "DICE_LUNG_eQTL_CD4+ TFH cells", "DICE_LUNG_eQTL_CD4+ THIFNR cells",
    "DICE_LUNG_eQTL_CD4+ TREG cells", "DICE_LUNG_eQTL_CD4+ TRM ITGAE high",
    "DICE_LUNG_eQTL_CD4+ TRM ITGAE low", "DICE_LUNG_eQTL_CD8+ GZMK+ cells",
    "DICE_LUNG_eQTL_CD8+ IFNR cells", "DICE_LUNG_eQTL_CD8+ KIR+ cells",
    "DICE_LUNG_eQTL_CD8+ MAIT cells", "DICE_LUNG_eQTL_CD8+ TCM cells",
    "DICE_LUNG_eQTL_CD8+ TEM cells", "DICE_LUNG_eQTL_CD8+ TRM cells",
    "DICE_LUNG_eQTL_Conventional DC1", "DICE_LUNG_eQTL_Conventional DC2",
    "DICE_LUNG_eQTL_DC3", "DICE_LUNG_eQTL_Macrophages",
    "DICE_LUNG_eQTL_Mature regulatory DC", "DICE_LUNG_eQTL_Memory B cells",
    "DICE_LUNG_eQTL_Naive B cells", "DICE_LUNG_eQTL_NK CD16+",
    "DICE_LUNG_eQTL_NK CD16-", "DICE_LUNG_eQTL_Plasma cells",
    "DICE_LUNG_eQTL_Plasmacytoid DC"
  ),
  cell_type = c(
    "B Atypical", "Mono CD14", "Mono CD16", "CD4 CTL 1", "CD4 CTL 2", "CD4 TCM",
    "CD4 TFH", "CD4 THIFNR", "CD4 Treg", "CD4 TRM ITGAE hi", "CD4 TRM ITGAE lo",
    "CD8 GZMK", "CD8 IFNR", "CD8 KIR", "CD8 MAIT", "CD8 TCM", "CD8 TEM", "CD8 TRM",
    "cDC1", "cDC2", "DC3", "Macrophages", "mregDC", "B Mem", "B Naive",
    "NK CD16 pos", "NK CD16 neg", "Plasma", "pDC"
  ),
  sample_size = c(
    58, 120, 120, 108, 111, 120, 97, 75, 119, 120, 120, 120, 76, 111, 112, 120, 120,
    120, 82, 116, 119, 118, 65, 117, 111, 120, 115, 49, 95
  ),
  stringsAsFactors = FALSE
)

query_format <- paste0(
  paste(
    c(
      "%CHROM", "%POS", "%REF", "%ALT", "%INFO/Gene", "%INFO/GeneSymbol",
      "%INFO/Pvalue", "%INFO/Beta", "%INFO/Statistic", "%INFO/SWAP"
    ),
    collapse = "\\t"
  ),
  "\\n"
)

write_contig_header <- function(fai_file, output_file) {
  fai <- data.table::fread(
    fai_file,
    header = FALSE,
    sep = "\t",
    select = c(1, 2),
    col.names = c("contig", "length"),
    showProgress = FALSE
  )
  writeLines(paste0("##contig=<ID=", fai$contig, ",length=", fai$length, ">"), output_file)
  return(invisible(NULL))
}

# Lift the raw VCF from GRCh37 to GRCh38 with bcftools, then stream the fields we
# need back into R. bcftools +liftover may swap REF/ALT to match the hg38
# reference, which is recorded in INFO/SWAP and handled below.
parse_lifted_vcf <- function(vcf_file, prefix) {
  rejected_file <- file.path(liftover_work_dir, paste0(prefix, ".rejected.vcf"))
  bcftools_command <- paste0(
    bcftools, " annotate --header-lines ", file.path(liftover_work_dir, "contigs.txt"),
    " -Ou ", vcf_file, " | ",
    bcftools, " +liftover --no-version -Ou -- -s ", hg19_fasta, " -f ", hg38_fasta,
    " -c ", lift_over_chain, " --reject ", rejected_file, " | ",
    bcftools, " query -f '", query_format, "'"
  )
  dat <- data.table::fread(
    cmd = bcftools_command,
    header = FALSE,
    sep = "\t",
    col.names = c("chr", "bp", "ref", "alt", "gene", "gene_symbol", "p", "beta", "statistic", "swap"),
    showProgress = TRUE
  )
  if (nrow(dat) == 0) {
    stop(paste0("bcftools produced no variants for ", vcf_file))
  }
  return(dat)
}

prepare_variants <- function(dat) {
  # Multi-allelic rows (SWAP = -1) get a new reference allele from +liftover and
  # do not map onto a single effect allele, so drop them.
  dat <- dat[!grepl(",", alt)]
  dat[, chr := sub("^chr", "", chr)]
  dat <- dat[chr %in% as.character(1:22)]
  dat[, gene := sub("\\..*", "", gene)]
  dat <- dat[gene != "" & gene != "."]
  dat[, `:=`(
    chr = as.integer(chr),
    bp = as.integer(bp),
    beta = as.numeric(beta),
    statistic = as.numeric(statistic),
    p = as.numeric(p)
  )]
  dat <- dat[!is.na(bp) & !is.na(beta) & !is.na(statistic) & !is.na(p)]
  # DICE VCFs carry the ANOVA statistic (beta / se), not se.
  dat[, se := abs(beta / statistic)]
  # When +liftover swaps REF/ALT (SWAP = 1), the ALT-relative beta changes sign.
  swapped <- as.character(dat$swap) == "1"
  dat[swapped, beta := -1 * beta]
  return(dat)
}

# The lung TSVs are already GRCh38 and carry beta/se directly, so this only maps
# them onto the columns the shared gene annotation step expects.
prepare_lung_variants <- function(dat) {
  dat <- data.table::as.data.table(dat)
  dat <- dat[!grepl(",", alt)]
  dat[, chr := sub("^chr", "", Chr)]
  dat <- dat[chr %in% as.character(1:22)]
  dat[, gene := sub("\\..*", "", gene)]
  dat <- dat[gene != "" & gene != "."]
  dat[, `:=`(
    chr = as.integer(chr),
    bp = as.integer(pos),
    gene_symbol = gene_name,
    p = as.numeric(pvalue),
    beta = as.numeric(beta),
    se = as.numeric(se)
  )]
  dat <- dat[!is.na(bp) & !is.na(beta) & !is.na(se) & !is.na(p)]
  return(dat)
}

join_gene_annotations <- function(dat, gene_map) {
  dat <- merge(dat, gene_map, by.x = "gene", by.y = "ENSEMBL_ID", all.x = TRUE)
  # DICE is cis, so every variant must sit on the gene's chromosome. The
  # hg19->hg38 chain maps terminal repeat regions across chromosomes (e.g. the
  # end of chr2 onto chr1/chr8), so drop anything on a different chromosome.
  dat[is.na(GeneChr), GeneChr := which.max(tabulate(chr)), by = gene]
  dat <- dat[chr == GeneChr]
  dat[, ProbeBp := BP_START]
  # Fall back to the most significant variant for genes absent from the map.
  dat[is.na(ProbeBp), ProbeBp := bp[which.min(p)], by = gene]
  dat[, Gene := ifelse(is.na(GENE_NAME) | GENE_NAME == "", gene_symbol, GENE_NAME)]
  dat[, Orientation := "+"]
  dat[, `:=`(BP_START = NULL, GENE_NAME = NULL, GeneChr = NULL)]
  return(dat)
}

write_esd_files <- function(dat, esd_out_dir) {
  dir.create(esd_out_dir, recursive = TRUE, showWarnings = FALSE)
  dat[, {
    data.table::fwrite(
      # liftover can collapse two source variants onto the same hg38 SNP; SMR
      # drops the duplicate, so drop it here (keeping the first) to match.
      .SD[!duplicated(SNP), .(Chr = CHR, SNP, Bp = BP, A1 = EA, A2 = OA, Freq = EAF, Beta = BETA, se = SE, p)],
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
  probes <- unique(dat[, .(Chr = CHR, ProbeID = gene, ProbeBp = ProbeBp, Gene = Gene, Orientation = Orientation)])
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

write_metadata_json <- function(metadata, json_file, tissue = "Whole Blood") {
  jsonlite::write_json(
    list(
      cell_type = metadata$cell_type,
      tissue = tissue,
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
  lung_files <- list.files(vcf_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
  if (length(args) > 0) {
    prefixes <- sub("\\.vcf\\.gz$", "", basename(vcf_files))
    vcf_files <- vcf_files[prefixes %in% args]
    lung_prefixes <- sub("\\.txt\\.gz$", "", basename(lung_files))
    lung_keys <- iconv(lung_prefixes, from = "UTF-8", to = "ASCII//TRANSLIT")
    lung_files <- lung_files[lung_prefixes %in% args | lung_keys %in% args]
  }
  if (length(vcf_files) == 0 && length(lung_files) == 0) {
    stop("No VCF or lung eQTL files found to process")
  }

  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(esd_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(flist_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(liftover_work_dir, recursive = TRUE, showWarnings = FALSE)
  write_contig_header(
    file.path(liftover_dir, "hg19.fa.fai"),
    file.path(liftover_work_dir, "contigs.txt")
  )

  gene_map <- data.table::fread(
    gene_map_file,
    select = c("ENSEMBL_ID", "GENE_NAME", "BP_START", "CHR"),
    showProgress = FALSE
  )
  gene_map <- unique(gene_map, by = "ENSEMBL_ID")
  data.table::setnames(gene_map, "CHR", "GeneChr")
  gene_map[, GeneChr := suppressWarnings(as.integer(GeneChr))]

  for (vcf_file in vcf_files) {
    prefix <- sub("\\.vcf\\.gz$", "", basename(vcf_file))
    metadata <- dice_metadata[dice_metadata$file_prefix == prefix, ]
    if (nrow(metadata) != 1) {
      stop(paste0("No cell type metadata defined for VCF prefix: ", prefix))
    }

    flist_file <- file.path(flist_dir, paste0(prefix, ".flist"))
    json_file <- file.path(flist_dir, paste0(prefix, ".json"))
    if (file.exists(flist_file) && file.exists(json_file)) {
      message(paste0("Skipping ", prefix, " (output already exists)"))
      next
    }

    message(paste0("Processing ", prefix))
    dat <- parse_lifted_vcf(vcf_file, prefix)
    dat <- prepare_variants(dat)
    dat <- join_gene_annotations(dat, gene_map)

    # Reuse the pipeline allele standardisation (EA/OA ordering, SNP naming).
    gwas <- data.frame(
      CHR = dat$chr,
      BP = dat$bp,
      EA = dat$alt,
      OA = dat$ref,
      EAF = NA_real_,
      BETA = dat$beta,
      SE = dat$se,
      gene = dat$gene,
      gene_symbol = dat$gene_symbol,
      p = dat$p,
      ProbeBp = dat$ProbeBp,
      Gene = dat$Gene,
      Orientation = dat$Orientation,
      stringsAsFactors = FALSE
    )
    rm(dat)
    gwas <- data.table::as.data.table(standardise_alleles(gwas))

    esd_out_dir <- file.path(esd_dir, prefix)
    write_esd_files(gwas, esd_out_dir)
    write_flist(gwas, flist_file, esd_out_dir)
    write_metadata_json(metadata, json_file)

    rm(gwas)
    gc()
  }

  for (lung_file in lung_files) {
    file_prefix <- sub("\\.txt\\.gz$", "", basename(lung_file))
    key <- iconv(file_prefix, from = "UTF-8", to = "ASCII//TRANSLIT")
    metadata <- dice_lung_metadata[dice_lung_metadata$file_prefix == key, ]
    if (nrow(metadata) != 1) {
      stop(paste0("No cell type metadata defined for lung file: ", file_prefix))
    }
    prefix <- paste0("DICE_LUNG_", gsub(" ", "_", metadata$cell_type))

    flist_file <- file.path(flist_dir, paste0(prefix, ".flist"))
    json_file <- file.path(flist_dir, paste0(prefix, ".json"))
    if (file.exists(flist_file) && file.exists(json_file)) {
      message(paste0("Skipping ", prefix, " (output already exists)"))
      next
    }

    message(paste0("Processing ", prefix))
    dat <- data.table::fread(
      cmd = paste0("zcat '", lung_file, "'"),
      select = c("gene", "Chr", "gene_name", "pos", "ref", "alt", "pvalue", "beta", "se"),
      showProgress = TRUE
    )
    dat <- prepare_lung_variants(dat)
    dat <- join_gene_annotations(dat, gene_map)

    # Reuse the pipeline allele standardisation (EA/OA ordering, SNP naming).
    gwas <- data.frame(
      CHR = dat$chr,
      BP = dat$bp,
      EA = dat$alt,
      OA = dat$ref,
      EAF = NA_real_,
      BETA = dat$beta,
      SE = dat$se,
      gene = dat$gene,
      gene_symbol = dat$gene_symbol,
      p = dat$p,
      ProbeBp = dat$ProbeBp,
      Gene = dat$Gene,
      Orientation = dat$Orientation,
      stringsAsFactors = FALSE
    )
    rm(dat)
    gwas <- data.table::as.data.table(standardise_alleles(gwas))

    esd_out_dir <- file.path(esd_dir, prefix)
    write_esd_files(gwas, esd_out_dir)
    write_flist(gwas, flist_file, esd_out_dir)
    write_metadata_json(metadata, json_file, tissue = "Lung")

    rm(gwas)
    gc()
  }

  message("Finished formatting DICE eQTLs. Run make_besd.sh to build BESD files.")
  return(invisible(NULL))
}

main()
