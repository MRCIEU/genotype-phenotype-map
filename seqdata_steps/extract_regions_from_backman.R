# Date: 30-09-2024
# Author: A.Hanson

# STEPS:
# 1. Identify GRCh38 ld regions containing rare variant top hits passing threshold (no clumping performed here)
# 2. Extract SNPs across region and store

pipeline_dir <- "/home/kf23130/Projects/021_GenoPhenoMap/genotype-phenotype-map" 
rawdata_dir <- "/local-scratch/data/ukb-seq"
output_dir <- "/local-scratch/projects/genotype-phenotype-map/data/study_sequencing"

metadata_file <- file.path(rawdata_dir,"downloads/backmanexwas/ukb-wes-bm-studies.tsv")

parser <- argparser::arg_parser('Extract genomic regions from Backman et al. summary statistics')
parser <- argparser::add_argument(parser, '--extracted_study_file', help = 'Summary stats file', type = 'character')
args <- argparser::parse_args(parser)

ld_regions <- vroom::vroom(file.path(pipeline_dir,"pipeline_steps/data/ld_regions_hg38.tsv"), show_col_types = F)

p_value_threshold <- 1.5e-4
maf_max_threshold <- 0.01
maf_min_threshold <- 0.000012 # Lower MAF theshold (for 10 expected observations of minor allele in 430,998 Europeans)
study_ancestry <- "EUR"

input_columns <- c("Name","chromosome","base_pair_location","effect_allele","other_allele","beta","standard_error","effect_allele_frequency","p_value")

main <- function(args) {
	
	message("Input file: ", args$extracted_study_file)

	metadata <- vroom::vroom(metadata_file) |>
		dplyr::filter(file_name == args$extracted_study_file)
	
	# Create output directory
	writeto <- file.path(output_dir, metadata$study_id)
	dir.create(writeto, showWarnings = F, recursive = T)
	
	message("Output dir:", writeto) 

	# Read in summary statistics	
	study <- vroom::vroom(metadata$file_path, show_col_types = F)
	
	# Check columns and convert odds ratios
	study <- check_gwas(study)

	# Extract varaints passing theshold (no clumping) and write to file
	top_vars <- find_hits(study)
	vroom::vroom_write(top_vars$Name, file.path(writeto,"top_vars.txt"))

	# Identify regions top variants fall in 
	vars_regions <- find_regions(top_vars)
	vroom::vroom_write(vars_regions, file.path(writeto,"extracted_snps.tsv"))

	# Extact all variants in region and write to file
	extract_regions(study, vars_regions)
	
}

check_gwas <- function(study){
	study_cols <- colnames(study)

	if(all(input_columns %in% study_cols)) {
		message("All required columns present")
	}else if("odds_ratio" %in% study_cols){
		message("Converting ORs to betas")
		study$beta <- log(study$odds_ratio)
	}else{
		message("Cannot find effect estimate column")
	}

	return(study)
}

find_hits <- function(study) {
	study_filt <- study |>
		dplyr::filter(chromosome %in% seq(1,22),
			      p_value <= p_value_threshold,
			     (effect_allele_frequency <= maf_max_threshold | 1-effect_allele_frequency <= maf_max_threshold),
			     (effect_allele_frequency >= maf_min_threshold | 1-effect_allele_frequency >= maf_min_threshold)) 	
	return(study_filt)
}

find_regions <- function(top_vars) {
	vars_regions <- data.frame(
		chr = numeric(),
		bp = numeric(),
		log_p = numeric(),
		ld_region = character(),
		file = character()
	)

	if(length(top_vars) == 0) {
		message(paste0(metadata$studyid, " (", metadata$accession, "): no top hits"))
		return(vars_regions)
	}

	ld_regions$string_region <- paste0(ld_regions$ancestry, "/", ld_regions$chr, "/", ld_regions$start, "_", ld_regions$end)

	vars_regions <- apply(top_vars, 1, function(variant){
		extraction_info <- data.frame(
			chr = as.numeric(variant["chromosome"]),
			bp = as.numeric(variant["base_pair_location"]),
			log_p = -log10(as.numeric(variant["p_value"])),
			ld_region = ld_regions |> dplyr::filter(chr == as.numeric(variant["chromosome"]), 
								start <= as.numeric(variant["base_pair_location"]), 
								end > as.numeric(variant["base_pair_location"]),
								ancestry == study_ancestry) |> pull(string_region)
		      )
		extraction_info$file = file.path(writeto, "original", paste0(ancestry, "_", extraction_info$chr, "_", extraction_info$bp, ".tsv.gz"))
		
		return(extraction_info)
	})

	vars_regions <- do.call("rbind", vars_regions)

	return(vars_regions)
}

extract_regions <- function(study, vars_regions){
	
	apply(vars_regions, 1, function(variant) {
		
		region_min <- as.numeric(strsplit(as.character(variant["ld_region"]), "_|/")[[1]][3])
		region_max <- as.numeric(strsplit(as.character(variant["ld_region"]), "_|/")[[1]][4])

		study_region <- study |>
                	dplyr::filter(chromosome == as.numeric(variant["chr"]),
                             (effect_allele_frequency <= maf_max_threshold | 1-effect_allele_frequency <= maf_max_threshold),
                             (effect_allele_frequency >= maf_min_threshold | 1-effect_allele_frequency >= maf_min_threshold),
			     base_pair_location >= region_min,
			     base_pair_location < region_max)

		extracted_vars <- data.frame(
			NAME = study_region$Name,
			CHR = study_region$chromosome,
			BP = study_region$base_pair_location,
			EA = study_region$effect_allele,
			OA = study_region$other_allele,
			EAF = study_region$effect_allele_frequency,
			BETA = study_region$beta,
			SE = study_region$standard_error,
			LP = -log10(study_region$p_value))
 
		# Write out to region file
		vroom::vroom_write(extracted_vars, as.character(variant["file"])) 
	})
}

main(args)
