# Date: 10-10-2024
# Author: A.Hanson

# STEPS:
# 1. Identify GRCh38 ld regions containing rare variant top hits passing threshold (no clumping performed here)
# 2. Extract SNPs across region and store

pipeline_dir <- "/home/kf23130/Projects/021_GenoPhenoMap/genotype-phenotype-map" 
rawdata_dir <- "/local-scratch/data/ukb-seq"
output_dir <- "/local-scratch/projects/genotype-phenotype-map/data/study_sequencing"

metadata_file <- file.path(rawdata_dir,"downloads/azexwas/ukb-wes-az-studies.tsv")

parser <- argparser::arg_parser('Extract genomic regions from AstraZeneca PheWas portal summary statistics')
parser <- argparser::add_argument(parser, '--extracted_study_file', help = 'Summary stats file', type = 'character')
args <- argparser::parse_args(parser)

ld_regions <- vroom::vroom(file.path(pipeline_dir,"pipeline_steps/data/ld_regions_hg38_updated.tsv"), show_col_types = F)

# Add chromosome ends to ld_regions
ld_split <- split(ld_regions, paste0(ld_regions$chr, ld_regions$ancestry))
ld_split <- lapply(ld_split, function(x){
	chr_start <- x[1,] |> dplyr::mutate(end = start, start = 0)
	chr_end <- x[nrow(x),] |> dplyr::mutate(start = end, end = 250000000)
	chr_join <- rbind(chr_start, x, chr_end)
	return(chr_join)
})

ld_regions <- do.call("rbind", ld_split) |> dplyr::arrange(ancestry, chr, start)
ld_regions$string_region <- paste0(ld_regions$ancestry, "/", ld_regions$chr, "/", ld_regions$start, "_", ld_regions$end)

p_value_threshold <- 1.5e-4
maf_max_threshold <- 0.01
maf_min_threshold <- 0.000012 # Lower MAF theshold (for 10 expected observations of minor allele in 430,998 Europeans)
study_ancestry <- "EUR"

input_columns <- c("Variant","Chromosome","Base pair location","Effect allele","Other allele","AAF","Effect size","Effect size standard error","p-value")

main <- function(args) {
	
	message("\nInput file: ", args$extracted_study_file)
	metadata <- vroom::vroom(metadata_file, show_col_types = F) |>
		dplyr::filter(file_name == args$extracted_study_file)
	
	# Create output directory
	writeto <- file.path(output_dir, metadata$study_id)
	dir.create(writeto, showWarnings = F, recursive = T)
	message("Output dir: ", writeto) 

	# Read in summary statistics	
	study <- vroom::vroom(metadata$file_path, show_col_types = F)
	
	# Check columns and convert odds ratios
	study <- check_gwas(study)

	# Extract varaints passing theshold (no clumping) and write to file
	top_vars <- find_hits(study)

	if(nrow(top_vars) == 0) {
		stop(paste0(metadata$study_id, " (", metadata$accession, "): no top hits"))
	} else {
		vroom::vroom_write(top_vars["Variant"], file.path(writeto,"top_vars.txt"), col_names = F)
		message("Writing top hits to: ", writeto, "/top_vars.txt")
	}

	# Identify regions top variants fall in 
	vars_regions <- find_regions(top_vars, writeto)
	vroom::vroom_write(vars_regions, file.path(writeto,"extracted_snps.tsv"))
	message("Writing top regions to: ", writeto, "/extracted_snps.tsv")

	# Extact all variants in region and write to file
	dir.create(file.path(writeto, "original"), showWarnings = F, recursive = T)
	extract_regions(study, vars_regions)
}

check_gwas <- function(study){

	# Add missing columns
	study$Chromosome <- as.numeric(gsub("-.*","",study$Variant))
	study$`Base pair location` <- as.numeric(sub("^.*-(.*)-.*-.*","\\1",study$Variant))
	study$`Effect allele` <- sub(".*-(.*$)","\\1",study$Variant)
	study$`Other allele` <- sub(".*-(.*)-.*$","\\1",study$Variant)

	# Update variant ID
	study$Variant <- gsub("-",":",study$Variant)

	study_cols <- colnames(study)
	
	if(all(input_columns %in% study_cols)) {
		message("All required columns present")
	}else if("Odds ratio" %in% study_cols){
		message("Converting ORs to betas, calculating SE and AAF")
		study$`Odds ratio` <- ifelse(study$`Odds ratio` == 0, 0.01, study$`Odds ratio`)
		study$`Effect size` <- log(study$`Odds ratio`)
		UCI <-	log(study$`Odds ratio UCI`)
		study$`Effect size standard error` <- (UCI - study$`Effect size`)/1.96
		study$`Effect size standard error` <- ifelse(study$`Odds ratio` == 0.01, 0, study$`Effect size standard error`)
		study$AAF <- study$`Control AAF`
	}else{
		message("Cannot find effect estimate column")
	}
	return(study)
}

find_hits <- function(study) {
	study_filt <- study |>
		dplyr::filter(Chromosome %in% seq(1,22),
			`p-value` <= p_value_threshold,
			((AAF <= maf_max_threshold & AAF >= maf_min_threshold) |
			(1-AAF <= maf_max_threshold & 1-AAF >= maf_min_threshold)))
	study_filt <- study_filt |> dplyr::arrange(Chromosome, `Base pair location`)
	return(study_filt)
}

find_regions <- function(top_vars, writeto) {
	vars_regions <- data.frame(
		chr = numeric(),
		bp = numeric(),
		log_p = numeric(),
		ld_region = character(),
		file = character())
	
	vars_regions <- apply(top_vars, 1, function(variant){	

		extraction_info <- data.frame(
			chr = as.numeric(variant["Chromosome"]),
			bp = as.numeric(variant["Base pair location"]),
			log_p = -log10(as.numeric(variant["p-value"])),
			ld_region = ld_regions |> 
				dplyr::filter(chr == as.numeric(variant["Chromosome"]),
					start <= as.numeric(variant["Base pair location"]), 
					end > as.numeric(variant["Base pair location"]),
					ancestry == study_ancestry) |> 
				dplyr::pull(string_region))
			
		extraction_info$file = file.path(writeto, "original", paste0(study_ancestry, "_", extraction_info$chr, "_", extraction_info$bp, ".tsv.gz"))
		return(extraction_info)
	})

	vars_regions <- do.call("rbind", vars_regions)
	return(vars_regions)
}

extract_regions <- function(study, vars_regions){

	invisible(apply(vars_regions, 1, function(variant) {
		
		region_min <- as.numeric(strsplit(as.character(variant["ld_region"]), "_|/")[[1]][3])
		region_max <- as.numeric(strsplit(as.character(variant["ld_region"]), "_|/")[[1]][4])
		message("Extracting from region: ", variant["chr"], ":", region_min, "-", region_max) 

		study_region <- study |> dplyr::filter(Chromosome == as.numeric(variant["chr"]),
			((AAF <= maf_max_threshold & AAF >= maf_min_threshold) |
            (1-AAF <= maf_max_threshold & 1-AAF >= maf_min_threshold)),			
			`Base pair location` >= region_min,
			`Base pair location` < region_max)
		
		extracted_vars <- data.frame(
			NAME = study_region$Variant,
			CHR = study_region$Chromosome,
			BP = study_region$`Base pair location`,
			EA = study_region$`Effect allele`,
			OA = study_region$`Other allele`,
			EAF = study_region$AAF,
			BETA = study_region$`Effect size`,
			SE = study_region$`Effect size standard error`,
			LP = -log10(study_region$`p-value`))
		
		# Write out to region file
        vroom::vroom_write(extracted_vars, file.path(as.character(variant["file"])))	
	}))
}

main(args)
