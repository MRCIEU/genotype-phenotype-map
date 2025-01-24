# STEPS:
# 1. Identify GRCh38 ld regions containing rare variant top hits passing threshold (no clumping performed here)
# 2. Extract SNPs across region and store

pipeline_dir <- "/home/kf23130/Projects/021_GenoPhenoMap/genotype-phenotype-map" 

# Input study_location and extracted_location from study metadata file
parser <- argparser::arg_parser('Extract genomic regions from sequencing GWAS summary statistics')
parser <- argparser::add_argument(parser, '--extracted_study_file', help = 'Summary stats file', type = 'character')
parser <- argparser::add_argument(parser, '--extracted_output_location', help = 'Extracted output location', type = 'character')
args <- argparser::parse_args(parser)

ld_regions <- vroom::vroom(file.path(pipeline_dir,"pipeline_steps/data/ld_blocks.tsv"), show_col_types = F)
ld_regions$string_region <- paste0(ld_regions$ancestry, "/", ld_regions$chr, "/", ld_regions$start, "_", ld_regions$stop)

p_value_threshold <- 1.5e-4
maf_max_threshold <- 0.01
maf_min_threshold <- 0.000012 # Lower MAF theshold (for 10 expected observations of minor allele in 431,000 Europeans (approx n UKB))
study_ancestry <- "EUR"

required_columns <- c("CHR","BP","EA","OA","EAF","BETA","SE","P")

main <- function(args) {
	input_file <- basename(args$extracted_study_file)
    accession <- basename(args$extracted_output_location)

	message("\nInput file: ", input_file)
	
	# Create output directory
	writeto <- args$extracted_output_location
	dir.create(writeto, showWarnings = F, recursive = T)
	message("Output dir: ", writeto) 

	# Read in summary statistics	
	study <- vroom::vroom(args$extracted_study_file, show_col_types = F)

	# Check columns and convert odds ratios
	study <- check_gwas(study)

    # Add formatted variant ID
    study$NAME <- paste(study$CHR, study$BP, study$OA, study$EA, sep = ":")

	# Extract varaints passing theshold (no clumping) and write to file
	top_vars <- find_hits(study)

	if(nrow(top_vars) == 0) {
		stop(paste0(input_file, " (", accession, "): no top hits"))
	} else {
		vroom::vroom_write(top_vars["NAME"], file.path(writeto,"top_vars.txt"), col_names = F)
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

	study_cols <- colnames(study)
	
	if(all(required_columns %in% study_cols)) {
		message("All required columns present")
	}else if(all(c("OR","CI_UPPER") %in% study_cols & all(required_columns[-6] %in% study_cols))){
		message("Converting OR to BETA, calculating SE from CIs")
		study$OR <- ifelse(study$OR == 0, 0.01, study$OR)
		study$BETA <- log(study$OR)
		UCI <-	log(study$CI_UPPER)
		study$SE <- (UCI - study$BETA)/1.96
		study$SE <- ifelse(study$OR == 0.01, 0, study$SE)
	}else{
		stop("Cannot find required columns: \n study must have BETA and SE or OR and CI_UPPER/CI_LOWER plus CHR, BP, EA, OA, EAF, P")
	}
	return(study)
}

find_hits <- function(study) {

    study_cols <- colnames(study)

    if(!all(c("CHR","P","EAF") %in% study_cols)){
        stop("Missing CHR, P or EAF column")
    }
    else{
        study_filt <- study |>
		    dplyr::filter(CHR %in% seq(1,22),
			    P <= p_value_threshold,
			    ((EAF <= maf_max_threshold & EAF >= maf_min_threshold) |
			    (1-EAF <= maf_max_threshold & 1-EAF >= maf_min_threshold)))
	    study_filt <- study_filt |> dplyr::arrange(CHR, BP)
    }
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
			chr = as.numeric(variant["CHR"]),
			bp = as.numeric(variant["BP"]),
			log_p = -log10(as.numeric(variant["P"])),
			ld_region = ld_regions |> 
				dplyr::filter(chr == as.numeric(variant["CHR"]),
					start <= as.numeric(variant["BP"]), 
					end > as.numeric(variant["BP"]),
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

		study_region <- study |> dplyr::filter(CHR == as.numeric(variant["chr"]),
			((EAF <= maf_max_threshold & EAF >= maf_min_threshold) |
            (1-EAF <= maf_max_threshold & 1-EAF >= maf_min_threshold)),			
			BP >= region_min,
			BP < region_max)
		
		extracted_vars <- data.frame(
			NAME = study_region$NAME,
			CHR = study_region$CHR,
			BP = study_region$BP,
			EA = study_region$EA,
			OA = study_region$OA,
			EAF = study_region$EAF,
			BETA = study_region$BETA,
			SE = study_region$SE,
			LP = -log10(study_region$P))
		
		# Write out to region file
        vroom::vroom_write(extracted_vars, file.path(as.character(variant["file"])))	
	}))
}

main(args)