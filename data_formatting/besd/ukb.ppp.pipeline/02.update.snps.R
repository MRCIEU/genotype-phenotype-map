library(data.table)
library(dplyr)
library(future)
library(furrr)

###setwd("/home/rj18633/scratch/gp.map/data/besd.formatting/test.out/v1/v3/processed")
###out.path <- "/home/rj18633/scratch/gp.map/data/besd.formatting/test.out/v1/v3/snps.updated"
 

# Use environment variables for directories
processed_dir <- Sys.getenv("PROCESSED_DIR")
setwd(processed_dir)

out.path <- Sys.getenv("SNPS_UPDATED_DIR")
dir.create(out.path, recursive = TRUE, showWarnings = FALSE)


# Set up parallel session
num_workers <- min(8, availableCores())  
plan(multisession, workers = num_workers)

sub_dirs <- list.dirs(path = ".", recursive = FALSE, full.names = TRUE)




# Function to standardize alleles
standardise_alleles <- function(pwas) {
  pwas$Original <- TRUE  # Add a column to track original rows
  
  # Identify rows to flip
  to_flip <- (pwas$A1 > pwas$A2)
  if (any(to_flip)) {
    pwas$Freq[to_flip] <- 1 - pwas$Freq[to_flip]
    pwas$Beta[to_flip] <- -1 * pwas$Beta[to_flip]
    
    temp <- pwas$A2[to_flip]
    pwas$A2[to_flip] <- pwas$A1[to_flip]
    pwas$A1[to_flip] <- temp
    
    pwas$Original[to_flip] <- FALSE  # Mark these rows as flipped
  }
  
  # Normalize SNP identifiers
  pwas$SNP <- toupper(paste0(pwas$SNP, "_", pwas$A1, "_", pwas$A2))
  # Remove duplicates, keeping non-flipped rows where available
  pwas <- pwas %>%
    dplyr::arrange(SNP, desc(Original)) %>%  # Prefer original rows
    dplyr::distinct(SNP, .keep_all = TRUE)  # Keep only the first occurrence
  
  pwas <- pwas %>% dplyr::select(-Original)  # Drop the Original column
  pwas <- dplyr::select(pwas, Chr, SNP, Bp, A1, A2, Freq, Beta, se, dplyr::everything())
  return(pwas)
}



process_file <- function(sub_dir) {
  gz_files <- list.files(path = sub_dir, pattern = "\\.gz$", full.names = TRUE)
  
  if (length(gz_files) > 0) {
    # Get the name of the subdirectory to create in the output path
    sub_dir_name <- basename(sub_dir)
    
    # Create a new subdirectory within the output path
    output_subdir <- file.path(out.path, sub_dir_name)
    if (!dir.exists(output_subdir)) {
      dir.create(output_subdir, recursive = TRUE)
    }
    
    # Process each .gz file in the subdirectory
    for (gz_file in gz_files) {
      data <- fread(gz_file)
      colnames(data) <- c("Chr", "SNP", "Bp", "A1", "A2", "Freq", "Beta", "se", "p")
      
      standardized_data <- standardise_alleles(data)
      
      # Construct the output file path
      output_file <- file.path(output_subdir, paste0(basename(gz_file)))
      
      fwrite(standardized_data, output_file, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
      
      message(paste("Processed:", gz_file, "->", output_file))
    }
  } else {
    message(paste("No .gz file found in", sub_dir))
  }
}
 

# process in parallel
future_map(sub_dirs, process_file)