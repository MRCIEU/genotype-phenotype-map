source('core/common.R')

#TODO: this is probably useless now, and is done in the python imputation step.  Delete later if not needed.

parser <- argparser::arg_parser("Standardising GWAS for finemap")
parser <- argparser::add_argument(parser,
                                  "--list_of_studies",
                                  help = "GWAS filename",
                                  type = "character"
)
parser <- argparser::add_argument(parser,
                                  "--ld_region",
                                  help = "List of RSIDs to filter on",
                                  type = "character"
)

args <- argparser::parse_args(parser)
args <- list(list_of_studies=list('/local-scratch/projects/genotype-phenotype-map/data/study/ukb-b-10003/EUR_6_32530029.tsv'),
             ld_region='/Users/wt23152/Documents/Projects/scratch/011/data/ld_block_matrices/EUR/6_31571218_32682663.tsv.gz'
)

start <- Sys.time()
ld_region <- vroom::vroom(args$ld_region, col_names=F, show_col_types = F)
colnames(ld_region) <- ld_region$X1
ld_region <- ld_region[, 2:(ncol(ld_region)-1)]
print(Sys.time() - start)

studies <- c('/Users/wt23152/Documents/Projects/scratch/011/data/study/ukb-b-10003/EUR_6_32530029.tsv')
for (study in studies) {
  study <- '/Users/wt23152/Documents/Projects/scratch/011/data/study/ukb-b-10003/EUR_6_32530029.tsv'
  chr <- 6
  dir <- '/Users/wt23152/Documents/Projects/scratch/011/data/study/ukb-b-10003/imputed/'
  imputation_gwas_filename <- paste0(dir, 'z_', basename(study), '_chr', chr, '.txt')
  imputation_ld_filename <- paste0(dir, 'z_', basename(study), '_chr', chr, '.ld')

  extracted_gwas <- vroom::vroom(study, show_col_types = F) |>
    dplyr::filter(!duplicated(RSID))
  extracted_gwas$Z <- convert_beta_and_se_to_z_score(extracted_gwas$BETA, extracted_gwas$SE)
  extracted_gwas <- dplyr::select(extracted_gwas, RSID, BP, EA, OA, Z) |>
    dplyr::rename(rsid = RSID, pos = BP, A0 = EA, A1 = OA)
  vroom::vroom_write(extracted_gwas, imputation_gwas_filename, delim = '\t')

  keep <- colnames(ld_region) %in% extracted_gwas$rsid
  ld_for_gwas <- ld_region[keep, keep]
  vroom::vroom_write(ld_for_gwas, imputation_ld_filename, col_names = F)
}
