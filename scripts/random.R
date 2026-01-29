source('../pipeline_steps/constants.R')
source('../pipeline_steps/gwas_calculations.R')
study <- 'ebi-a-GCST90014003'

imputed <- vroom::vroom('/local-scratch/projects/genotype-phenotype-map/data/study/ebi-a-GCST90014003/imputed/EUR_1_923421.tsv.gz')

library(ggplot2)

# Assuming 'imputed' columns include at least: CHR, POS, P, IMPUTED
# If column names are lowercase, adjust accordingly

# Take just the first chromosome segment (assuming 'CHR' or 'chr' column)
if ("CHR" %in% colnames(imputed)) {
  chr_col <- "CHR"
} else if ("chr" %in% colnames(imputed)) {
  chr_col <- "chr"
} else {
  stop("Cannot find chromosome column in data.")
}

# For plotting, P must be numeric, and IMPUTED must be logical or factor
df <- imputed[as.character(imputed[[chr_col]]) == "1", ]
if (!"IMPUTED" %in% names(df)) stop("IMPUTED column not found in data")

# Add about 200 more random IMPUTED = TRUE rows
imputed_true_rows <- df[df$IMPUTED == TRUE | df$IMPUTED == "TRUE" | df$IMPUTED == 1, ]
if (nrow(imputed_true_rows) > 0) {
  n_to_add <- min(400, nrow(imputed_true_rows))
  random_indices <- sample(nrow(imputed_true_rows), n_to_add)
  additional_rows <- imputed_true_rows[random_indices, ]
  df <- dplyr::bind_rows(df, additional_rows)
}

# Coerce IMPUTED to logical or factor for coloring
df$IMPUTED <- as.factor(df$IMPUTED)
cb_blue <- "#0072B2"
cb_orange <- "#E69F00"

# Manhattan plot for chromosome 1 only
p <- ggplot(df, aes(x = BP, y = -log10(P), color = IMPUTED)) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_manual(values = c(cb_blue, cb_orange)) +
  labs(
    title = "Manhattan plot: Chromosome 1",
    x = "Genomic Position",
    y = expression(-log[10](italic(P))),
    color = "Imputed"
  ) +
  theme_minimal()

ggsave('imputed_manhattan_plot.png', p, width = 12, height = 4, units = "in")

study <- 'ebi-a-GCST90014003'
# study <- 'ebi-a-GCST90025965'
finemapped_1 <- vroom::vroom(glue::glue('/local-scratch/projects/genotype-phenotype-map/data/study/{study}/finemapped/EUR_1_16103_1170341_1.tsv.gz'))
finemapped_2 <- vroom::vroom(glue::glue('/local-scratch/projects/genotype-phenotype-map/data/study/{study}/finemapped/EUR_1_16103_1170341_2.tsv.gz'))

finemapped_1 <- dplyr::mutate(finemapped_1, Z = convert_lbf_to_abs_z(LBF, SE), dataset = "Finemap 1")
finemapped_2 <- dplyr::mutate(finemapped_2, Z = convert_lbf_to_abs_z(LBF, SE), dataset = "Finemap 2")

# Combine the two datasets
finemapped_combined <- dplyr::bind_rows(finemapped_1, finemapped_2)

# Manhattan plot with both datasets overlaid
p2 <- ggplot(finemapped_combined, aes(x = BP, y = Z, color = dataset)) +
  geom_point(alpha = 0.7, size = 1) +
  scale_color_manual(values = c(cb_blue, cb_orange)) +
  labs(
    title = "Manhattan plot: Finemapped GWAS Comparison",
    x = "Genomic Position",
    y = "Absolute Z-score",
    color = "Dataset"
  ) +
  theme_minimal()

ggsave('finemapped_manhattan_plot.png', p2, width = 12, height = 4, units = "in")