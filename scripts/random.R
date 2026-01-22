source('../pipeline_steps/constants.R')
imputed <- vroom::vroom('/local-scratch/projects/genotype-phenotype-map/data/study/ebi-a-GCST90012877/imputed/EUR_1_1053552.tsv.gz')

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

# Coerce IMPUTED to logical or factor for coloring
df$IMPUTED <- as.factor(df$IMPUTED)

# Manhattan plot for chromosome 1 only
ggplot(df, aes(x = BP, y = -log10(P), color = IMPUTED)) +
  geom_point(alpha = 0.7, size = 1) +
  labs(
    title = "Manhattan plot: Chromosome 1",
    x = "Genomic Position",
    y = expression(-log[10](italic(P))),
    color = "Imputed"
  ) +
  theme_minimal()
