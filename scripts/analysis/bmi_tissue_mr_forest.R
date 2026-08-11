library(ggforestplot)
library(tidyverse)
library(data.table)
library(ggpubr)
library(ggforce)
library(RColorBrewer)
library(rlang)
library(patchwork)

forestplot1 <- function(df, name = name, estimate = estimate, se = se, pvalue = NULL,
                        colour = NULL, shape = NULL, logodds = FALSE, psignif = 0.05,
                        ci = 0.95, ...) {
  stopifnot(is.data.frame(df))
  stopifnot(is.logical(logodds))
  name <- enquo(name)
  estimate <- enquo(estimate)
  se <- enquo(se)
  pvalue <- enquo(pvalue)
  colour <- enquo(colour)
  shape <- enquo(shape)
  const <- stats::qnorm(1 - (1 - ci) / 2)
  df <- df %>% dplyr::mutate(
    `:=`(!!name, factor(!!name, levels = !!name %>% unique() %>% rev(), ordered = TRUE)),
    .xmin = !!estimate - const * !!se,
    .xmax = !!estimate + const * !!se,
    .filled = TRUE,
    .label = sprintf("%.2f", !!estimate)
  )
  if (logodds) {
    df <- df %>% mutate(
      .xmin = exp(.data$.xmin),
      .xmax = exp(.data$.xmax),
      `:=`(!!estimate, exp(!!estimate))
    )
  }
  if (!quo_is_null(pvalue)) {
    df <- df %>% dplyr::mutate(.filled = !!pvalue < !!psignif)
  }
  g <- ggplot2::ggplot(df, aes(x = !!estimate, y = !!name))
  if (logodds) {
    g <- g + scale_x_continuous(trans = "log10", breaks = scales::log_breaks(n = 7))
  }
  g <- g + theme_forest() + scale_colour_ng_d() + scale_fill_ng_d() +
    geom_vline(
      xintercept = ifelse(test = logodds, yes = 1, no = 0),
      linetype = "solid", size = 0.4, colour = "black"
    )
  g <- g + geom_effect(
    ggplot2::aes(
      xmin = .data$.xmin, xmax = .data$.xmax,
      colour = !!colour, shape = !!shape, filled = .data$.filled
    ),
    position = ggstance::position_dodgev(height = 0.5)
  ) +
    ggplot2::scale_shape_manual(values = c(23L, 21L, 22L, 24L, 25L)) +
    guides(
      colour = guide_legend(reverse = TRUE),
      shape = guide_legend(reverse = TRUE)
    )
  args <- list(...)
  if ("title" %in% names(args)) {
    g <- g + labs(title = args$title)
  }
  if ("subtitle" %in% names(args)) {
    g <- g + labs(subtitle = args$subtitle)
  }
  if ("caption" %in% names(args)) {
    g <- g + labs(caption = args$caption)
  }
  if ("xlab" %in% names(args)) {
    g <- g + labs(x = args$xlab)
  }
  if (!"ylab" %in% names(args)) {
    args$ylab <- ""
  }
  g <- g + labs(y = args$ylab)
  if ("xlim" %in% names(args)) {
    g <- g + coord_cartesian(xlim = args$xlim)
  }
  if ("ylim" %in% names(args)) {
    g <- g + ylim(args$ylim)
  }
  return(g)
}


#####
source("../../pipeline_steps/constants.R")
out_dir <- if (nzchar(results_dir)) {
  file.path(results_dir, "1.0.0", "analysis")
} else {
  getwd()
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

in_dir <- if (any(grepl("results\\.txt$", list.files(out_dir)))) out_dir else "."
files <- grep("results\\.txt$", list.files(in_dir), value = TRUE)

#############

file <- files[1]
message("Processing: ", basename(file))

# read and clean
input <- fread(file.path(in_dir, file), header = TRUE)
input$outcome <- sub(" \\|\\|.*", "", input$outcome)
input1 <- as_tibble(input)

## trying this
trait1 <- "Brain (Meta)"

# ensure exposure is character first
input1 <- input1 %>%
  mutate(exposure = as.character(exposure))

# set factor ordering
input1 <- input1 %>%
  mutate(
    exposure = factor(
      exposure,
      levels = c(trait1, setdiff(sort(unique(exposure)), trait1))
    )
  )

# plot
p <- forestplot1(
  df = input1,
  name = outcome,
  estimate = b,
  logodds = TRUE,
  colour = exposure,
  pvalue = pval,
  psignif = 0.05,
  title = "TPMR estimates for BMI\nBrain (Meta) vs Adipose Subcutaneous",
  xlab = "Odds ratio per 1-SD change in risk factor (95% CI)"
) +
  scale_color_manual(
    values = setNames(
      c("#E69F00", "#0072B2"),
      c(trait1, setdiff(sort(unique(input1$exposure)), trait1))
    )
  ) +
  guides(shape = guide_legend(override.aes = list(linetype = "blank"))) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    text = element_text(colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(colour = "black"),
    panel.border = element_rect(fill = NA, colour = "grey20", size = rel(1)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    axis.ticks = element_line(size = 0.25, colour = "grey20")
  )

output_file <- file.path(out_dir, "tpmr_brain_adipose.tiff")
ggsave(
  filename = output_file,
  plot = p,
  width = 6.5,
  height = 5,
  units = "in",
  dpi = 300,
  create.dir = TRUE
)
message("Saved forest plot to: ", output_file)
