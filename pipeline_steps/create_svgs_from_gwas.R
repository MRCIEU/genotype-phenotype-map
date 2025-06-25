create_svg_for_ld_block <- function(study) {
  dir.create(glue::glue('{study$extracted_location}/svgs/extractions'), showWarnings = F, recursive = T)

  plot_height <- 200
  plot_width <- 1000
  
  ld_block <- vroom::vroom(study[['file']]) |>
    dplyr::mutate(CHR = as.numeric(CHR), Z = Z) |>
    dplyr::filter(!is.na(CHR)) |>
    dplyr::arrange(CHR, BP)
  
  ld_data <- ld_block |>
    dplyr::filter(!is.na(Z)) |>
    dplyr::mutate(bin = floor(BP/10000)) |>
    dplyr::group_by(CHR, bin) |>
    dplyr::summarise(
      BP = mean(BP),
      Z = max(Z),
      .groups = "drop"
    )
  
  ld_area <- ld_data |>
    dplyr::mutate(bin = floor(BP/10000)) |>
    dplyr::group_by(CHR, bin) |>
    dplyr::summarise(
      BP = mean(BP),
      Z = max(Z),
      .groups = "drop"
    )

  # Add start and end rows with Z = 0 for each chromosome
  ld_area <- ld_area |>
    dplyr::group_by(CHR) |>
    dplyr::arrange(BP, .by_group = TRUE) |>
    dplyr::group_modify(~{
      min_bp <- min(.x$BP)
      max_bp <- max(.x$BP)
      
      dplyr::bind_rows(
        dplyr::tibble(CHR = .x$CHR[1], BP = min_bp, Z = 0),
        .x,
        dplyr::tibble(CHR = .x$CHR[1], BP = max_bp, Z = 0)
      )
    }) |>
    dplyr::ungroup()
  
  p <- ggplot2::ggplot(ld_data, ggplot2::aes(x = BP, y = Z)) +
    ggplot2::geom_area(
      data = ld_area,
      ggplot2::aes(group = CHR),
      position = "identity"
    ) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max(10, max(ld_data$Z, na.rm=TRUE)))) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      plot.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(0, "lines"),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    ) +
    ggplot2::coord_cartesian(clip = "off")
  
  block_name <- basename(study[['file']]) |> stringr::str_replace("\\.tsv\\.gz$", "")
  
  file_name <- glue::glue("{data_dir}{study$study}/svgs/extractions/{block_name}.svg")
  ggplot2::ggsave(
    file_name,
    p,
    width = plot_width/72,
    height = plot_height/72,
    units = "in",
    limitsize = FALSE
  )
}

create_svgs_from_gwas <- function(study, gwas) {
  plot_height <- 500
  plot_width <- 1250

  original_gwas <- gwas |>
    dplyr::mutate(CHR = as.numeric(CHR)) |>
    dplyr::filter(!is.na(CHR)) |>
    dplyr::arrange(CHR, BP)

  chr_ranges <- original_gwas |>
    dplyr::group_by(CHR) |>
    dplyr::summarise(
      bp_min = min(BP),
      bp_max = max(BP),
      .groups = "drop"
  ) |>
  dplyr::mutate(
      chr_size = bp_max - bp_min,
      # Adding a small buffer between chromosomes (1 unit for example)
      bp_cumulative_start = cumsum(c(0, chr_size[-dplyr::n()] + 1)) +
      (1:dplyr::n() - 1), # Additional 1 for separation if needed
      bp_cumulative_end = bp_cumulative_start + chr_size
  )

  gwas <- original_gwas |>
    dplyr::filter(!is.na(LP)) |>
    dplyr::mutate(bin = floor(BP/100000)) |>
    dplyr::group_by(CHR, bin) |>
    dplyr::summarise(
      BP = mean(BP),
      LP = max(LP),
      .groups = "drop"
    ) |>
    dplyr::arrange(CHR, BP) |>
    dplyr::left_join(chr_ranges |> dplyr::select(CHR, bp_cumulative_start, bp_min), by = "CHR") |>
    dplyr::mutate(BP_cumulative = BP - bp_min + bp_cumulative_start) |>
    # Add a flag for non-adjacent points
    dplyr::mutate(
      is_adjacent = CHR == dplyr::lead(CHR) &
       abs(BP - dplyr::lead(BP)) < 100000 * 2 # Allow for 2x bin size gap
    )

  # Create a second dataset with larger bins for the area
  gwas_area <- gwas |>
    dplyr::mutate(bin = floor(BP/1000000)) |> # Create 1Mb bins for the area
    dplyr::group_by(CHR, bin) |>
    dplyr::summarise(
      BP = mean(BP),
      LP = max(LP),
      .groups = "drop"
    ) |>
    dplyr::arrange(CHR, BP) |>
    dplyr::left_join(chr_ranges |> dplyr::select(CHR, bp_cumulative_start, bp_min), by = "CHR") |>
    dplyr::mutate(BP_cumulative = BP - bp_min + bp_cumulative_start)

  # Add start and end rows with Z = 0 for each chromosome
  gwas_area <- gwas_area |>
    dplyr::group_by(CHR) |>
    dplyr::arrange(BP_cumulative, .by_group = TRUE) |>
    dplyr::group_modify(~{
      min_bp_cum <- min(.x$BP_cumulative)
      max_bp_cum <- max(.x$BP_cumulative)

      dplyr::bind_rows(
      dplyr::tibble(
          CHR = .x$CHR[1],
          BP = min(.x$BP),
          LP = 0,
          BP_cumulative = min_bp_cum
      ),
      .x,
      dplyr::tibble(
          CHR = .x$CHR[1],
          BP = max(.x$BP),
          LP = 0,
          BP_cumulative = max_bp_cum
      )
      )
  }) |>
  dplyr::ungroup()

  # Create optimized Manhattan plot with lines
  p <- ggplot2::ggplot(mapping = ggplot2::aes(x = BP_cumulative, y = LP)) +
    # Add area with larger bins
    ggplot2::geom_area(
      data = gwas_area,
      ggplot2::aes(group = CHR, fill = "#7F7F7F"), # Removed 'color = color' from here
      position = "identity",
      alpha = 1 # Keep alpha at 1 for solid fill
    ) +
    # Add detailed lines
    ggplot2::geom_line(data = gwas, ggplot2::aes(group = CHR, color = "#7F7F7F"), linewidth = 0.8) +
    ggplot2::scale_color_manual(values = c("#7F7F7F")) +
    ggplot2::scale_fill_manual(values = c("#7F7F7F")) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max(10, max(gwas$LP, na.rm=TRUE)))) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = "none",
      plot.margin = ggplot2::margin(0, 0, 0, 0),
      plot.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.spacing = ggplot2::unit(0, "lines"),
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank()
    ) +
    ggplot2::coord_cartesian(clip = "off")


  ggplot2::ggsave(
    glue::glue("{study$extracted_location}/svgs/full.svg"),
    p,
    width = plot_width/72,  # Convert to inches for ggsave
    height = plot_height/72,  # Convert to inches for ggsave
    units = "in",  # Specify inches
    limitsize = FALSE
  )

  for (chr_num in unique(original_gwas$CHR)) {
    chr_data <- gwas |>
      dplyr::filter(CHR == chr_num)

    chr_area_data <- gwas_area |>
      dplyr::filter(CHR == chr_num)

    p_chr <- ggplot2::ggplot(mapping = ggplot2::aes(x = BP_cumulative, y = LP)) +
      ggplot2::geom_area(
        data = chr_area_data,
        ggplot2::aes(group = CHR, fill = "#7F7F7F", color = "#7F7F7F"),
        position = "identity",
        alpha = 1
      ) +
      ggplot2::geom_line(data = chr_data, ggplot2::aes(group = CHR, color = "#7F7F7F"), linewidth = 0.8) +
      ggplot2::scale_color_manual(values = c("#7F7F7F")) +
      ggplot2::scale_fill_manual(values = c("#7F7F7F")) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max(10, max(gwas$LP, na.rm=TRUE)))) +
      ggplot2::theme_void() +
      ggplot2::theme(
        legend.position = "none",
        plot.margin = ggplot2::margin(0, 0, 0, 0),
        plot.background = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        panel.spacing = ggplot2::unit(0, "lines"),
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        axis.text = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank()
      ) +
      ggplot2::coord_cartesian(clip = "off")

    svg_location <- glue::glue('{study$extracted_location}/svgs/chr{chr_num}.svg')
    ggplot2::ggsave(
      svg_location,
      p_chr,
      width = plot_width/72,
      height = plot_height/72,
      units = "in",
      limitsize = FALSE
    )
  }
  old_wd <- getwd()
  setwd(glue::glue('{study$extracted_location}/svgs'))
  files_to_zip <- Sys.glob('*.svg')
  utils::zip('full.zip', files_to_zip)
  file.remove(files_to_zip)
  setwd(old_wd)

  total_cumulative_bp <- max(chr_ranges$bp_cumulative_end)

  x_axis_metadata <- chr_ranges |>
    dplyr::mutate(
      pixel_start = bp_to_pixel(bp_cumulative_start, total_cumulative_bp, plot_width),
      pixel_end = bp_to_pixel(bp_cumulative_end, total_cumulative_bp, plot_width)
    ) |>
    dplyr::select(CHR, bp_start = bp_cumulative_start, bp_end = bp_cumulative_end, pixel_start, pixel_end)

  metadata <- list(
    x_axis = x_axis_metadata,
    y_axis = list(
      max_lp = max(10, max(gwas$LP)),
      min_lp = min(gwas$LP),
      min_lp_pixel = lp_to_pixel(min(gwas$LP), min(gwas$LP), max(10, max(gwas$LP)), plot_height), 
      max_lp_pixel = lp_to_pixel(max(10, max(gwas$LP)), min(gwas$LP), max(10, max(gwas$LP)), plot_height)
    ),
    svg_width = plot_width,
    svg_height = plot_height
  )

  jsonlite::write_json(metadata, glue::glue("{study$extracted_location}/svgs/full.json"), pretty = TRUE, auto_unbox = TRUE)
}

bp_to_pixel <- function(cumulative_bp, total_cumulative_bp, svg_width = 1250) {
  pixel <- (cumulative_bp / total_cumulative_bp) * svg_width
  return(pixel)
}

lp_to_pixel <- function(lp, min_lp, max_lp, height_pixels = 500, margin_ratio = 0.0) {
  position_ratio <- (lp - min_lp) / (max_lp - min_lp)
  margin <- height_pixels * margin_ratio
  plot_height <- height_pixels - (2 * margin)

  return(height_pixels - margin - (position_ratio * plot_height))
}
