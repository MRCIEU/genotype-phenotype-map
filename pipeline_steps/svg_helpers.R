source("../pipeline_steps/gwas_calculations.R")

create_svg_for_ld_block <- function(gwas, file_name, ld_block, is_sparse = FALSE) {
  # if (file.exists(file_name)) return()
  if (is_sparse || nrow(gwas) < 1000) {
    bin_size <- 2000
  } else {
    bin_size <- 20000
  }

  plot_height <- 200
  plot_width <- 1000

  # TODO: change this to just data/ld_blocks.tsv after backfill
  ld_blocks <- vroom::vroom(glue::glue("../pipeline_steps/data/ld_blocks.tsv"), show_col_types = F)
  ld_info <- construct_ld_block(ld_blocks$ancestry, ld_blocks$chr, ld_blocks$start, ld_blocks$stop)
  specific_ld_info <- ld_info[ld_info$block == ld_block, ]
  start_bp <- specific_ld_info$start
  end_bp <- specific_ld_info$stop

  ld_block <- gwas |>
    dplyr::mutate(CHR = as.numeric(CHR)) |>
    dplyr::filter(!is.na(CHR)) |>
    dplyr::arrange(CHR, BP)

  if (!"Z" %in% names(ld_block)) {
    ld_block <- dplyr::mutate(ld_block, Z = convert_lbf_to_abs_z(LBF, SE))
  }
  ld_block$Z <- abs(ld_block$Z)
  ld_block <- dplyr::mutate(ld_block, Z = ifelse(Z == 0, 0.001, Z))

  ld_data <- ld_block |>
    dplyr::filter(!is.na(Z)) |>
    dplyr::mutate(bin = floor(BP / bin_size)) |>
    dplyr::group_by(CHR, bin) |>
    dplyr::summarise(
      BP = mean(BP),
      Z = max(Z),
      .groups = "drop"
    )

  ld_area <- ld_data |>
    dplyr::mutate(bin = floor(BP / bin_size)) |>
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
    dplyr::group_modify(~ {
      min_bp <- min(.x$BP)
      max_bp <- max(.x$BP)

      dplyr::bind_rows(
        dplyr::tibble(CHR = .x$CHR[1], BP = start_bp, Z = 0),
        dplyr::tibble(CHR = .x$CHR[1], BP = min_bp, Z = 0.001),
        .x,
        dplyr::tibble(CHR = .x$CHR[1], BP = max_bp, Z = 0.001),
        dplyr::tibble(CHR = .x$CHR[1], BP = end_bp, Z = 0)
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
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max(10, max(ld_data$Z, na.rm = TRUE)))) +
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

  dir.create(dirname(file_name), showWarnings = F, recursive = T)
  ggplot2::ggsave(
    file_name,
    p,
    width = plot_width / 72,
    height = plot_height / 72,
    units = "in",
    limitsize = FALSE
  )
  return()
}

create_svgs_from_gwas <- function(study, gwas) {
  plot_height <- 500
  plot_width <- 1250

  max_lp_y_axis <- min(300, max(10, max(gwas$LP, na.rm = TRUE)))

  # Pad GWAS with BP 1 at the start of each chromosome, so we get consistent CHR sizes
  padded_chrs <- gwas[1, ][rep(1, 22), ]
  padded_chrs$CHR <- 1:22
  padded_chrs$BP <- 1
  padded_chrs$LP <- 0

  gwas <- rbind(padded_chrs, gwas)

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
      bp_cumulative_start = cumsum(c(0, chr_size[-dplyr::n()] + 1)) +
        (seq_len(dplyr::n()) - 1),
      bp_cumulative_end = bp_cumulative_start + chr_size
    )

  gwas <- original_gwas |>
    dplyr::filter(!is.na(LP)) |>
    dplyr::mutate(bin = floor(BP / 100000)) |>
    dplyr::group_by(CHR, bin) |>
    dplyr::summarise(
      BP = mean(BP),
      LP = max(LP),
      .groups = "drop"
    ) |>
    dplyr::arrange(CHR, BP) |>
    dplyr::left_join(
      chr_ranges |> dplyr::select(CHR, bp_cumulative_start, bp_min), by = "CHR"
    ) |>
    dplyr::mutate(BP_cumulative = BP - bp_min + bp_cumulative_start) |>
    # Add a flag for non-adjacent points
    dplyr::mutate(
      is_adjacent = CHR == dplyr::lead(CHR) &
        abs(BP - dplyr::lead(BP)) < 100000 * 2 # Allow for 2x bin size gap
    )

  # Create a second dataset with larger bins for the area
  gwas_area <- gwas |>
    dplyr::mutate(bin = floor(BP / 1000000)) |> # Create 1Mb bins for the area
    dplyr::group_by(CHR, bin) |>
    dplyr::summarise(
      BP = mean(BP),
      LP = max(LP),
      .groups = "drop"
    ) |>
    dplyr::arrange(CHR, BP) |>
    dplyr::left_join(
      chr_ranges |> dplyr::select(CHR, bp_cumulative_start, bp_min), by = "CHR"
    ) |>
    dplyr::mutate(BP_cumulative = BP - bp_min + bp_cumulative_start)

  # Add start and end rows with Z = 0 for each chromosome
  gwas_area <- gwas_area |>
    dplyr::group_by(CHR) |>
    dplyr::arrange(BP_cumulative, .by_group = TRUE) |>
    dplyr::group_modify(~ {
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
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max_lp_y_axis)) +
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
    width = plot_width / 72, # Convert to inches for ggsave
    height = plot_height / 72, # Convert to inches for ggsave
    units = "in", # Specify inches
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
      ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, max_lp_y_axis)) +
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

    svg_location <- glue::glue("{study$extracted_location}/svgs/chr{chr_num}.svg")
    ggplot2::ggsave(
      svg_location,
      p_chr,
      width = plot_width / 72,
      height = plot_height / 72,
      units = "in",
      limitsize = FALSE
    )
  }
  old_wd <- getwd()
  setwd(glue::glue("{study$extracted_location}/svgs"))
  files_to_zip <- Sys.glob("*.svg")
  utils::zip("full.zip", files_to_zip)
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

  jsonlite::write_json(
    metadata,
    glue::glue("{study$extracted_location}/svgs/full.json"),
    pretty = TRUE,
    auto_unbox = TRUE
  )
  return()
}

prepare_svg_files_for_use <- function(studies_db_file, do_all = FALSE) {
  original_wd <- getwd()
  dir.create(glue::glue("{svg_dir}/traits"), showWarnings = F, recursive = T)
  dir.create(glue::glue("{svg_dir}/groups"), showWarnings = F, recursive = T)
  dir.create(glue::glue("{svg_dir}/extractions"), showWarnings = F, recursive = T)

  # find out which new studies and extractions are missing from the svg directory
  current_studies_conn <- duckdb::dbConnect(duckdb::duckdb(), studies_db_file, read_only = TRUE)

  studies_query <- "SELECT * FROM studies WHERE data_type = 'phenotype' AND variant_type = 'common'"
  current_studies <- DBI::dbGetQuery(current_studies_conn, studies_query)
  existing_full_svgs <- list.files(glue::glue("{svg_dir}/traits"), full.names = TRUE)
  current_studies$new_full_svgs <- glue::glue("{svg_dir}/traits/{current_studies$id}_svgs.zip")
  new_studies <- current_studies[!current_studies$new_full_svgs %in% existing_full_svgs, ]

  if (do_all) {
    new_studies <- current_studies
  }
  message(glue::glue("Preparing {nrow(new_studies)} full trait svgs"))

  lapply(seq_len(nrow(new_studies)), function(i) {
    study <- new_studies[i, , drop = FALSE]
    file.link(
      glue::glue("{extracted_study_dir}{study$study_name}/svgs/full.zip"),
      glue::glue("{svg_dir}/traits/{study$trait_id}_svgs.zip")
    )
    file.link(
      glue::glue("{extracted_study_dir}{study$study_name}/svgs/full.json"),
      glue::glue("{svg_dir}/traits/{study$trait_id}_metadata.json")
    )
    return()
  })

  study_extractions_query <- "SELECT study_extractions.id, study_extractions.svg_file
    FROM study_extractions LEFT JOIN studies ON study_extractions.study_id = studies.id
    WHERE studies.variant_type = 'common'
  "
  current_study_extractions <- DBI::dbGetQuery(current_studies_conn, study_extractions_query)
  current_study_extractions$new_svg_file <- glue::glue("{svg_dir}/extractions/{current_study_extractions$id}.svg")

  existing_extraction_svgs <- list.files(glue::glue("{svg_dir}/extractions"), full.names = TRUE)
  missing_extraction_svgs <- current_study_extractions[
    !current_study_extractions$new_svg_file %in% existing_extraction_svgs,
  ]

  if (do_all) {
    missing_extraction_svgs <- current_study_extractions
  }
  message(glue::glue("Preparing {nrow(missing_extraction_svgs)} extraction svgs"))

  copied_extraction_svgs <- lapply(seq_len(nrow(missing_extraction_svgs)), function(i) {
    extraction <- missing_extraction_svgs[i, , drop = FALSE]
    file.link(glue::glue("{data_dir}{extraction$svg_file}"), extraction$new_svg_file)
    return()
  })

  coloc_groups_query <- "SELECT * FROM coloc_groups"
  coloc_groups <- DBI::dbGetQuery(current_studies_conn, coloc_groups_query)
  unique_coloc_group_ids <- unique(coloc_groups$coloc_group_id)
  message(glue::glue("Preparing {length(unique_coloc_group_ids)} coloc group svgs"))

  setwd(glue::glue("{svg_dir}/groups"))
  file.remove(Sys.glob("*.zip"))

  # zip the svg extractions by coloc group
  completed <- parallel::mclapply(unique_coloc_group_ids, mc.cores = 10, function(group_id) {
    tryCatch({
      zip_file <- glue::glue("coloc_group_{group_id}_svgs.zip")
      # if (file.exists(zip_file)) return()

      specific_coloc_group <- coloc_groups |>
        dplyr::filter(coloc_group_id == group_id)

      extraction_ids_in_group <- specific_coloc_group$study_extraction_id
      svg_files <- paste0(svg_dir, "/extractions/", extraction_ids_in_group, ".svg")

      zip::zipr(glue::glue("coloc_group_{group_id}_svgs.zip"), svg_files)
      return()
    },
    error = function(e) {
      message(glue::glue('Error zipping group {group_id}: {e$message}, {paste(svg_files, collapse = ", ")}'))
      return()
    })
    return()
  })
  print(warnings())
  setwd(original_wd)
  return()
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
