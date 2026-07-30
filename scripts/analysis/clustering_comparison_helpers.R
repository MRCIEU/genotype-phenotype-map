analyze_graph_structure_advanced <- function(g, original_graph = NULL) {
  if (is.null(g) || igraph::vcount(g) < 1L) {
    return(list(
      n_components = 0L,
      mean_weight = NA_real_,
      mean_transitivity = NA_real_,
      mean_conductance = NA_real_,
      median_size = NA_real_
    ))
  }

  graph_components <- igraph::components(g)
  if (graph_components$no < 1L) {
    return(list(
      n_components = 0L,
      mean_weight = NA_real_,
      mean_transitivity = NA_real_,
      mean_conductance = NA_real_,
      median_size = NA_real_
    ))
  }

  component_metrics <- lapply(seq_len(graph_components$no), function(component_id) {
    vertex_indices <- which(graph_components$membership == component_id)
    component_graph <- igraph::induced_subgraph(g, vertex_indices)

    mean_weight <- if (igraph::ecount(component_graph) > 0L) {
      mean(igraph::E(component_graph)$weight)
    } else {
      0
    }
    component_transitivity <- igraph::transitivity(component_graph, type = "global")
    if (is.na(component_transitivity)) component_transitivity <- 0

    component_conductance <- NA_real_
    if (!is.null(original_graph)) {
      vertex_names <- igraph::V(g)[vertex_indices]$name
      internal_weight <- sum(igraph::E(component_graph)$weight)
      incident_edges <- igraph::incident_edges(original_graph, v = vertex_names, mode = "all")
      all_edge_ids <- unique(unlist(incident_edges))
      total_original_weight <- sum(igraph::E(original_graph)[all_edge_ids]$weight)
      component_conductance <- if (total_original_weight > 0) {
        (total_original_weight - internal_weight) / total_original_weight
      } else {
        0
      }
    }

    return(data.frame(
      weight = mean_weight,
      transitivity = component_transitivity,
      conductance = component_conductance,
      size = igraph::vcount(component_graph)
    ))
  })

  component_metrics <- dplyr::bind_rows(component_metrics)
  return(list(
    n_components = graph_components$no,
    mean_weight = mean(component_metrics$weight, na.rm = TRUE),
    mean_transitivity = mean(component_metrics$transitivity, na.rm = TRUE),
    mean_conductance = mean(component_metrics$conductance, na.rm = TRUE),
    median_size = stats::median(component_metrics$size)
  ))
}

undirected_pair_id <- function(uid_a, uid_b) {
  return(
    ifelse(uid_a < uid_b, paste(uid_a, uid_b, sep = "||"), paste(uid_b, uid_a, sep = "||"))
  )
}

build_truthset_ti_meta <- function(truthset) {
  return(
    truthset |>
      dplyr::group_by(ti_uid) |>
      dplyr::summarise(
        pairwisecoloc_support = any(dplyr::coalesce(as.logical(pairwisecoloc_support), FALSE)),
        production_gpmap_support = any(dplyr::coalesce(as.logical(gpmap_support), FALSE)),
        production_preinfomap_support = any(dplyr::coalesce(as.logical(preinfomap_support), FALSE)),
        .groups = "drop"
      ) |>
      dplyr::mutate(
        support_tier = dplyr::case_when(
          production_gpmap_support & production_preinfomap_support ~ "pairwise_preinfomap_and_gpmap",
          !production_gpmap_support & production_preinfomap_support ~ "pairwise_and_preinfomap_only",
          production_gpmap_support & !production_preinfomap_support ~ "pairwise_and_gpmap_only",
          TRUE ~ "pairwise_only"
        ),
        gap_target = pairwisecoloc_support & !production_gpmap_support
      )
  )
}

evaluate_mode_ti_recovery <- function(input_row, truthset_uid_mappings) {
  clusters <- vroom::vroom(
    input_row$cluster_file,
    delim = "\t",
    col_select = dplyr::any_of(c("unique_study_id", "component")),
    show_col_types = FALSE
  )
  if (!all(c("unique_study_id", "component") %in% names(clusters))) {
    stop("Missing required cluster columns in ", input_row$cluster_file)
  }

  component_lookup <- stats::setNames(clusters$component, clusters$unique_study_id)
  block_mappings <- truthset_uid_mappings |>
    dplyr::filter(.data$ld_block == input_row$ld_block)

  pair_recovery <- block_mappings |>
    dplyr::mutate(
      molecular_component = component_lookup[molecular_unique_study_id],
      phenotype_component = component_lookup[phenotype_unique_study_id],
      cluster_support = !is.na(molecular_component) &
        !is.na(phenotype_component) &
        molecular_component == phenotype_component
    ) |>
    dplyr::group_by(ti_uid) |>
    dplyr::summarise(cluster_support = any(cluster_support), .groups = "drop")

  return(
    block_mappings |>
      dplyr::distinct(ti_uid) |>
      dplyr::left_join(pair_recovery, by = "ti_uid") |>
      dplyr::mutate(
        cluster_support = dplyr::coalesce(cluster_support, FALSE),
        ld_block = input_row$ld_block,
        mode = input_row$mode,
        mode_label = input_row$mode_label
      )
  )
}

evaluate_mode_coloc_burden <- function(
    input_row,
    coloc_file,
    h4_threshold,
    study_types,
    truthset_uid_mappings) {
  if (!file.exists(coloc_file)) {
    return(data.frame(
      ld_block = input_row$ld_block,
      mode = input_row$mode,
      mode_label = input_row$mode_label,
      n_h4_mol_pheno_pairs = 0L,
      n_same_cluster_h4_mol_pheno = 0L,
      n_truthset_supported_same_cluster = 0L,
      n_non_truthset_same_cluster = 0L,
      n_pipeline_false_positive = 0L,
      n_pipeline_false_negative = 0L
    ))
  }

  clusters <- vroom::vroom(
    input_row$cluster_file,
    delim = "\t",
    col_select = dplyr::any_of(c("unique_study_id", "component")),
    show_col_types = FALSE
  )
  coloc_results <- vroom::vroom(
    coloc_file,
    delim = "\t",
    show_col_types = FALSE
  )
  if (nrow(coloc_results) == 0L || nrow(clusters) == 0L) {
    return(data.frame(
      ld_block = input_row$ld_block,
      mode = input_row$mode,
      mode_label = input_row$mode_label,
      n_h4_mol_pheno_pairs = 0L,
      n_same_cluster_h4_mol_pheno = 0L,
      n_truthset_supported_same_cluster = 0L,
      n_non_truthset_same_cluster = 0L,
      n_pipeline_false_positive = 0L,
      n_pipeline_false_negative = 0L
    ))
  }

  study_to_component <- stats::setNames(clusters$component, clusters$unique_study_id)
  coloc_all <- coloc_results |>
    dplyr::filter(dplyr::coalesce(ignore, FALSE) == FALSE) |>
    dplyr::left_join(study_types, by = c("study_a" = "study_name")) |>
    dplyr::rename(type_a = data_type) |>
    dplyr::left_join(study_types, by = c("study_b" = "study_name")) |>
    dplyr::rename(type_b = data_type) |>
    dplyr::mutate(
      component_a = study_to_component[unique_study_a],
      component_b = study_to_component[unique_study_b],
      same_cluster = !is.na(component_a) & !is.na(component_b) & component_a == component_b,
      pair_id = undirected_pair_id(unique_study_a, unique_study_b)
    )

  coloc_h4 <- coloc_all |>
    dplyr::filter(!is.na(h4), h4 >= h4_threshold) |>
    dplyr::filter(
      (type_a == "phenotype" & type_b != "phenotype") |
        (type_b == "phenotype" & type_a != "phenotype")
    )

  truthset_pair_ids <- truthset_uid_mappings |>
    dplyr::filter(.data$ld_block == input_row$ld_block) |>
    dplyr::mutate(
      pair_id = undirected_pair_id(molecular_unique_study_id, phenotype_unique_study_id)
    ) |>
    dplyr::pull(pair_id) |>
    unique()

  same_cluster_pairs <- coloc_h4 |>
    dplyr::filter(same_cluster)

  clustered_result <- readRDS(input_row$rds_file)
  pruned_edges <- clustered_result$pruned_edges
  n_pipeline_false_positive <- 0L
  if (!is.null(pruned_edges) && nrow(pruned_edges) > 0L) {
    pruned_pairs <- data.frame(
      study_a = c(pruned_edges$V1, pruned_edges$V2),
      study_b = c(pruned_edges$V2, pruned_edges$V1)
    )
    matches <- dplyr::inner_join(
      coloc_h4,
      pruned_pairs,
      by = c("unique_study_a" = "study_a", "unique_study_b" = "study_b")
    )
    n_pipeline_false_positive <- nrow(matches)
  }

  n_pipeline_false_negative <- sum(
    coloc_all$same_cluster & !is.na(coloc_all$h4) & coloc_all$h4 < h4_threshold,
    na.rm = TRUE
  )

  return(data.frame(
    ld_block = input_row$ld_block,
    mode = input_row$mode,
    mode_label = input_row$mode_label,
    n_h4_mol_pheno_pairs = nrow(coloc_h4),
    n_same_cluster_h4_mol_pheno = nrow(same_cluster_pairs),
    n_truthset_supported_same_cluster = sum(same_cluster_pairs$pair_id %in% truthset_pair_ids),
    n_non_truthset_same_cluster = sum(!same_cluster_pairs$pair_id %in% truthset_pair_ids),
    n_pipeline_false_positive = n_pipeline_false_positive,
    n_pipeline_false_negative = n_pipeline_false_negative
  ))
}

summarise_clustered_graphs <- function(clustered_results) {
  if (
    is.null(clustered_results) ||
      is.null(clustered_results$unpruned_graph) ||
      is.null(clustered_results$pruned_graph)
  ) {
    stop("Clustered result must contain unpruned_graph and pruned_graph")
  }

  unpruned <- analyze_graph_structure_advanced(clustered_results$unpruned_graph)
  pruned <- analyze_graph_structure_advanced(
    clustered_results$pruned_graph,
    original_graph = clustered_results$unpruned_graph
  )

  return(data.frame(
    n_components_unpruned = unpruned$n_components,
    n_components_pruned = pruned$n_components,
    mean_weight_unpruned = unpruned$mean_weight,
    mean_weight_pruned = pruned$mean_weight,
    mean_transitivity_unpruned = unpruned$mean_transitivity,
    mean_transitivity_pruned = pruned$mean_transitivity,
    mean_conductance_pruned = pruned$mean_conductance,
    median_size_unpruned = unpruned$median_size,
    median_size_pruned = pruned$median_size
  ))
}
