source('constants.R')
library(dplyr)
dplyr.summarise.inform = FALSE
bp_range <- 50000
subgraph_density_theshold <- 0.6
min_internal_degree_percentage <- 0.05

parser <- argparser::arg_parser('Colocalise studies per region')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block that the studies are in', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Coloc result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--worker_guid', help = 'Worker GUID', type = 'character', default = NA)
parser <- argparser::add_argument(parser, '--force_clustering', help = 'Force clustering', type = 'logical', default = FALSE)

args <- argparser::parse_args(parser)

main <- function() {
  start_time <- Sys.time()
  if (!is.na(args$worker_guid)) {
    update_directories_for_worker(args$worker_guid)
  }

  ld_info <- ld_block_dirs(args$ld_block)
  coloc_results_file <- glue::glue('{ld_info$ld_block_data}/coloc_pairwise_results.tsv.gz')
  if (file.exists(coloc_results_file)) {
    coloc_results <- vroom::vroom(coloc_results_file, delim = '\t', show_col_types = F)
  } else {
    coloc_results <- data.frame()
  }

  block <- vroom::vroom(glue::glue('{pipeline_metadata_dir}updated_ld_blocks_to_colocalise.tsv'), show_col_types=F) |>
    dplyr::filter(data_dir == ld_info$ld_block_data)

  finemapped_file <- glue::glue('{ld_info$ld_block_data}/finemapped_studies.tsv')
  if (file.exists(finemapped_file)) {
    finemapped_studies <- vroom::vroom(finemapped_file, col_types = finemapped_column_types, show_col_types = F) |>
      dplyr::filter(variant_type == variant_types$common & min_p <= p_value_threshold & !ignore) |>
      dplyr::arrange(unique_study_id)
  }

  if (!is.na(args$worker_guid)) {
    existing_finemapped_studies_file <- glue::glue('{data_dir}/ld_blocks/{args$ld_block}/finemapped_studies.tsv')
    existing_finemapped_studies <- vroom::vroom(existing_finemapped_studies_file, col_types = finemapped_column_types, show_col_types=F)
    finemapped_studies <- dplyr::bind_rows(finemapped_studies, existing_finemapped_studies) |>
      dplyr::filter(variant_type == variant_types$common & min_p <= p_value_threshold & !ignore) |>
      dplyr::mutate(file = sub('/local-scratch/projects/genotype-phenotype-map/data/', data_dir, file)) |>
      dplyr::arrange(unique_study_id)
  }

  if (!file.exists(finemapped_file) || nrow(block) == 0 || nrow(finemapped_studies) == 0) {
    message(glue::glue('{args$ld_block}: Nothing to coloc, skipping.'))
    vroom::vroom_write(data.frame(), args$completed_output_file)
    return()
  }

  finemapped_studies <- dplyr::arrange(finemapped_studies, unique_study_id)

  # Get all possible pairs within bp_range, excluding already calculated pairs
  study_pairs <- get_study_pairs_to_coloc(finemapped_studies, coloc_results, args$worker_guid)
  message(glue::glue('{args$ld_block}: Found {nrow(study_pairs)} study pairs to coloc in {diff_time_taken(start_time)}'))

  if (nrow(study_pairs) == 0) {
    message(glue::glue('{args$ld_block}: No study pairs to coloc, skipping.'))
    if (args$force_clustering == FALSE) {
      vroom::vroom_write(data.frame(), args$completed_output_file)
      return()
    }
  }

  if (!is.na(args$worker_guid)) {
    finemapped_subset <- dplyr::filter(finemapped_studies,
      unique_study_id %in% study_pairs$unique_study_a | unique_study_id %in% study_pairs$unique_study_b
    )
  } else {
    finemapped_subset <- finemapped_studies
  }

  studies_to_colocalise <- lapply(finemapped_subset$file, function(file) {
    if (file.info(file)$size == 0) {
      message(glue::glue('{file} is empty, delete.'))
      return(NULL)
    }

    gwas <- vroom::vroom(file, delim = '\t',
      show_col_types = F,
      col_select = c('SNP', 'LBF'),
      altrep = F
    )
    return(gwas)
  })
  names(studies_to_colocalise) <- finemapped_subset$unique_study_id
  message(glue::glue('{args$ld_block}: Loaded {length(studies_to_colocalise)} studies in {diff_time_taken(start_time)}'))

  new_coloc_results <- lapply(seq_len(nrow(study_pairs)), function(i) {
    pair <- study_pairs[i,]
    first_gwas <- studies_to_colocalise[[pair$unique_study_a]]
    second_gwas <- studies_to_colocalise[[pair$unique_study_b]]

    tryCatch({
      result <- pairwise_coloc_analysis(first_gwas, second_gwas)
    }, error = function(e) {
      message(glue::glue('Error colocating {pair$unique_study_a} and {pair$unique_study_b}: {e}'))
      stop(glue::glue('Error colocating {pair$unique_study_a} and {pair$unique_study_b}: {e}'))
    })
    if (is.null(result)) {
      result <- data.frame(
        ld_block = args$ld_block,
        unique_study_a = pair$unique_study_a,
        study_a = NA,
        unique_study_b = pair$unique_study_b,
        study_b = NA,
        ignore = T,
        spurious = F,
        bp_distance = NA,
        nsnps = NA,
        hit1 = NA,
        hit2 = NA,
        idx1 = NA,
        idx2 = NA,
        PP.H0.abf = NA,
        PP.H1.abf = NA,
        PP.H2.abf = NA,
        PP.H3.abf = NA,
        PP.H4.abf = NA,
        h4 = NA
      )
      return(result)
    }

    result <- dplyr::bind_cols(pair, result)
    return(result)
  })

  new_coloc_results <- dplyr::bind_rows(new_coloc_results[!sapply(new_coloc_results, is.null)])
  message(glue::glue('{args$ld_block}: Colocated {nrow(new_coloc_results)} study pairs in {diff_time_taken(start_time)}'))

  if ((!is.null(nrow(new_coloc_results)) && nrow(new_coloc_results) > 0) || args$force_clustering == TRUE) {
    coloc_results <- dplyr::bind_rows(coloc_results, new_coloc_results) |>
      dplyr::distinct(unique_study_a, unique_study_b, .keep_all = TRUE) |>
      dplyr::mutate(ld_block = args$ld_block)

    results <- prune_poor_finemapping_results(finemapped_studies, coloc_results)
    finemapped_studies <- results$finemapped_studies
    coloc_results <- results$coloc_results
    vroom::vroom_write(finemapped_studies, finemapped_file)

    coloc_results_with_min_p <- coloc_results |>
      dplyr::mutate(
        min_p_study_a = finemapped_studies[match(unique_study_a, finemapped_studies$unique_study_id),]$min_p,
        min_p_study_b = finemapped_studies[match(unique_study_b, finemapped_studies$unique_study_id),]$min_p
      )

    clustered_results <- cluster_coloc_results(coloc_results_with_min_p, start_time)

    coloc_results <- mark_spurious_colocs(coloc_results, clustered_results$pruned_edges)
    
    message(glue::glue('{args$ld_block}: Marked {sum(coloc_results$spurious)} spurious colocs'))
    vroom::vroom_write(coloc_results, coloc_results_file)

    additional_data_per_cluster <- find_snp_and_connectedness_per_cluster(clustered_results$groups, studies_to_colocalise, coloc_results_with_min_p)

    clustered_results$groups <- clustered_results$groups |>
      dplyr::mutate(ld_block = args$ld_block) |>
      dplyr::left_join(additional_data_per_cluster, by = "component") |>
      dplyr::arrange(component)

    if (nrow(clustered_results$groups) == 0) {
      clustered_results$groups <- data.frame(
        unique_study_id = character(),
        component = integer(),
        ld_block = character(),
        snp = character(),
        h4_connectedness = numeric(),
        h3_connectedness = numeric()
      )
    }

    clustered_results$groups <- dplyr::bind_rows(clustered_results$groups)
    saveRDS(clustered_results, glue::glue('{ld_info$ld_block_data}/igraph_clustered_results.rds'))
    vroom::vroom_write(coloc_results, coloc_results_file)
    vroom::vroom_write(clustered_results$groups, glue::glue('{ld_info$ld_block_data}/coloc_clustered_results.tsv.gz'))
  }

  vroom::vroom_write(data.frame(), args$completed_output_file)
}

#' get_study_pairs_to_coloc takes a list of studies and a list of existing coloc results,
#' and returns a list of study pairs to colocalise
#' @param studies: list of finemapped studies to colocalise
#' @param existing_results: list of existing coloc results
#' @param worker_guid: worker GUID
#' @returns list of study pairs to colocalise, with the bp distance between the studies
#' @import data.table
#' @export
get_study_pairs_to_coloc <- function(studies, existing_results, worker_guid) {
  studies <- dplyr::mutate(studies, id = dplyr::row_number())
  studies <- data.table::as.data.table(studies)

  pairs_filtered <- studies[ studies, on = .(id < id), allow.cartesian = TRUE ][
    , bp_distance := abs(i.bp - bp) ][
    bp_distance <= bp_range | i.study == study ][
    , .(
      unique_study_a = unique_study_id,
      study_a = study,
      unique_study_b = i.unique_study_id,
      study_b = i.study,
      bp_distance = bp_distance,
      ignore = F,
      spurious = F
    )
  ] |>
    tibble::as_tibble()

  if (!is.na(worker_guid)) {
    pairs_filtered <- pairs_filtered |>
      dplyr::filter(study_a == worker_guid | study_b == worker_guid)
  }

  existing_results <- data.table::as.data.table(existing_results)
  if (nrow(existing_results) > 0) {
    pairs_filtered <- dplyr::anti_join(pairs_filtered, existing_results, by = c("unique_study_a", "unique_study_b")) 
  }

  return(pairs_filtered)
}

#' run_coloc_analysis takes two already harmonised gwases, and runs coloc on the results
#' @param first_gwas: first gwas to be run through coloc.  This is the gwas that results will be based off
#' @param second_gwas: second gwas to be run through coloc
#' @returns tibble of coloc results (h0 - h4)
#' @import coloc
#' @import tibble
#' @export
pairwise_coloc_analysis <- function(first_gwas, second_gwas) {
  harmonised_gwases <- harmonise_gwases(first_gwas, second_gwas)
  if (length(harmonised_gwases) == 0) return(NULL)
  
  first_gwas <- harmonised_gwases[[1]]
  second_gwas <- harmonised_gwases[[2]]

  if (nrow(first_gwas) < 50 || nrow(second_gwas) < 50) return(NULL)

  first_lbf <- first_gwas$LBF
  names(first_lbf) <- first_gwas$SNP
  second_lbf <- second_gwas$LBF
  names(second_lbf) <- second_gwas$SNP

  result <- coloc::coloc.bf_bf(bf1 = first_lbf, bf2 = second_lbf)
  result$summary$h4 <- result$summary$PP.H4.abf
  coloc_results <- result$summary
  return(coloc_results)
}

harmonise_gwases <- function(...) {
  gwases <- list(...)
  snpids <- Reduce(intersect, lapply(gwases, function(gwas) as.character(gwas$SNP)))
  if (length(snpids) <= 1) return(list())
  snpids <- sort(snpids)

  gwases <- lapply(gwases, function(gwas) {
    gwas <- as.data.frame(gwas)
    gwas <- gwas[!duplicated(gwas$SNP), ]
    gwas <- gwas[gwas$SNP %in% snpids, ]
    gwas <- gwas[match(snpids, gwas$SNP), ]
    return(gwas)
  })

  stopifnot(identical(gwases[[1]]$SNP, gwases[[2]]$SNP))
  return(gwases)
}

# Check for within study colocalising finemapped regions (requires pairwise coloc to between credible sets)
prune_poor_finemapping_results <- function(finemapped_studies, coloc_results) {
  bad_coloc_results <- coloc_results |>
    dplyr::filter(study_a == study_b & h4 >= posterior_prob_threshold) |>
      dplyr::mutate(
        min_p_study_a = finemapped_studies[match(unique_study_a, finemapped_studies$unique_study_id),]$min_p,
        min_p_study_b = finemapped_studies[match(unique_study_b, finemapped_studies$unique_study_id),]$min_p,
        bad_unique_study = ifelse(min_p_study_a > min_p_study_b, unique_study_a, unique_study_b)
      ) |>
      dplyr::select(unique_study_a, unique_study_b, study_a, study_b, h4, min_p_study_a, min_p_study_b, bad_unique_study)
  
  if (nrow(bad_coloc_results) == 0) {
    return(list(finemapped_studies = finemapped_studies, coloc_results = coloc_results))
  }

  poor_finemapped_studies <- unique(bad_coloc_results$bad_unique_study)

  finemapped_studies <- finemapped_studies |>
    dplyr::mutate(ignore = ignore | unique_study_id %in% poor_finemapped_studies)
  
  coloc_results <- coloc_results |> 
    dplyr::mutate(ignore = ignore | unique_study_a %in% poor_finemapped_studies | unique_study_b %in% poor_finemapped_studies)

  message(glue::glue('{args$ld_block}: New number of studies to ignore: {length(poor_finemapped_studies)}'))
  return(list(finemapped_studies = finemapped_studies, coloc_results = coloc_results))
}

cluster_coloc_results <- function(coloc_results, start_time) {
  coloc_results <- coloc_results |> dplyr::filter(!ignore)
  message(glue::glue('{args$ld_block}: Clustering {nrow(coloc_results)} coloc results starting {diff_time_taken(start_time)}'))

  h4_adj_mx <- make_adjacency_matrix(coloc_results)

  h4_graph <- igraph::graph_from_adjacency_matrix(h4_adj_mx, mode="undirected", weighted=TRUE, diag=FALSE)

  # Initial pruning of singleton studies (vertices) with no connections from the overall graph
  graph_components <- igraph::components(h4_graph)
  vert_out_initial <- igraph::V(h4_graph)[graph_components$membership %in% which(graph_components$csize == 1)]
  h4_graph <- igraph::delete_vertices(h4_graph, vert_out_initial)

  message(glue::glue('{args$ld_block}: Pruning graph starting {diff_time_taken(start_time)}'))

  h4_graph_components <- igraph::components(h4_graph)
  recombined_graph <- h4_graph
  subgraph_to_cluster_and_prune <- NULL
  modified_pruned_graph <- NULL
  pruned_studies <- c()
  pruned_edges <- data.frame(V1 = character(), V2 = character())

  if(h4_graph_components$no == 1) {
    message(glue::glue('Only one component found, assigning all vertices to component 1'))
    igraph::V(h4_graph)$component <- 1
  } else {
    message(glue::glue('{args$ld_block}: Clustering {h4_graph_components$no} components'))
    subgraphs <- lapply(1:h4_graph_components$no, function(x){
      igraph::induced_subgraph(h4_graph, vids = igraph::V(h4_graph)[h4_graph_components$membership == x])
    })
    subgraph_density <- sapply(subgraphs, function(x){igraph::edge_density(x)})

    if(all(subgraph_density >= subgraph_density_theshold)){
      # If all components are dense enough, just assign their component IDs
      igraph::V(h4_graph)$component <- h4_graph_components$membership
      recombined_graph <- h4_graph
    } else {
      subgraph_to_keep <- igraph::disjoint_union(subgraphs[subgraph_density >= subgraph_density_theshold])
      subgraph_to_cluster_and_prune <- igraph::disjoint_union(subgraphs[subgraph_density < subgraph_density_theshold])
      modified_pruned_graph <- subgraph_to_cluster_and_prune

      pruned_edge_ids <- c()
      
      if (igraph::vcount(modified_pruned_graph) > 0) {
        clustered_graph <- igraph::cluster_infomap(modified_pruned_graph)
        clustered_memberships <- igraph::membership(clustered_graph)


        # FIrst, remove all edges between communities
        edge_list_names <- igraph::as_edgelist(modified_pruned_graph)
        if (nrow(edge_list_names) > 0) {
          for (i in 1:nrow(edge_list_names)) {
            v1_name <- edge_list_names[i, 1]
            v2_name <- edge_list_names[i, 2]
            comm_v1 <- clustered_memberships[v1_name]
            comm_v2 <- clustered_memberships[v2_name]

            if (comm_v1 != comm_v2) {
              pruned_edges <- rbind(pruned_edges, 
                data.frame(V1 = edge_list_names[i, 1], V2 = edge_list_names[i, 2]))
              pruned_edge_ids <- c(pruned_edge_ids, i)
            }
          }
        }

        # Second, remove orphaned vertices after edge pruning
        current_degrees <- igraph::degree(modified_pruned_graph)
        pruned_studies <- c(pruned_studies, names(current_degrees[current_degrees == 0]))

        # Then, remove all studies with degree less than the threshold
        for (comm_id in unique(clustered_memberships)) {
          current_community_vertices <- names(clustered_memberships[clustered_memberships == comm_id])

          if (length(current_community_vertices) == 0) {
            next
          }

          community_graph <- igraph::induced_subgraph(subgraph_to_cluster_and_prune, vids = current_community_vertices)
          internal_degrees <- igraph::degree(community_graph)
          max_degree_in_community <- max(internal_degrees)
          threshold_degree <- max_degree_in_community * min_internal_degree_percentage

          pruned_studies <- c(pruned_studies, names(internal_degrees[internal_degrees < threshold_degree]))
        }
        pruned_studies <- unique(pruned_studies)

        message(glue::glue('Removing {nrow(pruned_edges)} edges'))
        if (nrow(pruned_edges) > 0) {
          modified_pruned_graph <- igraph::delete_edges(
            modified_pruned_graph,
            edges = pruned_edge_ids
          )
        }

        message(glue::glue('Removing {length(pruned_studies)} vertices'))
        if (length(pruned_studies) > 0) {
          modified_pruned_graph <- igraph::delete_vertices(modified_pruned_graph, pruned_studies)
        }
      }

      # Adding edges from pruned studies, to mark as spurious later
      for (study in pruned_studies) {
        study_edges <- igraph::incident_edges(subgraph_to_cluster_and_prune, study)
        if (length(study_edges) > 0) {
          for (edge in study_edges) {
            edge_vertices <- igraph::ends(subgraph_to_cluster_and_prune, edge)
            pruned_edges <- rbind(pruned_edges, 
              data.frame(V1 = edge_vertices[1], V2 = edge_vertices[2]))
          }
        }
      }
      pruned_edges <- unique(pruned_edges)

      # Recombine the 'kept' components and the 'modified pruned' components
      recombined_graph <- igraph::disjoint_union(subgraph_to_keep, modified_pruned_graph)
      recombined_graph_components <- igraph::components(recombined_graph)
      igraph::V(recombined_graph)$component <- recombined_graph_components$membership
    }
  }

  message(glue::glue('Outputting results in {diff_time_taken(start_time)}'))
  coloc_groups <- data.frame(unique_study_id = igraph::vertex_attr(recombined_graph)$name, component = igraph::V(recombined_graph)$component)

  message(glue::glue('Summary of pruning:'))
  message(glue::glue('  Initial vertices (after H4 adj matrix): {igraph::vcount(h4_graph)}'))
  message(glue::glue('  Initial edges (after H4 adj matrix): {igraph::ecount(h4_graph)}'))
  message(glue::glue('  Final vertices after all pruning: {igraph::vcount(recombined_graph)}'))
  message(glue::glue('  Final edges after all pruning: {igraph::ecount(recombined_graph)}'))
  message(glue::glue('  Total vertices removed from clustering: {length(pruned_studies)}'))
  message(glue::glue('  Total edges removed from clustering: {nrow(pruned_edges)}'))

  return(list(
    unpruned_graph = h4_graph,
    pruned_graph = recombined_graph,
    groups = coloc_groups,
    pruned_edges = pruned_edges
  ))
}

mark_spurious_colocs <- function(coloc_results, pruned_edges) {
  coloc_results$spurious <- FALSE
  
  if (nrow(pruned_edges) == 0) {
    return(coloc_results)
  coloc_results$unique_study_a %in% group_studies}
  
  pruned_pairs <- data.frame(
    study_a = c(pruned_edges$V1, pruned_edges$V2),
    study_b = c(pruned_edges$V2, pruned_edges$V1)
  )
  
  matches <- merge(coloc_results, pruned_pairs, 
                   by.x = c("unique_study_a", "unique_study_b"),
                   by.y = c("study_a", "study_b"),
                   all.x = FALSE)
  
  if (nrow(matches) > 0) {
    coloc_results$pair_id <- paste(coloc_results$unique_study_a, coloc_results$unique_study_b, sep = "_")
    matches$pair_id <- paste(matches$unique_study_a, matches$unique_study_b, sep = "_")
    
    coloc_results$spurious[coloc_results$pair_id %in% matches$pair_id] <- TRUE
    coloc_results$pair_id <- NULL
  }
  
  return(coloc_results)
}

#' make_adjacency_matrix takes a dataframe of coloc results and returns a symmetrical adjacency matrix
#' @param coloc_results: dataframe of coloc results
#' @returns symmetrical adjacency matrix
make_adjacency_matrix <- function(coloc_results) {
  pvals <- pmax(coloc_results$min_p_study_a, coloc_results$min_p_study_b)
  weights <- weights_from_pvals(pvals) * ifelse(coloc_results$h4 >= posterior_prob_h4_threshold, 1, 0)
  coloc_results <- coloc_results |> dplyr::mutate(weight = weights)

  h4_adj_mx <- rbind(
    coloc_results |> 
      dplyr::select(a = unique_study_a, b = unique_study_b, weight = weight) |>
      dplyr::group_by(a, b) |>
      dplyr::summarise(weight = max(weight), .groups = 'keep'),
    coloc_results |> 
      dplyr::select(a = unique_study_b, b = unique_study_a, weight = weight) |>
      dplyr::group_by(a, b) |>
      dplyr::summarise(weight = max(weight), .groups = 'keep')
  ) |>
    dplyr::distinct() |>
    tidyr::pivot_wider(names_from = a, values_from = weight) |>
    tibble::column_to_rownames("b") |>
    as.matrix()
  
  h4_adj_mx <- h4_adj_mx[,rownames(h4_adj_mx)]
  diag(h4_adj_mx) <- 0
  h4_adj_mx[is.na(h4_adj_mx)] <- 0
  h4_adj_mx[h4_adj_mx < posterior_prob_h4_threshold] <- 0

  return(h4_adj_mx)
}


# Generate weights from p-values
# More significant than threshold2 gets weight 1
# Less significant than threshold1 gets weight 0
# Attenuation parameter is how rapidly weight drops between the two thresholds
# The lowest_weight parameter scales the weights so that everything <= threshold1 has that lowest_weight
weights_from_pvals <- function(pvals, threshold1=1e-4, threshold2=5e-8, attenuation=4, lowest_weight=0.5) {
  w <- -log10(pvals) %>%
    pmax(., -log10(threshold1)) %>%
    pmin(., -log10(threshold2)) %>%
    {. - (-log10(threshold1))} %>%
    {./(-log10(threshold2) - (-log10(threshold1)))}
  w <- w^attenuation * (1-lowest_weight) + (lowest_weight)
  return(w)
}

find_snp_and_connectedness_per_cluster <- function(coloc_groups, studies_to_colocalise, coloc_results) {
  components <- unique(coloc_groups$component)
  snp_and_connectedness_per_cluster <- lapply(components, function(component) {
    group_studies <- coloc_groups$unique_study_id[coloc_groups$component == component]
    
    all_group_data <- lapply(group_studies, function(study_id) {
      return(studies_to_colocalise[[study_id]])
    })
    all_group_data <- do.call(rbind, all_group_data)

    snp <- all_group_data |>
      dplyr::group_by(SNP) |>
      dplyr::summarise(cumulative_lbf = sum(LBF)) |>
      dplyr::slice(which.max(cumulative_lbf)) |>
      dplyr::pull(SNP)

    coloc_results_per_cluster <- coloc_results[
      coloc_results$unique_study_a %in% group_studies &
      coloc_results$unique_study_b %in% group_studies,
    ]

    total_pairs <- nrow(coloc_results_per_cluster)
    total_pairs_with_h4 <- sum(coloc_results_per_cluster$PP.H4.abf >= posterior_prob_h4_threshold, na.rm = TRUE)
    total_pairs_with_h3 <- sum(coloc_results_per_cluster$PP.H3.abf >= posterior_prob_h4_threshold, na.rm = TRUE)
    h4_connectedness <- total_pairs_with_h4 / total_pairs
    h3_connectedness <- total_pairs_with_h3 / total_pairs

    return(data.frame(component = component, snp = snp, h4_connectedness = h4_connectedness, h3_connectedness = h3_connectedness))
  }) |> dplyr::bind_rows()

  return(snp_and_connectedness_per_cluster)
}

main()
