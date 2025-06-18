source('constants.R')
bp_range <- 10000
h4_threshold <- 0.8

parser <- argparser::arg_parser('Colocalise studies per region')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block that the studies are in', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Coloc result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--worker_guid', help = 'Worker GUID', type = 'character', default = NA)
args <- argparser::parse_args(parser)

main <- function() {
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
      dplyr::filter(variant_type == variant_types$common & min_p <= p_value_threshold) |>
      dplyr::arrange(unique_study_id)
  }

  if (!is.na(args$worker_guid)) {
    existing_finemapped_studies_file <- glue::glue('{data_dir}/ld_blocks/{args$ld_block}/finemapped_studies.tsv')
    existing_finemapped_studies <- vroom::vroom(existing_finemapped_studies_file, col_types = finemapped_column_types, show_col_types=F)
    finemapped_studies <- dplyr::bind_rows(finemapped_studies, existing_finemapped_studies) |>
      dplyr::filter(variant_type == variant_types$common & min_p <= p_value_threshold) |>
      dplyr::mutate(file = sub('/local-scratch/projects/genotype-phenotype-map/data/', data_dir, file)) |>
      dplyr::arrange(unique_study_id)
  }

  if (!file.exists(finemapped_file) || nrow(block) == 0 || nrow(finemapped_studies) == 0) {
    message(glue::glue('Nothing to coloc in LD region {ld_info$ld_block_data}, skipping.'))
    vroom::vroom_write(data.frame(), args$completed_output_file)
    return()
  }

  finemapped_studies <- dplyr::arrange(finemapped_studies, unique_study_id)

  # Get all possible pairs within bp_range, excluding already calculated pairs
  study_pairs <- get_study_pairs_to_coloc(finemapped_studies, coloc_results, args$worker_guid)
  message(glue::glue('Found {nrow(study_pairs)} study pairs to coloc for {args$ld_block}'))
  
  if (nrow(study_pairs) == 0) {
    message(glue::glue('No study pairs to coloc in LD region {ld_info$ld_block_data}, skipping.'))
    vroom::vroom_write(data.frame(), args$completed_output_file)
    return()
  }

  studies_to_colocalise <- lapply(finemapped_studies$file,
    function(file) vroom::vroom(file, delim = '\t',
      show_col_types = F,
      col_select = c('SNP', 'BP', 'BETA', 'SE', 'P', 'EAF'),
      altrep = F
    )
  )
  names(studies_to_colocalise) <- finemapped_studies$unique_study_id

  new_coloc_results <- lapply(1:nrow(study_pairs), function(i) {
    start_time <- Sys.time()
    pair <- study_pairs[i,]
    first_gwas <- studies_to_colocalise[[pair$unique_study_a]]
    second_gwas <- studies_to_colocalise[[pair$unique_study_b]]
    first_study <- finemapped_studies[finemapped_studies$unique_study_id == pair$unique_study_a,]
    second_study <- finemapped_studies[finemapped_studies$unique_study_id == pair$unique_study_b,]

    tryCatch({
      result <- pairwise_coloc_analysis(first_gwas, second_gwas, first_study, second_study)
      if (is.null(result)) return(NULL)
    }, error = function(e) {
      stop(glue::glue('Error colocating {pair$unique_study_a} and {pair$unique_study_b}: {e}'))
    })

    result <- dplyr::bind_cols(pair, result)
    time_taken <- Sys.time() - start_time
    message(glue::glue('coloc {result$h4} in {time_taken}'))
    return(result)
  })

  new_coloc_results <- dplyr::bind_rows(new_coloc_results[!sapply(new_coloc_results, is.null)])
  if (nrow(new_coloc_results) > 0) {
    coloc_results <- dplyr::bind_rows(coloc_results, new_coloc_results)
    vroom::vroom_write(coloc_results, coloc_results_file)
  }

  clustered_results <- cluster_coloc_results(coloc_results)
  # vroom::vroom_write(clustered_results$groups, glue::glue('{ld_info$ld_block_data}/coloc_clustered_results.tsv.gz'))
  # Do we also want to write out the igraph objects and pruned studies?
  
  vroom::vroom_write(data.frame(), args$completed_output_file)
}

get_study_pairs_to_coloc <- function(studies, existing_results, worker_guid) {
  pairs <- data.frame(t(utils::combn(studies$unique_study_id, 2))) |>
    dplyr::rename(unique_study_a = X1, unique_study_b = X2) |>
    dplyr::mutate(study_a = sub('_.*', '', unique_study_a), study_b = sub('_.*', '', unique_study_b))
  bp_distances <- utils::combn(studies$bp, 2, function(x) {abs(x[1] - x[2])})
  pairs$bp_distance <- bp_distances
  
  if (!is.na(worker_guid)) {
    pairs <- pairs |>
      dplyr::filter(study_a == worker_guid | study_b == worker_guid)
  }

  pairs <- pairs |>
    dplyr::filter(bp_distance <= bp_range)
  
  if (nrow(existing_results) > 0) {
    pairs <- dplyr::anti_join(pairs, existing_results, by = c("unique_study_a", "unique_study_b")) 
  }

  return(pairs)
}

#' run_coloc_analysis takes two already harmonised gwases, and runs coloc on the results
#' @param first_gwas: first gwas to be run through coloc.  This is the gwas that results will be based off
#' @param second_gwas: second gwas to be run through coloc
#' @returns tibble of coloc results (h0 - h4)
#' @import coloc
#' @import tibble
#' @export
pairwise_coloc_analysis <- function(first_gwas, second_gwas, first_study, second_study) {
  first_gwas <- dplyr::filter(first_gwas, EAF > 0 & EAF < 1) |>
    dplyr::mutate(P = ifelse(P == 0, .Machine$double.xmin, P))
  second_gwas <- dplyr::filter(second_gwas, EAF > 0 & EAF < 1) |>
    dplyr::mutate(P = ifelse(P == 0, .Machine$double.xmin, P))

  harmonised_gwases <- harmonise_gwases(first_gwas, second_gwas)
  if (length(harmonised_gwases) == 0) return(NULL)
  
  first_gwas <- harmonised_gwases[[1]]
  second_gwas <- harmonised_gwases[[2]]

  if (nrow(first_gwas) < 50 || nrow(second_gwas) < 50) return(NULL)

  first_type <- if (first_study$category == study_categories$continuous) "quant" else "cc"
  second_type <- if (second_study$category == study_categories$continuous) "quant" else "cc"

  first_coloc_dataset <- list(
    N = as.numeric(first_study$sample_size),
    type = first_type,
    snp = first_gwas$SNP,
    beta = first_gwas$BETA,
    varbeta = first_gwas$SE^2,
    pvalues = first_gwas$P,
    position = first_gwas$BP,
    MAF = ifelse(first_gwas$EAF < 0.5, first_gwas$EAF, 1-first_gwas$EAF)
  )

  second_coloc_dataset <- list(
    N = as.numeric(second_study$sample_size),
    type = second_type,
    snp = second_gwas$SNP,
    beta = second_gwas$BETA,
    varbeta = second_gwas$SE^2,
    pvalues = second_gwas$P,
    position = second_gwas$BP,
    MAF = ifelse(second_gwas$EAF < 0.5, second_gwas$EAF, 1-second_gwas$EAF)
  )

  result <- coloc::coloc.abf(dataset1 = first_coloc_dataset, dataset2 = second_coloc_dataset)
  coloc_results <- tibble::tribble(
    ~h0, ~h1, ~h2, ~h3, ~h4,
    result$summary[2], result$summary[3], result$summary[4], result$summary[5], result$summary[6]
  )

  return(coloc_results)
}

harmonise_gwases <- function(...) {
  gwases <- list(...)
  snpids <- Reduce(intersect, lapply(gwases, function(gwas) gwas$SNP))
  if (length(snpids) <= 1) return(list())

  gwases <- lapply(gwases, function(gwas) {
    dplyr::filter(gwas, SNP %in% snpids & !duplicated(SNP)) |>
      dplyr::arrange(SNP)
  })

  return(gwases)
}

cluster_coloc_results <- function(coloc_results) {
  h4_adj_mx <- make_adjacency_matrix(coloc_results = coloc_results)

  # Generate graph from adjacency matrix
  g <- igraph::graph_from_adjacency_matrix(h4_adj_mx, mode="undirected", weighted=TRUE, diag=FALSE)
  # Extract components (groups of linked studies)
  g_comps <- igraph::components(g)
  # Prune singleton studies (vertices) with no connections
  vert_out <- igraph::V(g)[g_comps$membership %in% which(g_comps$csize == 1)]
  g2 <- igraph::delete_vertices(g, vert_out)
  # Rederive components
  g2_comps <- igraph::components(g2)
  # Split components into subgraphs
  g2_subgraphs <- lapply(seq(1:g2_comps$no), function(x){
    sg <- igraph::induced_subgraph(g2, vids = igraph::V(g2)[g2_comps$membership == x])
    igraph::V(sg)$group <- x
    return(sg)
  })
  
  # Calculate edge betweenenss for each subgraph and remove edges to maximise modularity
  g2_pruned <- lapply(g2_subgraphs, function(x){
    eb <- igraph::cluster_edge_betweenness(x)
    if(length(unique(eb$membership)) == 1){
      igraph::V(x)$ebc_group <- igraph::V(x)$group
      igraph::V(x)$component <- igraph::V(x)$group
      return(x)
    } else {
      # Delete edges removed to achieve the maximum modularity score
      x <- igraph::delete_edges(x, edges = eb$removed.edges[seq(1:which.max(eb$modularity))])
      # Resulting study cluster membership
      igraph::V(x)$ebc_group <- paste(igraph::V(x)$group, eb$membership, sep = ".")
      # Prune resulting single studies
      grp_single <- as.numeric(which(table(eb$membership) == 1))
      vert_out2 <- igraph::V(x)[eb$membership %in% grp_single]
      x <- igraph::delete_vertices(x, vert_out2)
      # Label modules unlinked after edge removal (components)
      x_comps <- igraph::components(x)
      igraph::V(x)$component <- paste(igraph::V(x)$group, x_comps$membership, sep = ".")
      return(x)
    }
  })
  
  # Recombine subgraphs
  g_out <- igraph::disjoint_union(g2_pruned)
  # Study group membership
  grps <- factor(igraph::V(g_out)$ebc_group)
  levels(grps) <- seq(1:length(unique(grps)))
  # Study component membership
  comp <- factor(igraph::V(g_out)$component)
  levels(comp) <- seq(1:length(unique(comp)))
  
  df_out <- data.frame(study_id = igraph::vertex_attr(g_out)$name, ebc_group = grps, component = comp)
  
  # Pruned studies with no module membership
  pruned <- setdiff(igraph::vertex_attr(g2)$name, igraph::vertex_attr(g_out)$name)
  
  # Output original igraph obj, clustered and pruned igraph object, cluster membership and studies pruned
  return(list(igraph_obj = g, igraph_obj_ebc = g_out, groups = df_out, pruned_studies = pruned))
}

make_adjacency_matrix <- function(coloc_results) {
  # Make symmetrical adjacency matrix
  h4_adj_mx <- rbind(
    coloc_results |> dplyr::select(a = unique_study_a, b = unique_study_b, h4 = h4 ),
    coloc_results |> dplyr::select(a = unique_study_b, b = unique_study_a, h4 = h4 )) |>
    tidyr::pivot_wider(names_from = a, values_from = h4) |>
    tibble::column_to_rownames("b") |>
    as.matrix()
  
  h4_adj_mx <- h4_adj_mx[,rownames(h4_adj_mx)]
  diag(h4_adj_mx) <- 0
  h4_adj_mx[is.na(h4_adj_mx)] <- 0
  h4_adj_mx[h4_adj_mx < h4_threshold] <- 0

  return(h4_adj_mx)
}

main()
