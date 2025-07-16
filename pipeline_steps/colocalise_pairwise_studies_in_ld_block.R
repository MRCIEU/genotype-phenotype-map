source('constants.R')
bp_range <- 10000
h4_threshold <- 0.8

parser <- argparser::arg_parser('Colocalise studies per region')
parser <- argparser::add_argument(parser, '--ld_block', help = 'LD block that the studies are in', type = 'character')
parser <- argparser::add_argument(parser, '--completed_output_file', help = 'Coloc result file to save', type = 'character')
parser <- argparser::add_argument(parser, '--worker_guid', help = 'Worker GUID', type = 'character', default = NA)
args <- argparser::parse_args(parser)

# For testing only
#coloc_results_file <- "/local-scratch/projects/genotype-phenotype-map/data/ld_blocks/EUR/22/35695831-37164130/coloc_pairwise_results.tsv.gz"
#finemapped_file <- "/local-scratch/projects/genotype-phenotype-map/data/ld_blocks/EUR/22/35695831-37164130/finemapped_studies.tsv"

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
    message(glue::glue('Nothing to coloc in LD region {ld_info$ld_block_data}, skipping.'))
    vroom::vroom_write(data.frame(), args$completed_output_file)
    return()
  }

  finemapped_studies <- dplyr::arrange(finemapped_studies, unique_study_id)

  # Get all possible pairs within bp_range, excluding already calculated pairs
  study_pairs <- get_study_pairs_to_coloc(finemapped_studies, coloc_results, args$worker_guid)
  message(glue::glue('Found {nrow(study_pairs)} study pairs to coloc for {args$ld_block} in {diff_time_taken(start_time)}'))

  if (nrow(study_pairs) == 0) {
    message(glue::glue('No study pairs to coloc in LD region {ld_info$ld_block_data}, skipping.'))
    vroom::vroom_write(data.frame(), args$completed_output_file)
    return()
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
  message(glue::glue('Loaded {length(studies_to_colocalise)} studies in {diff_time_taken(start_time)}'))

  new_coloc_results <- lapply(1:nrow(study_pairs), function(i) {
    start_time <- Sys.time()
    pair <- study_pairs[i,]
    first_gwas <- studies_to_colocalise[[pair$unique_study_a]]
    second_gwas <- studies_to_colocalise[[pair$unique_study_b]]
    first_study <- finemapped_subset[finemapped_subset$unique_study_id == pair$unique_study_a,]
    second_study <- finemapped_subset[finemapped_subset$unique_study_id == pair$unique_study_b,]

    tryCatch({
      result <- pairwise_coloc_analysis(first_gwas, second_gwas, first_study, second_study)
      if (is.null(result)) return(NULL)
    }, error = function(e) {
      stop(glue::glue('Error colocating {pair$unique_study_a} and {pair$unique_study_b}: {e}'))
    })

    result <- dplyr::bind_cols(pair, result)
    return(result)
  })

  message(glue::glue('Colocated {nrow(new_coloc_results)} study pairs in {diff_time_taken(start_time)}'))
  print(warnings())

  new_coloc_results <- dplyr::bind_rows(new_coloc_results[!sapply(new_coloc_results, is.null)])
  if (nrow(new_coloc_results) > 0) {
    coloc_results <- dplyr::bind_rows(coloc_results, new_coloc_results) |>
      dplyr::distinct(unique_study_a, unique_study_b, .keep_all = TRUE) |>
      dplyr::mutate(ld_block = args$ld_block)

    # Check for within study colocalising finemapped regions (requires pairwise coloc to between credible sets)
    poor_finemapping <- nrow(coloc_results |> dplyr::filter(study_a == study_b & h4 >= 0.5)) != 0

    if (poor_finemapping) {
      ignore_regions <- prune_finemapped(finemapped_studies, coloc_results)
      # Remove dodgy finemapped regions from coloc_results and mark as ignored in finemapped_studies
      coloc_results <- coloc_results |>
        dplyr::mutate(ignore = ignore | unique_study_a %in% ignore_regions | unique_study_b %in% ignore_regions)

      finemapped_studies <- finemapped_studies |>
        dplyr::mutate(ignore = ignore | unique_study_id %in% ignore_regions)

      message(glue::glue('Number of studies to ignore: {sum(finemapped_studies$ignore)}'))
      message(glue::glue('Number of coloc results to ignore: {sum(coloc_results$ignore)}'))
      vroom::vroom_write(finemapped_studies, finemapped_file)
    } 
    vroom::vroom_write(coloc_results, coloc_results_file)

    clustered_results <- cluster_coloc_results(coloc_results, start_time)
    snp_per_cluster <- find_snp_per_cluster(clustered_results$groups, studies_to_colocalise)
    clustered_results$groups <- clustered_results$groups |>
      dplyr::mutate(ld_block = args$ld_block) |>
      dplyr::left_join(snp_per_cluster, by = "component") |>
      dplyr::arrange(component)

    vroom::vroom_write(clustered_results$groups, glue::glue('{ld_info$ld_block_data}/coloc_clustered_results.tsv.gz'))
    saveRDS(clustered_results$pruned_studies, glue::glue('{ld_info$ld_block_data}/coloc_pruned_studies.rds'))
    saveRDS(clustered_results$pruned_igraph_obj, glue::glue('{ld_info$ld_block_data}/coloc_igraph_obj.rds')) # Return the pruned object instead?
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

  existing_results <- data.table::as.data.table(existing_results)

  pairs_filtered <- studies[ studies, on = .(id < id), allow.cartesian = TRUE ][
    , bp_distance := abs(i.bp - bp) ][
    bp_distance <= bp_range | i.study == study ][
    , .(
      unique_study_a = unique_study_id,
      study_a = study,
      unique_study_b = i.unique_study_id,
      study_b = i.study,
      bp_distance = bp_distance,
      ignore = F
    )
  ] |>
    tibble::as_tibble()

  if (!is.na(worker_guid)) {
    pairs_filtered <- pairs_filtered |>
      dplyr::filter(study_a == worker_guid | study_b == worker_guid)
  }

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
pairwise_coloc_analysis <- function(first_gwas, second_gwas, first_study, second_study) {
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
  snpids <- Reduce(intersect, lapply(gwases, function(gwas) gwas$SNP))
  if (length(snpids) <= 1) return(list())

  gwases <- lapply(gwases, function(gwas) {
    gwas <- data.table::as.data.table(gwas)
    gwas[SNP %in% snpids & !duplicated(SNP), ][order(SNP)]
    return(gwas)
  })

  return(gwases)
}

## Need to check the performance of this when all finemapped regions per study are added to pairwise coloc results
## Check matching of finemapped_studies and coloc_results unique IDs with updated convention (currently mismatched for some regions)
prune_finemapped <- function(finemapped_studies, coloc_results){
  # Split colocalising finemapped region by study
  finemapped_colocs <- coloc_results |> 
    dplyr::filter(study_a == study_b & h4 >= 0.5) |>
    dplyr::mutate(
      min_p_study_a = finemapped_studies[match(unique_study_a, finemapped_studies$unique_study_id),]$min_p,
      min_p_study_b = finemapped_studies[match(unique_study_b, finemapped_studies$unique_study_id),]$min_p
    ) |>
    dplyr::group_by(study_a) |>
    dplyr::group_split()
  
  remove_regions <- lapply(finemapped_colocs, function(study_regions){
    # Store minimum p-value per finemapped region
    pvals <- data.frame(
      study = c(study_regions$unique_study_a, study_regions$unique_study_b),
      min_p = c(study_regions$min_p_study_a, study_regions$min_p_study_b)) |> unique()
    
    regions <- pvals$study
    n_regions <- length(regions)
    n_links <- nrow(study_regions)
    n_max_links <- n_regions*(n_regions - 1)/2
    all_linked <- n_links == n_max_links
    
    # If all problematic finemapped regions colocalise with each other (maximum connectivity) then retain the region with the smallest min_p
    # Else, if several groups of colocalsing finemapped regions exist for a study, or regions show different connectivity with each other,
    # retain the region with the smalled min_p in each group by ittertively pruning regions that are most connected (splitting ties by min_p)
    
    if(all_linked == TRUE){
      keep <- pvals$study[which.min(pvals$min_p)]
    } else {
      links <- data.frame(
        table(study = c(study_regions$unique_study_a, study_regions$unique_study_b))) |>
        dplyr::left_join(pvals, by = "study") |>
        dplyr::arrange(desc(Freq), desc(min_p)) |> unique()
      
      study_regions_pruned <- study_regions
      
      while(any(links$Freq > 1)){
        study_regions_pruned <- study_regions_pruned |> 
          dplyr::filter(!(unique_study_a %in% links$study[1]) & !(unique_study_b %in% links$study[1]))
        
        if(nrow(study_regions_pruned) == 0){
          keep <- links$study[-1]
          break
        }
        
        links <- data.frame(
          table(study = c(study_regions_pruned$unique_study_a, study_regions_pruned$unique_study_b))) |>
          dplyr::left_join(pvals, by = "study") |>
          dplyr::arrange(desc(Freq), desc(min_p)) |> unique()
      }
      
      if(nrow(study_regions_pruned) != 0){
        keep <- apply(study_regions_pruned, 1, function(region_pair){
          min_ps <- region_pair[c("min_p_study_a","min_p_study_b")]
          regions <- region_pair[c("unique_study_a","unique_study_b")]
          regions[which.min(min_ps)]
        })
      }
    }
    
    remove <- regions[!(regions %in% keep)]
    return(remove)
  })
  
  remove_regions <- unlist(remove_regions)
  return(remove_regions)
}

cluster_coloc_results <- function(coloc_results, start_time) {
  coloc_results <- coloc_results |> dplyr::filter(!ignore)
  message(glue::glue('Clustering {nrow(coloc_results)} coloc results starting {diff_time_taken(start_time)}'))
  h4_adj_mx <- make_adjacency_matrix(coloc_results = coloc_results)
  
  # Generate graph from adjacency matrix
  g <- igraph::graph_from_adjacency_matrix(h4_adj_mx, mode="undirected", weighted=TRUE, diag=FALSE)
  # Extract components (groups of linked studies)
  g_comps <- igraph::components(g)
  # Prune singleton studies (vertices) with no connections
  vert_out <- igraph::V(g)[g_comps$membership %in% which(g_comps$csize == 1)]
  g2 <- igraph::delete_vertices(g, vert_out)
  
  message(glue::glue('Pruning graph starting {diff_time_taken(start_time)}'))
  
  # Calculate edge betweenenss and remove edges to maximise modularity
  eb <- igraph::cluster_edge_betweenness(g2)
  
  g2_pruned <- g2
  if(length(unique(eb$membership)) == 1){
    igraph::V(g2_pruned)$component <- 1
  } else {
    # Delete edges removed to achieve the maximum modularity score
    g2_pruned <- igraph::delete_edges(g2_pruned, edges = eb$removed.edges[seq(1:which.max(eb$modularity))])
    # Recalculate components (disconnected modules)
    g2_pruned_comps <- igraph::components(g2_pruned)
    igraph::V(g2_pruned)$component <- g2_pruned_comps$membership
    # Prune new resulting singleton studies (vertices) with no connections
    vert_out2 <- igraph::V(g2_pruned)[g2_pruned_comps$membership %in% which(g2_pruned_comps$csize == 1)]
    g2_pruned <- igraph::delete_vertices(g2_pruned, vert_out2)
  }
  
  message(glue::glue('Outputting results in {diff_time_taken(start_time)}'))
  coloc_groups <- data.frame(study_id = igraph::vertex_attr(g2_pruned)$name, component = igraph::V(g2_pruned)$component)
  
  # Pruned studies with no module membership after edge betweenness clustering (this does not include singleton studies that show no H4 > threshold removed above)
  pruned <- setdiff(igraph::vertex_attr(g2)$name, igraph::vertex_attr(g2_pruned)$name)
  
  # Output original igraph object, pruned igraph object, cluster membership and studies pruned
  return(list(igraph_obj = g, pruned_igraph_obj = g2_pruned, groups = coloc_groups, pruned_studies = pruned))
}

make_adjacency_matrix <- function(coloc_results) {
  # Make symmetrical adjacency matrix
  h4_adj_mx <- rbind(
    coloc_results |> 
      dplyr::select(a = unique_study_a, b = unique_study_b, h4 = h4) |>
      dplyr::group_by(a, b) |>
      dplyr::summarise(h4 = max(h4)),
    coloc_results |> 
      dplyr::select(a = unique_study_b, b = unique_study_a, h4 = h4) |>
      dplyr::group_by(a, b) |>
      dplyr::summarise(h4 = max(h4))
  ) |>
    tidyr::pivot_wider(names_from = a, values_from = h4) |>
    tibble::column_to_rownames("b") |>
    as.matrix()
  
  h4_adj_mx <- h4_adj_mx[,rownames(h4_adj_mx)]
  diag(h4_adj_mx) <- 0
  h4_adj_mx[is.na(h4_adj_mx)] <- 0
  h4_adj_mx[h4_adj_mx < h4_threshold] <- 0

  return(h4_adj_mx)
}

find_snp_per_cluster <- function(coloc_groups, studies_to_colocalise) {
  components <- unique(coloc_groups$component)
  snp_per_cluster <- sapply(components, function(component) {
    group_studies <- coloc_groups$study_id[coloc_groups$component == component]
    
    all_group_data <- lapply(group_studies, function(study_id) {
      return(studies_to_colocalise[[study_id]])
    })
    all_group_data <- do.call(rbind, all_group_data)

    snp <- all_group_data |>
      dplyr::group_by(SNP) |>
      dplyr::summarise(cumulative_lbf = sum(LBF)) |>
      dplyr::slice(which.max(cumulative_lbf)) |>
      dplyr::pull(SNP)

    return(snp)
  })

  return(data.frame(component = components, snp = snp_per_cluster))
}

main()
