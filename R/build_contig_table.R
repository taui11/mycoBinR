#' Title
#'
#' @param paths_df TODO
#' @param API_KEY TODO
#'
#' @returns TODO
#' @export
#'
#' @examples TODO
build_contig_table <- function(paths_df, API_KEY = API_KEY) {
    report <- load_assembly_info(paths_df)
    a <- load_dna_sequences(paths_df)
    nombres <- names(a)
    cov <- load_coverage_data(paths_df)
    gc <- compute_gc_content(a)
    buscos <- load_busco_data(paths_df)

    data <- build_initial_df(report, nombres, cov, gc)
    data <- enrich_with_graph_info(data, report)
    data <- add_busco_data(data, buscos)
    data <- cluster_data(data)
    data <- flag_busco_completeness(data)

    taxonomy_list <- process_taxonomy(paths_df$path[paths_df$name == "taxonomy"], API_KEY = API_KEY_GLOBAL)
    consensus_ids <- compute_consensus_clusters(taxonomy_list)
    data <- assign_consensus_clusters(data, consensus_ids)
    classification <- extract_taxonomic_labels(data, taxonomy_list)
    data <- cbind(data, classification)
    data <- summarize_taxonomy(data)
    return(data)
}
