#' Build contig-level data table from input paths and taxonomy
#'
#' This function integrates multiple sources of contig-level information including
#' assembly reports, DNA sequences, coverage data, GC content, BUSCO results, and taxonomy.
#' It performs clustering and classification steps to enrich the contig table with biological
#' and taxonomic context.
#'
#' @param paths_df Data frame. Must contain a column named `path` with file paths to input data,
#'        and a column `name` to identify the type of each file (e.g., "taxonomy"). Typically from create_filepaths_df.
#' @param api_key Character string or NULL. NCBI Entrez API key needed for fetching taxonomy.
#'
#' @returns A data frame with one row per contig, enriched with computed features, clustering,
#'          BUSCO information, taxonomic classification, and summary labels.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' paths_df <- data.frame(
#'   name = c("assembly", "coverage", "busco", "taxonomy"),
#'   path = c("assembly_report.txt", "coverage.tsv", "busco.txt", "taxonomy.tsv")
#' )
#' build_contig_table(paths_df, API_KEY = "your_ncbi_api_key")
#' }
build_contig_table <- function(paths_df, api_key = NULL) {
    report <- load_assembly_info(paths_df)
    dna_seq <- load_dna_sequences(paths_df)
    contig_names <- names(dna_seq)
    cov <- load_coverage_data(paths_df)
    gc <- compute_gc_content(dna_seq)
    buscos <- load_busco_data(paths_df)

    data <- build_initial_df(report, contig_names, cov, gc)
    data <- enrich_with_graph_info(data, report)
    data <- add_busco_data(data, buscos)
    data <- cluster_data(data)
    data <- flag_busco_completeness(data)

    taxonomy_list <- process_taxonomy(paths_df$path[paths_df$name == "taxonomy"], api_key)
    consensus_ids <- compute_consensus_clusters(taxonomy_list)
    data <- assign_consensus_clusters(data, consensus_ids)
    classification <- extract_taxonomic_labels(data, taxonomy_list)
    data <- cbind(data, classification)
    data <- summarize_taxonomy(data)
    return(data)
}
