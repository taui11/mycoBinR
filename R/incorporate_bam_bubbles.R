#' Title
#'
#' @param DATA TODO
#' @param telomeres TODO
#' @param paths_df TODO
#'
#' @return TODO
#' @export
#'
#' @examples \dontrun{TODO}
incorporate_bam_bubbles <- function(DATA, telomeres, paths_df) {
    include_ranges <- create_include_ranges(DATA, telomeres)
    to_net <- load_supplementary_reads(paths_df$path[paths_df$name == "bam_file"])
    to_net <- filter_by_include_ranges(to_net, include_ranges, DATA)
    shared_matrix <- calculate_shared_matrix(to_net)
    m1_final <- detect_cluster(shared_matrix, DATA)
    backup_DATA <- DATA
    DATA <- refine_taxonomic_assignments(m1_final, DATA)

    knitr::kable(table(backup_DATA$tax_cons, paste("new", DATA$tax_cons, sep = "_")))

    return(DATA)
}
