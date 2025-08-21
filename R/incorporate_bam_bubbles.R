#' Title
#'
#' @param DATA Data frame. The main data dable containing contig information. Row names must be contig names. Typically from build_contig_table.
#' @param telomeres Data frame. Telomere information per contig, typically produced by `parse_telomeres()`. Must include columns for left and right motifs, scores, and telomere completeness.
#' @param paths_df Data frame. Must contain columns `name` and `path`. Typically from create_filepaths_df.
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
