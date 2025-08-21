#' Incorporate supplementary read connections to refine taxonomic assignments
#'
#' This function integrates supplementary read alignments from BAM files to detect contigs that are likely connected via sequencing reads.
#' It builds inclusion ranges based on telomere motifs, filters reads to those mapping within these ranges, calculates shared read support between contig ends, detects clusters of connected contigs,
#' and updates the taxonomic assignments (`tax_cons`) of contigs to reflect cluster-level consensus.
#'
#' @param DATA Data frame. The main data dable containing contig information. Row names must be contig names. Typically from build_contig_table.
#' @param telomeres Data frame. Telomere information per contig, typically produced by `parse_telomeres()`. Must include columns for left and right motifs, scores, and telomere completeness.
#' @param paths_df Data frame. Must contain columns `name` and `path`. Typically from create_filepaths_df.
#'
#' @return A data frame identical to `DATA`, but with updated taxonomic assignments (`tax_cons`) where clusters of contigs supported by supplementary reads have a consensus assignment.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' paths_df <- data.frame(
#'     name = c("assembly_info", "coverage", "busco", "taxonomy", "bam_file"),
#'     path = c("assembly_info.txt", "coverage.tsv", "busco.txt", "taxonomy.tsv", "reads.bam"),
#'     stringsAsFactors = FALSE
#' )
#' DATA <- build_contig_table(paths_df, API_KEY = "your_ncbi_api_key")
#' telomeres <- parse_telomeres(DATA, paths_df)
#' DATA <- incorporate_bam_bubbles(DATA, telomeres, paths_df)
#' }
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
