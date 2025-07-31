#' Export binned FASTA and BED files based on taxonomic categories
#'
#' This function exports contig-level sequence and annotation data into
#' separate FASTA and BED files, grouped by taxonomic bin. It organizes the data
#' based on a classification scheme and writes each bin’s sequences and regions
#' to separate files for downstream analysis or visualization.
#'
#' @param DATA Data frame. Contains taxonomic and positional data used to create categories.
#' @param paths_df Data frame. Must include columns \code{name} and \code{path} specifying locations of input and output directories/files. Typically from create_filepaths_df.
#' @param SAMPLE Character string. Sample or dataset identifier used for naming output files.
#'
#' @returns NULL. Wites binned FASTA and BED files to specified directories.
#' @export
#'
#' @examples
#' \dontrun{
#' paths <- data.frame(
#'     name = c("assembly", "fasta_folder", "bed_folder"),
#'     path = c("path/to/assembly.fasta", "path/to/fasta_bins", "path/to/bed_bins"),
#'     stringsAsFactors = FALSE
#' )
#' export_binned(DATA, paths, "001")
#' }
export_binned <- function(DATA, paths_df, SAMPLE) {
    fasta <- ape::read.FASTA(paths_df$path[paths_df$name == "assembly"])

    fasta_folder <- paths_df$path[paths_df$name == "fasta_folder"]
    bed_folder <- paths_df$path[paths_df$name == "bed_folder"]

    create_output_dirs(fasta_folder, bed_folder)
    categories <- build_categories_df(DATA)

    export_fasta_bins(fasta, categories, fasta_folder, SAMPLE)
    export_bed_bins(categories, bed_folder, SAMPLE)
}
