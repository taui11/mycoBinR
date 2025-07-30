#' Title
#'
#' @param DATA
#' @param paths_df
#' @param SAMPLE
#'
#' @returns
#' @export
#'
#' @examples
export_binned <- function(DATA, paths_df, SAMPLE) {
    fasta <- read.FASTA(paths_df$path[paths_df$name == "assembly"])

    fasta_folder <- paths_df$path[paths_df$name == "fasta_folder"]
    bed_folder <- paths_df$path[paths_df$name == "bed_folder"]

    create_output_dirs(fasta_folder, bed_folder)
    categories <- build_categories_df(DATA)

    export_fasta_bins(fasta, categories, fasta_folder, SAMPLE)
    export_bed_bins(categories, bed_folder, SAMPLE)
}
