#' Title
#'
#' @param fasta_folder TODO
#' @param bed_folder TODO
#'
#' @returns TODO
#' @keywords internal
create_output_dirs <- function(fasta_folder, bed_folder) {
    dir.create(fasta_folder, recursive = TRUE, showWarnings = FALSE)
    dir.create(bed_folder, showWarnings = FALSE)
}


#' Title
#'
#' @param data TODO
#'
#' @returns TODO
#' @keywords internal
build_categories_df <- function(data) {
    data$bin_label <- paste0("bin_", data$tax_cons)
    data.frame(
        chrom = rownames(data),
        chromStart = 0,
        chromEnd = data$len,
        bin = factor(data$bin_label)
    )
}

#' Title
#'
#' @param fasta TODO
#' @param categories TODO
#' @param fasta_folder TODO
#' @param sample_name TODO
#'
#' @returns TODO
#' @keywords internal
export_fasta_bins <- function(fasta, categories, fasta_folder, sample_name) {
    for (bin_level in levels(categories$bin)) {
        contigs <- categories$chrom[categories$bin == bin_level]
        out_fasta <- file.path(fasta_folder, paste0(sample_name, "_", bin_level, ".fa"))
        write.FASTA(fasta[contigs], out_fasta)
    }
}

#' Title
#'
#' @param categories TODO
#' @param bed_folder TODO
#' @param sample_name TODO
#'
#' @returns TODO
#' @keywords internal
export_bed_bins <- function(categories, bed_folder, sample_name) {
    for (bin_level in levels(categories$bin)) {
        bed_data <- categories[categories$bin == bin_level, 1:3]
        out_bed <- file.path(bed_folder, paste0(sample_name, "_", bin_level, ".bed"))
        write.table(bed_data, out_bed, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    }
}
