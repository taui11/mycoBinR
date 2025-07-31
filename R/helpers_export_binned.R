#' Create output directories
#'
#' @param fasta_folder Character string. Path to the directory where binned FASTA files will be saved.
#' @param bed_folder Character string. Path to the directory where binned BED files will be saved.
#'
#' @return NULL. Creates two new directories.
#' @keywords internal
create_output_dirs <- function(fasta_folder, bed_folder) {
    dir.create(fasta_folder, recursive = TRUE, showWarnings = FALSE)
    dir.create(bed_folder, showWarnings = FALSE)
}


#' Building data frame for taxonomic categories
#'
#' @param data Data frame. Contains at leas tax_cons and len column. Typically from incorporate_bam_bubbles.
#'
#' @return A data frame with columns:
#'      \describe{
#'          \item{chrom}{Chromosome/contig names (rownames of input).}
#'          \item{chromStart}{Start position, always 0.}
#'          \item{chromEnd}{End position, taken from the `len` column of input data.}
#'          \item{bin}{Factor vector representing bins labeled by taxonomic consensus.}
#'      }
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

#' Export FASTA sequences grouped by taxonomic bins
#'
#' @param fasta List of class "DNAbin" or "AAbin".
#' @param categories Data frame. Contains at least a bin column. Typically from build_categories_df.
#' @param fasta_folder Character string. Path to the binned FASTA directory.
#' @param sample_name Character string. Sample/Title identifier of the data (e.g. "001").
#'
#' @return NULL. Saves binned FASTA files in the fasta_folder directory.
#' @keywords internal
export_fasta_bins <- function(fasta, categories, fasta_folder, sample_name) {
    for (bin_level in levels(categories$bin)) {
        contigs <- categories$chrom[categories$bin == bin_level]
        out_fasta <- file.path(fasta_folder, paste0(sample_name, "_", bin_level, ".fa"))
        ape::write.FASTA(fasta[contigs], out_fasta)
    }
}

#' EXPORT BED files grouped by taxonomic bins
#'
#' @param categories Data frame. Contains at least a bin column. Typically from build_categories_df.
#' @param bed_folder Character string. Path to the binned BED directory.
#' @param sample_name Character string. Sample/Title identifier of the data (e.g. "001").
#'
#' @return NULL. Saves binned BED files in the bed_folder directory.
#' @keywords internal
export_bed_bins <- function(categories, bed_folder, sample_name) {
    for (bin_level in levels(categories$bin)) {
        bed_data <- categories[categories$bin == bin_level, 1:3]
        out_bed <- file.path(bed_folder, paste0(sample_name, "_", bin_level, ".bed"))
        utils::write.table(bed_data, out_bed, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    }
}
