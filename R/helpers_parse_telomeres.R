if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("chr", "score", "motif"))
}

#' Process telomere motif hits for one chromosome side
#'
#' @param file Data frame. Contains columns motif and chrscores, typically from a parsed telomere file.
#' @param side Character string. Indicates whether the hits are from the "left" or "right" telomeric region.
#'
#' @return A tibble with columns: `motif`, `side`, `contig`, `score`, and `rc` (reverse complement).
#' @keywords internal
process_telomere_side <- function(file, side) {
    df <- file
    tibble::tibble(
        motif = df$motif,
        side = side,
        chr = strsplit(df$chrscores, ";")
    ) %>%
        tidyr::unnest(chr) %>%
        tidyr::separate_wider_delim(chr, names = c("contig", "score"), delim = "|") %>%
        dplyr::filter(as.numeric(score) >= 50) %>%
        dplyr::mutate(rc = sapply(motif, function(x) as.character(Biostrings::reverseComplement(Biostrings::DNAString(x)))))
}

#' Building telomere data frame with left and right telomeres
#'
#' @param contigs Character vector. List of contig names.
#' @param left Data frame. Output from process_telomere_side of the left side.
#' @param right Data frame. Output from process_telomere_side of the right side.
#'
#' @return A data frame with telomeric motif and score information for each contig side.
#' @keywords internal
build_telomere_dataframe <- function(contigs, left, right) {
    telomeres_out <- lapply(contigs, function(contig) {
        l <- left[left$contig == contig, ]
        r <- right[right$contig == contig, ]

        data.frame(
            contig = contig,
            left_motif <- if (nrow(l) > 0) l$motif[1] else "0",
            left_score <- if (nrow(l) > 0) l$score[1] else "0",
            right_motif <- if (nrow(r) > 0) r$motif[1] else "0",
            right_rc    <- if (nrow(r) > 0) r$rc[1]    else "0",
            right_score <- if (nrow(r) > 0) r$score[1] else "0",
            stringsAsFactors = FALSE
        )
    })

    df <- do.call(rbind, telomeres_out)
    rownames(df) <- make.unique(df$contig)
    df$contig <- NULL
    colnames(df) <- c("leftmotif", "leftscore", "rightmotif", "rightrc", "rightscore")
    df$leftscore <- as.numeric(df$leftscore)
    df$rightscore <- as.numeric(df$rightscore)

    return(df)
}

#' Classifying telomeres based on motif match strength
#'
#' @param df Data frame. Output from build_telomere_dataframe, containing left and right telomeric motifs, scores, and reverse complements.
#'
#' @return The same data frame with an added column `telcomp` indicating:
#' \describe{
#'   \item{0}{No confident telomere motif}
#'   \item{1}{Partial match (one side matches known telomere motifs)}
#'   \item{2}{Confident match (both sides match same or reverse-complement motifs)}
#' }
#' @keywords internal
classify_telomeres <- function(df) {
    df$telcomp <- 0

    both_match <- with(df, leftmotif != "0" & (leftmotif == rightmotif | leftmotif == rightrc))
    df$telcomp[both_match] <- 2

    motifs <- unique(c("CCCTAA", "TTAGGG", df$leftmotif[both_match], df$rightmotif[both_match], df$rightrc[both_match]))

    partial_match <- with(df, telcomp != 2 & (leftmotif %in% motifs | rightmotif %in% motifs))
    df$telcomp[partial_match] <- 1

    weak_signal <- with(df, telcomp == 0 & (leftmotif != "0" | rightmotif != "0"))
    df$telcomp[weak_signal]

    return(df)
}
