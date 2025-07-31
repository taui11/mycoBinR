if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("chr", "score", "motif"))
}

#' Title
#'
#' @param file TODO
#' @param side TODO
#'
#' @return TODO
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

#' Title
#'
#' @param contigs TODO
#' @param left TODO
#' @param right TODO
#'
#' @return TODO
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

#' Title
#'
#' @param df TODO
#'
#' @return TODO
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
