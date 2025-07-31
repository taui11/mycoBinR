#' Title
#'
#' @param DATA TODO
#' @param paths_df TODO
#'
#' @return TODO
#' @export
#'
#' @examples \dontrun{TODO}
parse_telomeres <- function(DATA, paths_df) {
    left_file <- utils::read.delim(paths_df$path[paths_df$name == "telomere_left"], header = TRUE)
    right_file <- utils::read.delim(paths_df$path[paths_df$name == "telomere_right"], header = TRUE)

    left <- process_telomere_side(left_file, "left")
    right <- process_telomere_side(right_file, "right")

    telomeres_df <- build_telomere_dataframe(rownames(DATA), left, right)
    telomeres_df <- classify_telomeres(telomeres_df)

    return(telomeres_df)
}
