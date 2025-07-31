if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("pos", "sl", "el", "sr", "er", "value"))
}

#' Title
#'
#' @param DATA TODO
#' @param telomeres TODO
#'
#' @return TODO
#' @keywords internal
create_include_ranges <- function(DATA, telomeres) {
    include_ranges <- data.frame(
        sl = 1,
        el = pmin(10000, DATA$len),
        sr = pmax(DATA$len - 10000, 1),
        er = DATA$len,
        row.names = rownames(DATA)
    )

    left_telomeres <- rownames(telomeres)[telomeres$leftmotif %in% c("TTAGGG", "CCCTAA")]
    right_telomeres <- rownames(telomeres)[telomeres$rightmotif %in% c("TTAGGG", "CCCTAA")]

    include_ranges[left_telomeres, c("sl", "el")] <- 0
    include_ranges[right_telomeres, c("sr", "er")] <- 0

    duplicated_start <- include_ranges$sr == 1
    include_ranges[duplicated_start, c("sr", "er")] <- 0

    overlaps <- include_ranges$el > include_ranges$sr & include_ranges$sr != 0
    include_ranges[overlaps, c("el", "sr", "er")] <- data.frame(
        el = include_ranges$er[overlaps],
        sr = 0,
        er = 0
    )

    return(include_ranges)
}

#' Title
#'
#' @param bam_file TODO
#'
#' @return TODO
#' @keywords internal
load_supplementary_reads <- function(bam_file) {
    bam_data <- Rsamtools::scanBam(Rsamtools::BamFile(bam_file))[[1]]

    supplementary_flags <- c(2048, 2064)
    supp_reads <- bam_data$flag %in% supplementary_flags
    supp_qnames <- bam_data$qname[supp_reads]

    is_supp_qname <- bam_data$qname %in% supp_qnames

    to_net <- data.frame(
        qname  = bam_data$qname[is_supp_qname],
        rname  = bam_data$rname[is_supp_qname],
        strand = bam_data$strand[is_supp_qname],
        pos    = bam_data$pos[is_supp_qname],
        mapq   = bam_data$mapq[is_supp_qname],
        tlen   = bam_data$isize[is_supp_qname]
    )

    rm(bam_data)
    return(to_net)
}

#' Title
#'
#' @param to_net TODO
#' @param include_ranges TODO
#' @param DATA TODO
#'
#' @return TODO
#' @keywords internal
filter_by_include_ranges <- function(to_net, include_ranges, DATA) {
    include_ranges <- include_ranges %>% tibble::rownames_to_column(var = "contig")

    to_net <- to_net %>%
        dplyr::left_join(include_ranges, by = c("rname" = "contig")) %>%
        dplyr::filter(
            (pos >= sl & pos <= el) |
                (pos >= sr & pos <= er)
        ) %>%
        dplyr::select(-sl, -el, -sr, -er)

    to_net$tlen <- DATA$len[to_net$rname]

    to_net$flank <- ifelse(
        to_net$pos >= (to_net$tlen / 2),
        paste0(to_net$rname, "_r"),
        paste0(to_net$rname, "_l")
    )

    to_net$flank <- factor(to_net$flank)

    return(to_net)
}

#' Title
#'
#' @param to_net TODO
#'
#' @return TODO
#' @keywords internal
calculate_shared_matrix <- function(to_net) {
    flank_levels <- levels(to_net$flank)
    shared_matrix <- matrix(
        0,
        nrow = length(flank_levels),
        ncol = length(flank_levels),
        dimnames = list(flank_levels, flank_levels)
    )

    tri <- upper.tri(shared_matrix, diag = FALSE)
    shared_matrix <- reshape2::melt(shared_matrix)[tri, ]
    foo_net <- table(to_net$qname, to_net$flank)

    nc <- parallel::detectCores(logical = FALSE) - 1
    cl <- parallel::makeCluster(nc)
    parallel::clusterExport(cl, "foo_net", envir = environment())

    shared_matrix$value <- parallel::parApply(
        cl,
        shared_matrix,
        1,
        FUN = function(x) {
            col1 <- as.character(x[1])
            col2 <- as.character(x[2])
            if (!(col1 %in% colnames(foo_net)) || !(col2 %in% colnames(foo_net))) return(NA)

            present1 <- foo_net[, col1] != 0
            present2 <- foo_net[, col2] != 0
            union_size <- sum(present1 | present2)
            if (union_size == 0) return(0)
            intersection_size <- sum(present1 & present2)
            intersection_size / union_size
        }
    )
    parallel::stopCluster(cl)
    return(shared_matrix[shared_matrix$value != 0, ])
}

#' Title
#'
#' @param shared_matrix TODO
#' @param DATA TODO
#'
#' @return TODO
#' @keywords internal
detect_cluster <- function(shared_matrix, DATA) {
    left_right_pairs <- data.frame(
        Var1 = paste0(rownames(DATA),"_l"),
        Var2 = paste0(rownames(DATA),"_r"),
        value = 1
    )

    shared_matrix <- rbind(shared_matrix, left_right_pairs)

    filtered_edges <- shared_matrix[shared_matrix$value > 0.7, ]

    g <- igraph::graph_from_data_frame(filtered_edges, directed = TRUE)

    igraph::E(g)$weight <- filtered_edges$value

    com1 <- igraph::cluster_walktrap(g, weights = igraph::E(g)$weight)

    m1 <- igraph::membership(com1)

    contig_names <- paste0("contig_", sapply(strsplit(names(m1), "_"), `[`, 2))

    m1_table <- table(contig_names, m1)
    m1_melted <- reshape2::melt(m1_table)

    m1_filtered <- subset(m1_melted, value == 2)

    rownames(m1_filtered) <- m1_filtered$contig_names

    m1_filtered$m1 <- factor(m1_filtered$m1)

    clusters_with_multiple <- names(which(table(m1_filtered$m1) != 1))
    m1_final <- subset(m1_filtered, m1 %in% clusters_with_multiple)

    m1_final$m1 <- factor(m1_final$m1)
    return(m1_final)
}

#' Title
#'
#' @param m1_final TODO
#' @param DATA TODO
#'
#' @return TODO
#' @keywords internal
refine_taxonomic_assignments <- function(m1_final, DATA) {
    for (group in levels(m1_final$m1)) {
        group_contigs <- rownames(m1_final)[m1_final$m1 == group]
        new_tax <- choose_best_taxon(group_contigs, DATA)
        DATA$tax_cons[group_contigs] <- new_tax
    }
    return(DATA)
}

#' Title
#'
#' @param group_contigs TODO
#' @param DATA TODO
#'
#' @return TODO
#' @keywords internal
choose_best_taxon <- function(group_contigs, DATA) {
    cluster_data <- DATA[group_contigs, ]
    candidate_taxa <- unique(cluster_data$tax_cons)

    scores <- sapply(candidate_taxa, function(taxon) {
        members <- DATA[DATA$tax_cons == taxon, ]
        cov_diff <- abs(members$cov - mean(cluster_data$cov))
        mean_diff <- mean(cov_diff)
        size <- nrow(members)
        list(score = mean_diff, size = size)
    })

    best <- names(which.min(sapply(scores, `[[`, "score")))
    return(best)
}
