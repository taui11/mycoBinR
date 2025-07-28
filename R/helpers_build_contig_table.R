#' Title
#'
#' @param paths_df TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
load_assembly_info <- function(paths_df) {
    report <- read.delim(paths_df$path[paths_df$name == "assembly_info"], header = TRUE)
    rownames(report) <- report[, 1]
    return(report)
}

#' Title
#'
#' @param paths_df TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
load_dna_sequences <- function(paths_df) {
    read.dna(paths_df$path[paths_df$name == "assembly"], format="fasta")
}

#' Title
#'
#' @param paths_df TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
load_coverage_data <- function(paths_df) {
    read.delim(paths_df$path[paths_df$name == "coverage"], row.names = 1)
}

#' Title
#'
#' @param paths_df TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
load_busco_data <- function(paths_df) {
    read_busco_table(paths_df$path[paths_df$name == "busco"])
}

#' Title
#'
#' @param a TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
compute_gc_content <- function(a) {
    sapply(a, function(seq) {
        class(seq) <- "DNAbin"
        GC.content(seq)
    })
}

#' Title
#'
#' @param file TODO
#' @param bacterial TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
read_busco_table <- function(file, bacterial = FALSE) {
    df <- read_busco_file(file)
    df <- clean_busco_entries(df, bacterial = bacterial)
    summarize_busco_by_squence(df)
}

#' Title
#'
#' @param file TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
read_busco_file <- function(file) {
    df <- read.delim2(file, sep = "\t", header = TRUE, skip = 2)
    df[df$Status != "Missing", ]
}

#' Title
#'
#' @param df TODO
#' @param bacterial TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
clean_busco_entries <- function(df, bacterial = FALSE) {
    df$Sequence <- sapply(strsplit(as.character(df$Sequence), ":"), `[`, 1)
    if (bacterial) {
        parts <- strsplit(as.character(df[, 3]), "_")
        df$Sequence <- sapply(parts, function(x) paste(head(x, length(x) - 9), collapse = "_"))
    }
    df$Sequence <- factor(df$Sequence)
    df[, 2] <- factor(as.character(df[, 2]))
    df
}

#' Title
#'
#' @param df TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
summarize_busco_by_squence <- function(df) {
    table(df[, 3], df[, 2])
}

#' Title
#'
#' @param report TODO
#' @param nombres TODO
#' @param cov TODO
#' @param gc TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
build_initial_df <- function(report, nombres, cov, gc) {
    data.frame(
        cov=cov[nombres,1],
        len=report[nombres,2],
        gc=gc[nombres],
        circ=report[nombres,4],
        rep=report[nombres,5],
        multiplicity=report[nombres,6],
        alt_group=report[nombres,7],
        path=report[nombres,8],
        complete_graph=0,
        Complete=0,
        Duplicated=0,
        Fragmented=0
    )
}

#' Title
#'
#' @param data TODO
#' @param report TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
enrich_with_graph_info <- function(data, report) {
    graph_paths <- report[rownames(data), 8]
    star_pattern <- grepl("\\*", graph_paths)
    data$complete_graph <- 0
    data$complete_graph[star_pattern] <- 1
    data$complete_graph[star_pattern & duplicated(graph_paths)] <- 2
    return(data)
}

#' Title
#'
#' @param data TODO
#' @param buscos TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
add_busco_data <- function(data, buscos) {
    busco_cols <- intersect(colnames(buscos), colnames(data))
    data[rownames(buscos), busco_cols] <- buscos[, busco_cols]
    return(data)
}

#' Title
#'
#' @param data TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
cluster_data <- function(data) {
    clusters <- Mclust(data[, c("cov", "gc")], G = 1:20)
    data$clas <- clusters$classification
    return(data)
}

#' Title
#'
#' @param data TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
flag_busco_completeness <- function(data) {
    data$bc1 <- data$Complete != 0
    data$bc_combined <- factor(as.numeric(data$Complete != 0) + as.numeric(data$Duplicated != 0))
    return(data)
}

#' Title
#'
#' @param TAX TODO
#' @param API_KEY TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
process_taxonomy <- function(TAX, API_KEY = NULL) {
    tax <- load_taxonomy_table(TAX)
    tax <- parse_protein_ids(tax)
    tax <- filter_valid_taxids(tax)
    taxonomy_list <- fetch_tax_classification(tax$taxid, API_KEY)
    tax_df <- build_taxonomy_dataframe(taxonomy_list)
    tax <- merge_tax_data(tax, tax_df)
    tax_filled <- fill_missing_taxonomy(tax)
    cl_list <- build_taxonomic_partitions(tax_filled)
    return(cl_list)
}

#' Title
#'
#' @param TAX TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
load_taxonomy_table <- function(TAX) {
    tax <- read.delim(TAX, header = FALSE, stringsAsFactors = FALSE)
    colnames(tax)<-c("prot","taxid","ev","path")
    dup_idx <- duplicated(tax$prot)
    tax$prot[dup_idx] <- paste0(tax$prot[dup_idx], "_dup")
    return(tax)
}

#' Title
#'
#' @param tax TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
parse_protein_ids <- function(tax) {
    parsed <- strsplit(sapply(strsplit(tax$prot, "\\|"), `[`, 2), split = "[:\\-]")
    tax <- data.frame(
        contig = sapply(parsed, `[`, 1),
        start = sapply(parsed, `[`, 2),
        end = sapply(parsed, `[`, 3),
        taxid = tax$taxid,
        ev = tax$ev,
        stringsAsFactors = FALSE
    )
    return(tax)
}

#' Title
#'
#' @param tax TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
filter_valid_taxids <- function(tax) {
    tax <- tax[tax$taxid != 0, ]
}

#' Title
#'
#' @param taxids TODO
#' @param API_KEY TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
fetch_tax_classification <- function(taxids, API_KEY = NULL) {
    #if (!is.null(API_KEY)) Sys.setenv(ENTREZ_KEY = API_KEY)
    Sys.setenv(ENTREZ_KEY = API_KEY)
    classification(unique(taxids), db = "ncbi")
}

#' Title
#'
#' @param taxonomy_list TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
build_taxonomy_dataframe <- function(taxonomy_list) {
    class(taxonomy_list) <- NULL
    tibble(
        names = names(taxonomy_list),
        taxonomy_list = taxonomy_list
    ) %>%
        unnest(cols = c(taxonomy_list)) %>%
        filter(rank %in% c("domain", "kingdom", "phylum", "class", "order", "family", "genus")) %>%
        select(-id) %>%
        tidyr::pivot_wider(names_from = rank, values_from = name) %>%
        select(name = names, domain, kingdom, phylum, class, order, family, genus) %>%
        as.data.frame() %>%
        {
            rownames(.) <- make.unique(as.character(.$name))
            .
        }
}

#' Title
#'
#' @param tax TODO
#' @param tax_df TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
merge_tax_data <- function(tax, tax_df) {
    cbind(tax, tax_df[as.character(tax$taxid), ])
}

#' Title
#'
#' @param tax TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
fill_missing_taxonomy <- function(tax) {
    tax2 <- tax
    levels <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus")

    for (i in seq_along(levels)) {
        parent <- if (i == 1) levels[1] else levels[i - 1]
        current <- levels[i]

        missing_idx <- is.na(tax2[[current]]) & !is.na(tax2[[parent]])
        tax2[[current]][missing_idx] <- tax2[[parent]][missing_idx]
    }
    return(tax2)
}

#' Title
#'
#' @param tax TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
build_taxonomic_partitions <- function(tax) {
    levels <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus")
    lapply(levels, function(level){
        taxonomy_list <- table(tax$contig, tax[[level]])
        taxonomy_list <- taxonomy_list / rowSums(taxonomy_list)
        taxonomy_list[is.na(taxonomy_list)] <- 0
        taxonomy_list <- cbind(taxonomy_list, unknown = 0)
        taxonomy_list[rowSums(taxonomy_list) == 0, "unknown"] <- 1
        as.cl_partition(taxonomy_list)
    })
}

#' Title
#'
#' @param cl_list TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
compute_consensus_clusters <- function(cl_list) {
    cl_ens <- cl_ensemble(list = cl_list)
    consensus <- cl_consensus(cl_ens, method="SE", control = list(verbose = TRUE, nruns = 10, maxiter = 100))
    cl_class_ids(consensus)
}

#' Title
#'
#' @param data TODO
#' @param consensus_ids TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
assign_consensus_clusters <- function(data, consensus_ids) {
    data$tax_cons <- 0
    data[names(consensus_ids), "tax_cons"] <- consensus_ids
    return(data)
}


#' Title
#'
#' @param data TODO
#' @param cl_list TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
extract_taxonomic_labels <- function(data, cl_list) {
    levels_tax <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus")
    classification <- matrix(
        "unknown",
        nrow = nrow(data),
        ncol = length(levels_tax),
        dimnames = list(rownames(data), levels_tax)
    )
    for (i in seq_along(levels_tax)) {
        cl <- cl_list[[i]]
        cluster_ids <- cl_class_ids(cl)
        rownames_matrix <- rownames(cl[[1]])
        colnames_matrix <- colnames(cl[[1]])
        assigned <- colnames_matrix[cluster_ids]
        classification[rownames_matrix, i] <- assigned
    }

    classification_df <- as.data.frame(classification, stringsAsFactors = FALSE)
    classification_df <- classification_df[rownames(data), , drop = FALSE]
    return(classification_df)
}

#' Title
#'
#' @param data TODO
#'
#' @returns TODO
#' @keywords  internal
#'
#' @examples TODO
summarize_taxonomy <- function(data) {
    data$summ_class <- ifelse(data$domain == "Bacteria", "Bacteria", data$class)
    data$summ_class[data$class == "Cyanophyceae"] <- "Cyanophyceae"

    data$tax_cons <- factor(data$tax_cons)
    recoded <- table(data$tax_cons, data$summ_class)

    new_levels <- sapply(levels(data$tax_cons), function(x) {
        top_class <- names(which.max(recoded[x, ]))
        paste(x, top_class, sep="_")
    })

    new_levels <- gsub("^\\d+_unknown$", "unknown", new_levels)
    levels(data$tax_cons)  <- new_levels

    return(data)
}
