if (getRversion() >= "2.15.1") {
    utils::globalVariables(c(".data", "id", "name", "domain", "kingdom", "phylum", "family", "genus", "rank", "names", "taxonomy_list"))
}

#' Loading assembly_info.txt file
#'
#' @param paths_df Data frame. Typically created with \code{create_filepaths_df}, containing columns \code{name} and \code{path} for locating input files.
#'
#' @return A data frame with the assembly info and the rownames as the contig names.
#' @keywords  internal
load_assembly_info <- function(paths_df) {
    report <- utils::read.delim(paths_df$path[paths_df$name == "assembly_info"], header = TRUE)
    rownames(report) <- report[, 1]
    return(report)
}

#' Loading assembly.fasta file
#'
#' @param paths_df Data frame. Typically created with \code{create_filepaths_df}, containing columns \code{name} and \code{path} for locating input files.
#'
#' @return A matrix of class \code{"DNAbin"} containing DNA sequences
#' @keywords  internal
load_dna_sequences <- function(paths_df) {
    ape::read.dna(paths_df$path[paths_df$name == "assembly"], format="fasta")
}

#' Loading mean_coverage.tsv file
#'
#' @param paths_df Data frame. Typically created with \code{create_filepaths_df}, containing columns \code{name} and \code{path} for locating input files.
#'
#' @return A data frame with the coverage data
#' @keywords  internal
load_coverage_data <- function(paths_df) {
    utils::read.delim(paths_df$path[paths_df$name == "coverage"], row.names = 1)
}

#' Loading full_table.tsv file from busco_output
#'
#' @param paths_df Data frame. Typically created with \code{create_filepaths_df}, containing columns \code{name} and \code{path} for locating input files.
#'
#' @return A data frame with BUSCO data
#' @keywords  internal
load_busco_data <- function(paths_df) {
    read_busco_table(paths_df$path[paths_df$name == "busco"])
}

#' Computing GC content
#'
#' @param dna_seq A matrix of class \code{"DNAbin"}. Containing DNA sequences to be analyzed.
#'
#' @return a single numeric value representing the GC content of that sequence
#' @keywords  internal
compute_gc_content <- function(dna_seq) {
    sapply(dna_seq, function(seq) {
        class(seq) <- "DNAbin"
        ape::GC.content(seq)
    })
}

#' Read and process BUSCO summary table
#'
#' @param file Character string. Path to the BUSCO summary file.
#' @param bacterial Logical. Whether the dataset is bacterial, affecting sequence name cleaning. Default is FALSE.
#'
#' @return A contingency table summarizing BUSCO counts by sequence and status.
#' @keywords  internal
read_busco_table <- function(file, bacterial = FALSE) {
    df <- read_busco_file(file)
    df <- clean_busco_entries(df, bacterial = bacterial)
    summarize_busco_by_squence(df)
}

#' Read raw BUSCO data file
#'
#' @param file Character string. Path to the BUSCO data file
#'
#' @return A filtered data frame of BUSCO entries excluding those with status "Missing".
#' @keywords  internal
read_busco_file <- function(file) {
    df <- utils::read.delim2(file, sep = "\t", header = TRUE, skip = 2)
    df[df$Status != "Missing", ]
}

#' Clean BUSCO sequence entries
#'
#' @param df Data frame. BUSCO data frame with at least a \code{Sequence} column.
#' @param bacterial Logical. Whether to apply bacterial-specific sequence name cleaning. Default is FALSE.
#'
#' @return A cleaned data frame with formatted \code{Sequence} and \code{Status} factors.
#' @keywords  internal
clean_busco_entries <- function(df, bacterial = FALSE) {
    df$Sequence <- sapply(strsplit(as.character(df$Sequence), ":"), `[`, 1)
    if (bacterial) {
        parts <- strsplit(as.character(df[, 3]), "_")
        df$Sequence <- sapply(parts, function(x) paste(utils::head(x, length(x) - 9), collapse = "_"))
    }
    df$Sequence <- factor(df$Sequence)
    df[, 2] <- factor(as.character(df[, 2]))
    df
}

#' Summarize BUSCO counts by sequence and status
#'
#' @param df Data frame. BUSCO data frame with sequence and status columns
#'
#' @return A contingency table (matrix) of counts indexed by sequence and status.
#' @keywords  internal
summarize_busco_by_squence <- function(df) {
    table(df[, 3], df[, 2])
}

#' Build initial contig-level data frame
#'
#' @param report Data frame. Contains the assembly info and the rownames as the contig names.
#' @param contig_names Character vector. Contains the contig IDs.
#' @param cov Data frame. Contains the covetage data.
#' @param gc Named numeric vector. Contains GC content values, names should match contig IDs.
#'
#' @return A Data frame with contig-level metrics and placeholder columns for later annotation.
#' @keywords  internal
build_initial_df <- function(report, contig_names, cov, gc) {
    data.frame(
        cov=cov[contig_names, 1],
        len=report[contig_names, 2],
        gc=gc[contig_names],
        circ=report[contig_names, 4],
        rep=report[contig_names, 5],
        multiplicity=report[contig_names, 6],
        alt_group=report[contig_names, 7],
        path=report[contig_names, 8],
        complete_graph = 0,
        Complete = 0,
        Duplicated = 0,
        Fragmented = 0
    )
}

#' Annotate contigs with graph completeness
#'
#' @param data Data frame. Contains a contig-level data frame from build_initial_df.
#' @param report Data frame. Contains the assembly info and the rownames as the contig names.
#'
#' @return A modified version of `data` with a `complete_graph` column where:
#'      \describe{
#'          \item{0}{Path does not contain a `"*"`}
#'          \item{1}{Path contains a `"*"`}
#'          \item{2}{Path contains a `"*"` and is duplicated}
#'      }
#' @keywords  internal
enrich_with_graph_info <- function(data, report) {
    graph_paths <- report[rownames(data), 8]
    star_pattern <- grepl("\\*", graph_paths)
    data$complete_graph <- 0
    data$complete_graph[star_pattern] <- 1
    data$complete_graph[star_pattern & duplicated(graph_paths)] <- 2
    return(data)
}

#' Adding the BUSCO data to data frame
#'
#' @param data Data frame. Contains a contig-level data frame from enrich_with_graph_info.
#' @param buscos Data frame Contains BUSCO data.
#'
#' @return A data frame with additional columns containing BUSCO metrics.
#' @keywords  internal
add_busco_data <- function(data, buscos) {
    busco_cols <- intersect(colnames(buscos), colnames(data))
    data[rownames(buscos), busco_cols] <- buscos[, busco_cols]
    return(data)
}

#' Creating and adding cluster classification based on coverage and GC content.
#'
#' @param data Data frame. Contains a contig-level data frame with columns `cov`and `gc`.
#'
#' @return A data frame with added `clas` column indicating cluster membership.
#' @keywords  internal
cluster_data <- function(data) {
    clusters <- mclust::Mclust(data[, c("cov", "gc")], G = 1:20)
    data$clas <- clusters$classification
    return(data)
}

#' Flag BUSCO completeness and duplication status
#'
#' @param data Data frame. Should contain at least `Complete` and a combined completeness-duplication factor.
#'
#' @return The input data frame with two additional columns:
#'      \describe{
#'          \item{bc1}{Logical vector indicating if `Complete` is non-zero (TRUE if complete).}
#'          \item{bc_combined}{Factor indicating combined BUSCO completeness and duplication status:
#'              \itemize{
#'                  \item 0 = neither complete nor duplicated
#'                  \item 1 = complete only
#'                  \item 2 = complete and duplicated
#'              }
#'          }
#'      }
#' @keywords  internal
flag_busco_completeness <- function(data) {
    data$bc1 <- data$Complete != 0
    data$bc_combined <- factor(as.numeric(data$Complete != 0) + as.numeric(data$Duplicated != 0))
    return(data)
}

#' Process taxonomy data through multiple steps
#'
#' @param tax_path Character string. Contains the path to the taxonomy table.
#' @param api_key Character string or NULL. NCBI API key needed for fetching taxonomy classification.
#'
#' @return A list of taxonomic partitions constructed from processed taxonomy data.
#' @keywords  internal
process_taxonomy <- function(tax_path, api_key = NULL) {
    tax <- load_taxonomy_table(tax_path)
    tax <- parse_protein_ids(tax)
    tax <- filter_valid_taxids(tax)
    taxonomy_list <- fetch_tax_classification(tax$taxid, api_key)
    tax_df <- build_taxonomy_dataframe(taxonomy_list)
    tax <- merge_tax_data(tax, tax_df)
    tax_filled <- fill_missing_taxonomy(tax)
    cl_list <- build_taxonomic_partitions(tax_filled)
    return(cl_list)
}

#' Load taxonomy table from a file
#'
#' @param tax_path Character string. File path to the taxonomy table
#'
#' @return A data frame with columns `prot`, `taxid`, `ev`, and `path`.
#'      Duplicated protein IDs are suffixed with `_dup`
#' @keywords  internal
load_taxonomy_table <- function(tax_path) {
    tax <- utils::read.delim(tax_path, header = FALSE, stringsAsFactors = FALSE)
    colnames(tax)<-c("prot","taxid","ev","path")
    dup_idx <- duplicated(tax$prot)
    tax$prot[dup_idx] <- paste0(tax$prot[dup_idx], "_dup")
    return(tax)
}

#' Parse protein identifiers into config, start, end, and taxid columns
#'
#' @param tax Data frame. Must contain a column `prot` with protein IDs and a `taxid` column. Typically from load_taxonomy_table
#'
#' @return A Data frame with columns `contig`, `start`, `end` and  `ev`
#' @keywords  internal
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

#' Filter taxonomy entries with valid taxonomic IDs
#'
#' @param tax Data frame. Must contain a `taxid` column. Typically from parse_protein_ids.
#'
#' @return Filtered data frame with only rows where `taxid` is not zero.
#' @keywords  internal
filter_valid_taxids <- function(tax) {
    tax <- tax[tax$taxid != 0, ]
}

#' Fetch taxonomy classification from NCBI database
#'
#' @param taxids Integer or character vector of taxonomic IDs.
#' @param api_key Character string or NULL. API key needed for NCBI Entrez.
#'
#' @return A named list of taxonomic classification.
#' @keywords  internal
fetch_tax_classification <- function(taxids, api_key = NULL) {
    #if (!is.null(api_key)) Sys.setenv(ENTREZ_KEY = api_key)
    Sys.setenv(ENTREZ_KEY = api_key)
    taxize::classification(unique(taxids), db = "ncbi")
}

#' Build a taxonomy data frame from classification list
#'
#' @param taxonomy_list Named list. Taxonomic classification data frame filtered by major ranks.
#'
#' @importFrom magrittr %>%
#'
#' @return Data frame with columns from taxon names and major taxonomic ranks (domain to genus).
#' @keywords  internal
build_taxonomy_dataframe <- function(taxonomy_list) {
    class(taxonomy_list) <- NULL
    tibble::tibble(
        names = names(taxonomy_list),
        taxonomy_list = taxonomy_list
    ) %>%
        tidyr::unnest(cols = c(taxonomy_list)) %>%
        dplyr::filter(rank %in% c("domain", "kingdom", "phylum", "class", "order", "family", "genus")) %>%
        dplyr::select(-id) %>%
        tidyr::pivot_wider(names_from = rank, values_from = name) %>%
        dplyr::select(name = .data$names, .data$domain, .data$kingdom, .data$phylum, .data$class, .data$order, .data$family, .data$genus) %>%
        as.data.frame() %>%
        {
            rownames(.) <- make.unique(as.character(.$name))
            .
        }
}

#' Merge original taxonomy data with classification data frame
#'
#' @param tax Data frame. Original taxonomy data including `taxid`.
#' @param tax_df Data frame. Classification data frame keyed by taxonomic ID.
#'
#' @return Combined data frame with columns from both inputs.
#' @keywords  internal
merge_tax_data <- function(tax, tax_df) {
    cbind(tax, tax_df[as.character(tax$taxid), ])
}

#' Fill missing taxonomy ranks by propagating from higher levels
#'
#' @param tax Data frame. Taxonomy data frame with columns for ranks: domain, kingdom, phylum, etc.
#'
#' @return Taxonomy data frame with missing ranks filled where possible.
#' @keywords  internal
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

#' Build taxonomic partitions for clustering
#'
#' @param tax Data frame. Taxonomy data frame with contig and rank columns.
#'
#' @return List of `cl_partition` objects, one per taxonomic rank.
#' @keywords  internal
build_taxonomic_partitions <- function(tax) {
    levels <- c("domain", "kingdom", "phylum", "class", "order", "family", "genus")
    lapply(levels, function(level){
        taxonomy_list <- table(tax$contig, tax[[level]])
        taxonomy_list <- taxonomy_list / rowSums(taxonomy_list)
        taxonomy_list[is.na(taxonomy_list)] <- 0
        taxonomy_list <- cbind(taxonomy_list, unknown = 0)
        taxonomy_list[rowSums(taxonomy_list) == 0, "unknown"] <- 1
        clue::as.cl_partition(taxonomy_list)
    })
}

#' Compute consensus cluster classification from cluster list
#'
#' @param cl_list List of cluster partitions.
#'
#' @return Vector of factor of consensus cluster class IDs.
#' @keywords  internal
compute_consensus_clusters <- function(cl_list) {
    cl_ens <- clue::cl_ensemble(list = cl_list)
    consensus <- clue::cl_consensus(cl_ens, method="SE", control = list(verbose = TRUE, nruns = 10, maxiter = 100))
    clue::cl_class_ids(consensus)
}

#' Assign consensus cluster IDs to the data
#'
#' @param data Data frame. Input data to which consensus clusters will be assigned.
#' @param consensus_ids Named vector or factor. Consensus cluster IDs indexed by row names of `data`
#'
#' @return The input data frame with an additional column `tax_cons` containing cluster assignments.
#' @keywords  internal
assign_consensus_clusters <- function(data, consensus_ids) {
    data$tax_cons <- 0
    data[names(consensus_ids), "tax_cons"] <- consensus_ids
    return(data)
}


#' Extract taxonomic labels from cluster partitions
#'
#' @param data Data frame. Input data with rownames matching cluster partitions.
#' @param cl_list List. A list of cluster partition objects, each representing taxonomic clustering at different ranks.
#'
#' @return A data frame with taxonomic classification labels for each input row, with columns for taxonomic ranks.
#' @keywords  internal
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
        cluster_ids <- clue::cl_class_ids(cl)
        rownames_matrix <- rownames(cl[[1]])
        colnames_matrix <- colnames(cl[[1]])
        assigned <- colnames_matrix[cluster_ids]
        classification[rownames_matrix, i] <- assigned
    }

    classification_df <- as.data.frame(classification, stringsAsFactors = FALSE)
    classification_df <- classification_df[rownames(data), , drop = FALSE]
    return(classification_df)
}

#' Summarize taxonomic classification with simplified labels
#'
#' @param data Data frame. Input data with taxonomic columns (`domain`, `class` and `tax_cons`)
#'
#' @return The input data frame with updated `summ_class` and recoded `tax_cons` factor levels.
#' @keywords  internal
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




