mycoBinR
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

# mycoBinR

<!-- badges: start -->

<!-- badges: end -->

**mycoBinR** helps you organize pipeline file paths and build tidy,
analysis-ready tables for fungal/metagenome assemblies. It integrates
contig metrics (coverage, GC, BUSCO, taxonomy), parses telomeres,
refines consensus taxonomy using read-connection “bubbles”, and exports
per-taxon **FASTA**/**BED** bins.

This README is written in the style of the
[*r-pkgs*](https://r-pkgs.org/) website: short motivation, reproducible
examples, and clear function sections.

------------------------------------------------------------------------

## Installation

You can install the development version of **mycoBinR** from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("taui11/mycoBinR")
```

------------------------------------------------------------------------

## Example

The `create_filepaths_df()` function builds a data frame listing file
names and absolute paths for a given project directory, project number,
and title. This *paths index* is used by the downstream functions.

``` r
library(mycoBinR)

filepaths <- create_filepaths_df(
  base_path = "/project/data_output",
  project_nr = "pr_01_",
  TITLE     = "001"
)
print(filepaths)
#>              name
#> 1        assembly
#> 2        coverage
#> 3           busco
#> 4        taxonomy
#> 5   assembly_info
#> 6    fasta_folder
#> 7      bed_folder
#> 8        bam_file
#> 9   telomere_left
#> 10 telomere_right
#>                                                                                                               path
#> 1                                                   /project/data_output/flye_output/pr_01_001_flye/assembly.fasta
#> 2                                          /project/data_output/coverm_output/pr_01_001_mean_cov/mean_coverage.tsv
#> 3                            /project/data_output/busco_output/pr_01_001_busco/run_ascomycota_odb12/full_table.tsv
#> 4                                  /project/data_output/diamond_output/pr_01_001_taxonomy_prots/001_taxonomy_prots
#> 5                                                /project/data_output/flye_output/pr_01_001_flye/assembly_info.txt
#> 6                                                     /project/data_output/binning_output/pr_01_001_binned/contigs
#> 7                                                        /project/data_output/binning_output/pr_01_001_binned/beds
#> 8  /project/data_output/coverm_output/pr_01_001_mean_cov/bam_cache/assembly.fasta.pr_01_001.hifireads.fastq.gz.bam
#> 9                                     /project/data_output/telfinder_output/pr_01_001_telomeres/NCR_left_score.txt
#> 10                                   /project/data_output/telfinder_output/pr_01_001_telomeres/NCR_right_score.txt
```

------------------------------------------------------------------------

## Function index

| Function | Purpose | Typical key inputs | Output |
|----|----|----|----|
| `create_filepaths_df()` | Declaratively define where pipeline outputs live. | `base_path`, `project_nr`, `TITLE` | tibble with `name`, `path` |
| `build_contig_table()` | Assemble a master contig table (coverage, GC, BUSCO, taxonomy, clustering). | `paths_df`, `api_key` (optional) | data frame, one row per contig |
| `parse_telomeres()` | Parse left/right telomere motifs; classify completeness. | `DATA`, `paths_df` | telomere annotations per contig |
| `incorporate_bam_bubbles()` | Use supplementary read connections to refine `tax_cons`. | `DATA`, `telomeres`, `paths_df` | updated contig table |
| `export_binned()` | Write per-taxon FASTA and BED files. | `DATA`, `paths_df`, `SAMPLE` | files on disk |

------------------------------------------------------------------------

## Typical workflow

``` r
library(mycoBinR)

# 1) Describe your pipeline outputs (only once per sample/run)
paths_df <- create_filepaths_df(
  base_path  = "/project/data_output",
  project_nr = "pr_01_",
  TITLE      = "001"
)

# 2) Build the contig-level table (coverage, GC, BUSCO, taxonomy, clustering)
DATA <- build_contig_table(paths_df, api_key = Sys.getenv("ENTREZ_KEY"))

# 3) Parse telomere motif calls and classify completeness
telos <- parse_telomeres(DATA, paths_df)

# 4) Refine consensus taxonomy using read-connection “bubbles”
DATA2 <- incorporate_bam_bubbles(DATA, telos, paths_df)

# 5) Export taxon-specific FASTA and BED bins
export_binned(DATA2, paths_df, SAMPLE = "001")
```

------------------------------------------------------------------------

## `build_contig_table()`

Constructs a tidy, per-contig data table by merging assembly statistics,
GC content, coverage values, BUSCO completeness, and taxonomy
assignments.  
This table is the central object used for downstream filtering, telomere
parsing, and binning.

Optionally, you can provide an NCBI **Entrez API key** to make taxonomy
lookups faster and avoid rate limits.

``` r
# Build the main contig-level table
DATA <- build_contig_table(
  paths_df,
  api_key = Sys.getenv("ENTREZ_KEY")  # optional but recommended
)

# View key summary columns
key_cols <- c(
  "contig",          # contig ID
  "length",          # contig length in bp
  "gc",              # GC content
  "cov_mean",        # mean coverage
  "busco_complete",  # BUSCO completeness fraction
  "tax_raw",         # raw taxonomy hit
  "tax_cons"         # consensus taxonomy
)

# Display only the relevant columns (if present)
print(DATA[intersect(names(DATA), key_cols)])
```

**Inputs**

- `paths_df`: from `create_filepaths_df()`, must point to assembly,
  coverage, BUSCO, taxonomy outputs, etc.  
- `api_key`: optional NCBI Entrez key for stable taxonomy lookups.

**Output**

- Data frame with per-contig metrics and taxonomic fields (e.g.,
  `tax_raw`, `tax_cons`).

------------------------------------------------------------------------

## `parse_telomeres()`

Parse left/right telomere motif hits and classify telomere completeness
for each contig.

``` r
telos <- parse_telomeres(
  DATA,
  paths_df
)

head(telos)
# dplyr::count(telos, telomere_class)
```

**Expected columns (typical)**  
`contig`, `leftmotif`, `leftscore`, `rightmotif`, `rightrc`,
`rightscore`, `telomere_class`.

------------------------------------------------------------------------

## `incorporate_bam_bubbles()`

Use supplementary alignments (from your BAM) to find connected contig
“bubbles” and propagate/refine `tax_cons` across those connections.

``` r
DATA2 <- incorporate_bam_bubbles(
  DATA,
  telomeres = telos,
  paths_df  = paths_df
)

# Inspect changes in consensus taxonomy
# dplyr::count(DATA$tax_cons)  |> dplyr::arrange(desc(n))
# dplyr::count(DATA2$tax_cons) |> dplyr::arrange(desc(n))
```

**Notes**

- Ensure your mapper emitted *supplementary alignments*; otherwise
  connectivity may be sparse.  
- `paths_df` must reference the correct BAM and any required indices.

------------------------------------------------------------------------

## `export_binned()`

Export one **FASTA** and one **BED** per **taxonomic bin**. Filenames
usually include your `SAMPLE` tag.

``` r
export_binned(
  DATA2,
  paths_df,
  SAMPLE = "001"  # appears in output filenames
)

# After running, inspect output directories referenced in paths_df
# (e.g., fasta_folder, bed_folder).
```

**Side effects**

- Writes files to disk; no return value. Designed for downstream
  visualization and analysis.

------------------------------------------------------------------------

## Troubleshooting

- **File not found**: print `paths_df` and confirm every expected `name`
  has a valid `path`.  
- **Taxonomy rate limits**: set `ENTREZ_KEY` in your shell profile (or
  within R) before running `build_contig_table()`.  
- **No effect from bubbles**: verify BAM has supplementary alignments;
  double-check `paths_df` entries for BAM/index.  
- **Reproducibility**: keep your directory layout stable so
  `create_filepaths_df()` can infer paths consistently.

------------------------------------------------------------------------

## Session info

``` r
sessionInfo()
#> R version 4.5.1 (2025-06-13)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Linux Mint 21.3
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so;  LAPACK version 3.10.0
#> 
#> locale:
#>  [1] LC_CTYPE=de_DE.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=de_DE.UTF-8    
#>  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=de_DE.UTF-8   
#>  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Berlin
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] mycoBinR_0.0.0.9000
#> 
#> loaded via a namespace (and not attached):
#>  [1] compiler_4.5.1    magrittr_2.0.3    fastmap_1.2.0     cli_3.6.5        
#>  [5] tools_4.5.1       htmltools_0.5.8.1 rstudioapi_0.17.1 yaml_2.3.10      
#>  [9] rmarkdown_2.29    knitr_1.50        xfun_0.52         digest_0.6.37    
#> [13] rlang_1.1.6       evaluate_1.0.4
```
