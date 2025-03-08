#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: run_GSEA.R
# Description:
#   This script performs Gene Set Enrichment Analysis (GSEA) on differential 
#   expression (DEG) results. The DEG data is provided as a CSV file, and the 
#   gene sets are retrieved using msigdbr based on a specified category and 
#   subcategory. The analysis is performed using clusterProfiler.
#
# Usage:
#   Rscript run_GSEA.R --deg-path <path_to_DEG_csv> --category <CATEGORY_SUBCATEGORY>
#         --suffix <output_suffix> --out <output_folder>
#
#   Example:
#     Rscript run_GSEA.R --deg-path data/deg_results.csv --category HALLMARK_Null \
#         --suffix v1 --out results/
#
# Output:
#   A CSV file with the GSEA results is saved in the specified output folder with a
#   filename constructed from the provided category and suffix.
# ------------------------------------------------------------------------------

# Load required libraries with suppressed startup messages -------------------
suppressPackageStartupMessages({
  library(here)
  set_here(path = normalizePath('..'))
  # AnnotationHub is set to use the local cache
  library(AnnotationHub)
  AnnotationHub::setAnnotationHubOption("LOCAL", TRUE)
  
  library(optparse)
  library(dplyr)
  library(stringr)
  library(msigdbr)
  library(clusterProfiler)
})

#' Run Gene Set Enrichment Analysis (GSEA)
#'
#' This function performs GSEA on DEG data. It:
#'   1. Filters out genes with missing log2 fold change values.
#'   2. Creates a named vector of log2 fold changes, sorted in descending order.
#'   3. If the subcategory is provided as "Null", it resets it to NULL to include
#'      all gene sets in the category.
#'   4. Retrieves gene set definitions from msigdbr for the specified category and
#'      subcategory.
#'   5. Runs GSEA using clusterProfiler and returns the result table.
#'
#' @param DEG Data frame. Differential expression results with rownames as gene
#'   identifiers and a column "log2FoldChange".
#' @param category Character. The main gene set category (e.g., "HALLMARK").
#' @param subcategory Character or NULL. The gene set subcategory. If "Null" is passed,
#'   it is treated as NULL.
#'
#' @return A data frame containing the GSEA result table.
#'
#' @examples
#' \dontrun{
#'   deg <- read.csv("data/deg_results.csv", row.names = 1)
#'   gsea_results <- run_GSEA(deg, "HALLMARK", NULL)
#' }
run_GSEA <- function(DEG, category, subcategory) {
  # Filter out rows with missing log2FoldChange
  DEG <- dplyr::filter(DEG, !is.na(log2FoldChange))
  
  # Create a ranked gene list using log2FoldChange values
  rankd_genes <- DEG$log2FoldChange
  names(rankd_genes) <- rownames(DEG)
  rankd_genes <- sort(rankd_genes, decreasing = TRUE)
  
  # If subcategory is 'Null', treat it as NULL
  if (subcategory == 'Null') {
    subcategory <- NULL
  }
  
  # Retrieve gene set definitions using msigdbr for the specified category/subcategory
  term_to_gene <- msigdbr(species = "Homo sapiens", category = category, subcategory = subcategory) %>%
    dplyr::select(gs_name, ensembl_gene)
  
  # Run GSEA using clusterProfiler
  gsea_res <- clusterProfiler::GSEA(geneList = rankd_genes, TERM2GENE = term_to_gene, pvalueCutoff = 1)
  
  # Return the GSEA result table
  return(gsea_res@result)
}

# Define command-line options -----------------------------------------------
option_list <- list(
  make_option(c("--deg-path"),
              type = "character",
              default = NULL,
              help = "Path to the DEG CSV file.",
              metavar = "character"),
  make_option(c("--category"),
              type = "character",
              default = NULL,
              help = "Gene set category and subcategory separated by an underscore (e.g., 'HALLMARK_Null').",
              metavar = "character"),
  make_option(c("--suffix"),
              type = "character",
              default = NULL,
              help = "Suffix to append to the output filename.",
              metavar = "character"),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "Output folder path.",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that all mandatory options are provided -----------------------------
if (any(unlist(lapply(option_list, function(option) { is.null(opt[[option@dest]]) })))) {
  print_help(opt_parser)
  stop("Error: All options are mandatory", call. = FALSE)
}

# Split the category parameter into category and subcategory components
category <- str_split_1(opt$category, '_')[1]
subcategory <- str_split_1(opt$category, '_')[2]
suffix <- opt$suffix

# Read DEG results from the provided CSV file
res <- tryCatch(
  expr = {
    deg <- read.csv(opt$`deg-path`, row.names = 1)
    run_GSEA(deg, category, subcategory)
  },
  error = function(e) {
    message("Caught an error!")
    print(e)
    c()
  }
)

# Create the output directory if it doesn't exist
if (!dir.exists(opt$out)) {
  dir.create(opt$out, recursive = TRUE)
}

# Write the GSEA results to a CSV file
output_file <- here::here(opt$out, paste0(opt$category, "_GSEA_", suffix, ".csv"))
write.csv(res, output_file, row.names = TRUE)

message("GSEA analysis complete. Results written to: ", output_file)
