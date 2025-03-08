#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: controlled_GSEA.R
# Description:
#   This script performs a controlled Gene Set Enrichment Analysis (GSEA) on a 
#   differential expression dataset. It uses the msigdbr package to obtain gene 
#   set definitions and clusterProfiler for running GSEA. The analysis is performed 
#   on a list of genes with their corresponding effect sizes (e.g. from an IPTW analysis),
#   and the results are saved as a CSV file.
#
# Usage:
#   Rscript controlled_GSEA.R --arm <chromosome_arm> --deg-path <path_to_DEG_csv> 
#         --category <category_subcategory> --out <output_folder>
#
#   The DEG CSV file should contain at least a column named "effect.size" and "name".
#   The "category" parameter should be provided in the format "CATEGORY_SUBCATEGORY"
#   (e.g., "HALLMARK_Null" if no subcategory is desired).
#
# Output:
#   A CSV file containing the GSEA results is written to the specified output folder.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  set_here(path = normalizePath('..'))
  # Set AnnotationHub to use local cache
  library(AnnotationHub)
  AnnotationHub::setAnnotationHubOption("LOCAL", TRUE)
  
  # Load required libraries
  library(optparse)
  library(dplyr)
  library(stringr)
  library(msigdbr)
  library(clusterProfiler)
})

#' Get Controlled GSEA Results
#'
#' This function performs a Gene Set Enrichment Analysis (GSEA) using controlled 
#' differential expression data. It:
#'   1. Filters and orders the gene effect data (arm.del.effect) based on effect sizes.
#'   2. Constructs a named vector of effect sizes (gene list) using gene names.
#'   3. If the subcategory is 'Null', it sets it to NULL to retrieve all gene sets in the category.
#'   4. Retrieves gene set definitions using msigdbr based on the provided category and subcategory.
#'   5. Runs GSEA using clusterProfiler on the gene list.
#'
#' @param arm Character. The chromosome arm (e.g., "X3p"). (Not used in the current function.)
#' @param arm.del.effect Data frame. Differential expression results with at least columns
#'        "effect.size" and "name".
#' @param category Character. Gene set category (e.g., "HALLMARK").
#' @param subcategory Character or NULL. Gene set subcategory. If 'Null' is provided, it will be set to NULL.
#'
#' @return A data frame containing the GSEA result table.
#'
#' @examples
#' \dontrun{
#'   deg <- read.csv("deg_results.csv")
#'   gsea_results <- get_controlled_GSEA("X3p", deg, "HALLMARK", NULL)
#' }
get_controlled_GSEA <- function(arm, arm.del.effect, category, subcategory) {
  # Filter out genes with missing effect sizes and sort descending
  arm.del.effect <- arm.del.effect %>% 
    dplyr::filter(!is.na(effect.size)) %>% 
    arrange(desc(effect.size))
  
  # Create a named numeric vector of effect sizes
  gene.list <- arm.del.effect$effect.size
  names(gene.list) <- arm.del.effect$name
  
  # If subcategory is "Null", set it to NULL
  if (subcategory == 'Null') {
    subcategory <- NULL
  }
  
  # Retrieve gene set definitions from msigdbr for the specified category and subcategory
  term_to_gene <- msigdbr(species = "Homo sapiens", category = category, subcategory = subcategory) %>%
    dplyr::select(gs_name, gene_symbol)
  
  # Run GSEA using clusterProfiler
  gsea_res <- clusterProfiler::GSEA(gene.list, TERM2GENE = term_to_gene, pvalueCutoff = 1)
  
  return(gsea_res@result)
}

# Define command-line options -----------------------------------------------
option_list <- list(
  make_option(c("-a", "--arm"),
              type = "character",
              default = NULL,
              help = "Chromosome arm name (e.g., 'X3p')",
              metavar = "character"),
  make_option(c("--deg-path"),
              type = "character",
              default = NULL,
              help = "Path to the CSV file containing differential expression (DEG) results",
              metavar = "character"),
  make_option(c("--category"),
              type = "character",
              default = NULL,
              help = "Gene set category and subcategory separated by an underscore (e.g., 'HALLMARK_Null')",
              metavar = "character"),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "Output folder path",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that all mandatory options are provided -----------------------------
if (any(unlist(lapply(option_list, function(option) { is.null(opt[[option@dest]]) })))) {
  print_help(opt_parser)
  stop("Error: All options are mandatory", call. = FALSE)
}

# Split the category parameter into category and subcategory
category <- str_split_1(opt$category, '_')[1]
subcategory <- str_split_1(opt$category, '_')[2]

# Read DEG results from the provided CSV file
deg_results <- read.csv(opt$`deg-path`)

# Compute the controlled GSEA results using the provided parameters
res <- get_controlled_GSEA(opt$arm, deg_results, category, subcategory)

# Create the output directory if it doesn't exist
if (!dir.exists(opt$out)) {
  dir.create(opt$out, recursive = TRUE)
}

# Write the GSEA results to CSV
output_file <- here::here(opt$out, paste0(opt$category, "_IPTW_controlled_GSEA.csv"))
write.csv(res, output_file, row.names = FALSE)

message("Controlled GSEA analysis complete. Results written to: ", output_file)
