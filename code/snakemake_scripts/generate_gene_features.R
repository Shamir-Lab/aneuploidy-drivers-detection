#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: generate_gene_features.R
# Description:
#   Merges various gene-level data sets—including mutation data, CNV-based
#   probabilities, Fisher test results, and co-deletion counts—for a given cancer
#   type and chromosome arm. The final merged gene features table is output as a CSV.
#
# Usage:
#   Rscript generate_gene_features.R --cancer-type <cancer_type> --arm <arm>
#       --mutation-data-path <path_to_mutation_csv>
#       --focal-given-broad-probabilities-path <path_to_csv>
#       --focal-given-no-broad-probabilities-path <path_to_csv>
#       --co-deletion-fisher-path <path_to_csv>
#       --focal-and-broad-path <path_to_csv>
#       --out <output_folder>
#
# Example:
#   Rscript generate_gene_features.R --cancer-type BRCA --arm X3p 
#       --mutation-data-path data/mutations.csv 
#       --focal-given-broad-probabilities-path data/focal_broad.csv 
#       --focal-given-no-broad-probabilities-path data/focal_no_broad.csv 
#       --co-deletion-fisher-path data/co_del_fisher.csv 
#       --focal-and-broad-path data/focal_and_broad.csv 
#       --out results/
# ------------------------------------------------------------------------------

# Load required libraries and source utility functions ----------------------
suppressPackageStartupMessages({
  library(here)
  set_here(path=normalizePath('..'))
  source(here::here('genome_annotation_utils.R'))
  
  library(optparse)
  library(purrr)
  library(tidyr)
  library(plyr)
  library(dplyr)
})

#' Generate Gene Features
#'
#' This function merges various gene-level datasets for a given cancer type and 
#' chromosome arm. Data sources include:
#'   - Basic gene annotation data from protein-coding genes on the specified arm.
#'   - Mutation data (with Fisher test p-values and odds ratios).
#'   - Focal CNV probabilities given broad events.
#'   - Focal CNV probabilities given absence of broad events.
#'   - Co-deletion Fisher test q-values.
#'   - Counts of genes with both focal and broad events.
#'
#' @param cancer.type Character. Cancer type identifier (e.g., "BRCA").
#' @param arm Character. Chromosome arm name (e.g., "X3p").
#' @param mutation.data Data frame. Mutation data with at least columns "pval" and "odds.ratio".
#' @param focal.given.broad.probabilities Data frame. Row-named probabilities for focal events 
#'   given broad CNVs (expected column name: "prob").
#' @param focal.given.no.broad.probabilities Data frame. Row-named probabilities for focal events 
#'   given absence of broad CNVs (expected column name: "prob").
#' @param co.deletion.fisher Data frame. Contains Fisher test results with columns "name" (gene)
#'   and "qval" (adjusted p-value).
#' @param focal.and.broad Data frame. Row-named counts of genes with both focal and broad events 
#'   (expected column name: "count").
#'
#' @return A data frame of merged gene features with one row per gene. Missing values are 
#'         replaced with 0.
#'
#' @examples
#' \dontrun{
#'   features <- generate_gene_features("BRCA", "X3p", mutation.df, 
#'               focal.broad, focal.no.broad, fisher.df, focal.and.broad)
#' }
generate_gene_features <- function(cancer.type, arm, mutation.data,
                                   focal.given.broad.probabilities,
                                   focal.given.no.broad.probabilities,
                                   co.deletion.fisher, focal.and.broad) {
  # Retrieve protein-coding genes for the specified arm
  arm.genes <- dplyr::filter(get_arm_genes(arm, GRCh.version = "37"),
                             transcript_biotype == "protein_coding")
  
  # Merge gene-level features from various sources
  parameter.by.gene <- list(
    # Basic gene annotation data (ensure unique gene names and IDs)
    basic.gana.data = arm.genes %>% 
      filter(!is.na(entrezgene_id)) %>% 
      filter(!duplicated(external_gene_name)) %>% 
      filter(!duplicated(ensembl_gene_id)) %>% 
      dplyr::select(name = external_gene_name, gene = ensembl_gene_id, start_position, end_position),
    
    # Mutation data: rename columns for clarity (note: "mutatuin.odds.ratio" corrected to "mutation.odds.ratio")
    mutation.data %>% dplyr::select(-name) %>% 
      dplyr::rename(mutation.fisher.pval = pval, mutation.odds.ratio = odds.ratio),
    
    # Focal CNV probabilities given broad events
    data.frame(gene = rownames(focal.given.broad.probabilities), 
               focal.given.broad.probabilities = focal.given.broad.probabilities$prob),
    
    # Focal CNV probabilities given absence of broad events
    data.frame(gene = rownames(focal.given.no.broad.probabilities), 
               focal.given.no.broad.probabilities = focal.given.no.broad.probabilities$prob),
    
    # Co-deletion Fisher test q-values
    dplyr::select(co.deletion.fisher, gene = name, co.deletion.fisher.qval = qval),
    
    # Count of genes with both focal and broad events
    data.frame(gene = rownames(focal.and.broad), 
               focal.and.broad = focal.and.broad$count)
  )
  
  # Merge all feature data frames by the 'gene' column
  all.gene.parms <- Reduce(function(x, y) merge(x, y, all.x = TRUE, by = "gene"), parameter.by.gene)
  
  # Replace NA values with 0
  all.gene.parms[is.na(all.gene.parms)] <- 0
  
  return(all.gene.parms)
}

# Define command-line options -----------------------------------------------
option_list <- list(
  make_option(c("-c", "--cancer-type"),
              type = "character",
              default = NULL,
              help = "Cancer type (e.g., BRCA)",
              metavar = "character"),
  make_option(c("-a", "--arm"),
              type = "character",
              default = NULL,
              help = "Chromosome arm name (e.g., 'X3p')",
              metavar = "character"),
  make_option(c("--mutation-data-path"),
              type = "character",
              default = NULL,
              help = "Path to CSV file containing mutation data",
              metavar = "character"),
  make_option(c("--focal-given-broad-probabilities-path"),
              type = "character",
              default = NULL,
              help = "Path to CSV file with focal given broad probabilities",
              metavar = "character"),
  make_option(c("--focal-given-no-broad-probabilities-path"),
              type = "character",
              default = NULL,
              help = "Path to CSV file with focal given no broad probabilities",
              metavar = "character"),
  make_option(c("--co-deletion-fisher-path"),
              type = "character",
              default = NULL,
              help = "Path to CSV file with co-deletion Fisher test results",
              metavar = "character"),
  make_option(c("--focal-and-broad-path"),
              type = "character",
              default = NULL,
              help = "Path to CSV file with focal and broad event counts",
              metavar = "character"),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "Output folder path",
              metavar = "character")
)

# Parse command-line options -------------------------------------------------
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate that all required options are provided ---------------------------
if (any(unlist(lapply(option_list, function(option) { is.null(opt[[option@dest]]) })))) {
  print_help(opt_parser)
  stop("Error: All options are mandatory", call. = FALSE)
}

# Read input CSV files and generate gene features --------------------------
mutation.data <- read.csv(opt$`mutation-data-path`)
focal.given.broad.probabilities <- read.csv(opt$`focal-given-broad-probabilities-path`, row.names = 1)
focal.given.no.broad.probabilities <- read.csv(opt$`focal-given-no-broad-probabilities-path`, row.names = 1)
co.deletion.fisher <- read.csv(opt$`co-deletion-fisher-path`)
focal.and.broad <- read.csv(opt$`focal-and-broad-path`, row.names = 1)

res <- generate_gene_features(opt$`cancer-type`, opt$arm,
                         mutation.data, focal.given.broad.probabilities,
                         focal.given.no.broad.probabilities,
                         co.deletion.fisher, focal.and.broad)


# Create the output directory if it does not exist --------------------------
if (!dir.exists(opt$out)) {
  message("Creating output directory: ", opt$out)
  dir.create(opt$out, recursive = TRUE)
}

# Write the merged gene features to a CSV file ------------------------------
output_file <- here::here(opt$out, "gene_features.csv")
write.csv(res, output_file, row.names = FALSE)
message("Gene features successfully written to: ", output_file)
