#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: calculate_differential_gene_expression.R
# Description:
#   Calculates differential gene expression between aneuploid samples and normal 
#   samples for a specified cancer type using HiSeq data. The analysis is performed
#   with DESeq2. Optionally, log fold change shrinkage is applied using the apeglm 
#   method.
#
# Usage:
#   Rscript calculate_differential_gene_expression.R --cancer-type <cancer_type> 
#         --aneuploid-samples-path <path_to_csv_of_aneuploid_samples> 
#         --hiSeq-data-path <path_to_HiSeq_data_RDS> 
#         [--shrink] --out <output_folder>
#
# Output:
#   A CSV file containing the differential expression results will be written to 
#   the specified output folder.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(stringr)
  library(DESeq2)
})

#' Calculate Differential Gene Expression
#'
#' This function calculates differential gene expression between aneuploid samples 
#' and normal samples from HiSeq data for a given cancer type. Aneuploid samples 
#' are provided as a character vector. Normal samples are selected from the HiSeq 
#' data based on the sample_type "Solid Tissue Normal". The function builds a 
#' DESeq2 dataset and performs differential expression analysis. Optionally, log 
#' fold change shrinkage is applied.
#'
#' @param aneuploid.samples Character vector. Sample identifiers for aneuploid samples.
#'   These should be provided in the format of 15-character sample IDs.
#' @param shrink Logical. If TRUE, apply log fold change shrinkage using apeglm.
#' @param HiSeq.data A SummarizedExperiment object containing HiSeq expression data.
#'
#' @return A DESeq2 results object (data frame) containing differential expression 
#'   statistics.
#'
#' @examples
#' \dontrun{
#'   hiSeq <- readRDS("path/to/hiSeq_data.rds")
#'   aneuploid_ids <- c("TCGA-XX-XXXXXX", "TCGA-YY-YYYYYY")
#'   res <- calculate_differential_gene_expression("BRCA", aneuploid_ids, TRUE, hiSeq)
#' }
calculate_differential_gene_expression <- function(aneuploid.samples, shrink, HiSeq.data) {
  # Truncate sample names in HiSeq data to 15 characters for consistency
  HiSeq.data@colData$sample <- str_sub(HiSeq.data@colData$sample, 1, 15)
  
  message("Extracting normal sample IDs from HiSeq data...")
  # Get sample IDs for normal samples (solid tissue normal)
  normal.sample.ids <- subset(HiSeq.data, select = colData(HiSeq.data)$sample_type == 'Solid Tissue Normal')$sample
  
  # Process the aneuploid sample IDs: replace dots with hyphens and truncate to 15 characters
  message("Processing aneuploid sample IDs...")
  aneuploid.samples <- gsub('\\.', '-', str_sub(aneuploid.samples, 1, 15))
  
  # Subset HiSeq data to include only normal and aneuploid samples
  message("Subsetting HiSeq data to include only relevant samples...")
  HiSeq.data <- subset(HiSeq.data, select = (colData(HiSeq.data)$sample %in% c(normal.sample.ids, aneuploid.samples)))
  HiSeq.data.filtered <- subset(HiSeq.data, select = !duplicated(colData(HiSeq.data)$sample))
  
  message("Extracting raw count expression matrix...")
  # Extract raw count expression matrix from HiSeq data
  expression.matrix <- assay(HiSeq.data.filtered, "raw_count")
  
  # Remove duplicated rows based on ensembl_gene_id and set row names accordingly
  expression.matrix <- expression.matrix[!duplicated(HiSeq.data@rowRanges$ensembl_gene_id), ]
  rownames(expression.matrix) <- HiSeq.data@rowRanges$ensembl_gene_id[!duplicated(HiSeq.data@rowRanges$ensembl_gene_id)]
  
  # Create annotation vector based on whether each sample is aneuploid (group1) or not (group2)
  message("Creating sample group annotations...")
  sample.annotation.vector <- ifelse(colData(HiSeq.data.filtered)$sample %in% aneuploid.samples, 'group2', 'group1')
  group <- factor(sample.annotation.vector)
  
  message("Building DESeq2 dataset...")
  ddsMat <- DESeqDataSetFromMatrix(countData = round(expression.matrix),
                                   colData = DataFrame(group),
                                   design = ~ group)
  
  message("Performing differential expression analysis...")
  de <- DESeq(ddsMat)
  
  if (shrink) {
    message("Applying log fold change shrinkage using apeglm...")
    res <- lfcShrink(de, coef = 2, type = "apeglm")
  } else {
    res <- results(de)
  }
  
  message("Differential expression analysis complete.")
  return(res)
}

# Define command-line options -----------------------------------------------
option_list <- list(
  make_option(c("-a", "--aneuploid-samples-path"),
              type = "character",
              default = NULL,
              help = "Path to CSV file containing a list of aneuploid sample IDs (in column 'x')",
              metavar = "character"),
  make_option(c("-i", "--hiSeq-data-path"),
              type = "character",
              default = NULL,
              help = "Path to the HiSeq data file (RDS file with a SummarizedExperiment object)",
              metavar = "character"),
  make_option(c("-s", "--shrink"),
              action = "store_true",
              default = FALSE,
              help = "Apply log fold change shrinkage (default: FALSE)"),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "Output folder path for the differential expression CSV file",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that all mandatory options are provided -----------------------------
if (any(unlist(lapply(option_list, function(option) { is.null(opt[[option@dest]]) })))) {
  print_help(opt_parser)
  stop("Error: All options are mandatory", call. = FALSE)
}

message("Reading aneuploid samples from: ", opt$`aneuploid-samples-path`)
aneuploid.samples <- read.csv(opt$`aneuploid-samples-path`, row.names = 1)$x

message("Loading HiSeq data from: ", opt$`hiSeq-data-path`)
HiSeq.data <- readRDS(opt$`hiSeq-data-path`)

# Calculate differential gene expression
message("Calculating differential gene expression for cancer type: ", opt$`cancer-type`)
res <- calculate_differential_gene_expression(aneuploid.samples, opt$shrink, HiSeq.data)

# Ensure the output directory exists
if (!dir.exists(opt$out)) {
  message("Creating output directory: ", opt$out)
  dir.create(opt$out, recursive = TRUE)
}

# Write results to CSV file
output_file <- file.path(opt$out, 'aneuploid_vs_normal_deg.csv')
message("Writing differential expression results to: ", output_file)
write.csv(as.data.frame(res), output_file, row.names = TRUE)

message("Differential gene expression analysis complete.")
