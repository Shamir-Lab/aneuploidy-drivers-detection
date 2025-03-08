#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: generate_aneuploidy_samples.R
# Description:
#   Loads arm deletion data for TCGA samples and identifies samples with aneuploidy 
#   (arm deletions) and those without for a specified cancer type and chromosome arm.
#
# Usage:
#   Rscript generate_aneuploidy_samples.R --cancer-type <cancer_type> 
#        --arm <arm> --arm-deletion-data-path <path_to_arm_deletion_csv> 
#        --out <output_folder>
#
# Output:
#   Two CSV files will be written to the output folder:
#     - aneuploid_samples.csv: Samples with arm deletion (aneuploid)
#     - non_aneuploid_samples.csv: Samples without arm deletion
# ------------------------------------------------------------------------------

# Load required libraries and set working directory -------------------------
suppressPackageStartupMessages({
  library(here)
  set_here(path=normalizePath('..'))
  source(here::here('utils.R'))
  
  # Load required packages
  library(optparse)
  library(dplyr)
  library(purrr)
  library(foreach)
})

#' Generate Aneuploidy Sample Lists
#'
#' Given a cancer type, chromosome arm, and arm deletion data, this function
#' identifies samples with a deletion in the specified arm ("aneuploid") and
#' those without such a deletion.
#'
#' @param cancer.type_ Character. The cancer type (e.g., "BRCA").
#' @param arm_ Character. The chromosome arm name (e.g., "X3p").
#' @param arm_deletion_data A data.frame with arm deletion data per sample.
#'
#' @return A list with two elements:
#'         - del.samples: vector of sample identifiers with arm deletion.
#'         - non.del.samples: vector of sample identifiers without arm deletion.
#'
#' @examples
#' \dontrun{
#'   data <- read.csv("path/to/arm_deletion_data.csv")
#'   res <- generate_mut_data("BRCA", "X3p", data)
#' }
generate_mut_data <- function(cancer.type_, arm_, arm_deletion_data) {
  # Use the provided arm deletion data instead of loading from a constant CSV
  anp.data <- arm_deletion_data
  
  # Identify samples with deletion in the specified arm (value == -1)
  del.samples <- dplyr::filter(anp.data, 
                               Type == cancer.type_,
                               !is.na(get({{ arm_ }})),
                               get({{ arm_ }}) == -1) %>% 
    .$sample
  
  # Identify samples without deletion in the specified arm (value != -1)
  non.del.samples <- dplyr::filter(anp.data, 
                                   Type == cancer.type_,
                                   !is.na(get({{ arm_ }})),
                                   get({{ arm_ }}) != -1) %>% 
    .$sample
  
  return(list(del.samples = del.samples, non.del.samples = non.del.samples))
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
  make_option(c("-d", "--arm-deletion-data-path"),
              type = "character",
              default = NULL,
              help = "Path to CSV with arm deletion data per sample",
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

# Validate that all mandatory options are provided --------------------------
if (any(unlist(lapply(option_list, function(option) {
  is.null(opt[[option@dest]])
})))) {
  print_help(opt_parser)
  stop("Error: All options are mandatory", call. = FALSE)
}

# Load arm deletion data from the specified CSV file --------------------------
arm_deletion_data <- load_ANP_data(opt$`arm-deletion-data-path`)$tcga.anp.data

# Generate sample lists for the specified cancer type and chromosome arm ------
res <- generate_mut_data(opt$`cancer-type`, opt$arm, arm_deletion_data)

# Ensure the output directory exists -----------------------------------------
if (!dir.exists(opt$out)) {
  message("Creating output directory: ", opt$out)
  dir.create(opt$out, recursive = TRUE)
}

# Write results to CSV files -----------------------------------------------
output_file_del <- here::here(opt$out, 'aneuploid_samples.csv')
output_file_non_del <- here::here(opt$out, 'non_aneuploid_samples.csv')

write.csv(res$del.samples, output_file_del, row.names = TRUE)
write.csv(res$non.del.samples, output_file_non_del, row.names = TRUE)

message("Aneuploid samples written to: ", output_file_del)
message("Non-aneuploid samples written to: ", output_file_non_del)
