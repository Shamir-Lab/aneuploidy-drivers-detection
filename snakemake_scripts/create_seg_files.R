#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: filter_seg_files.R
# Description:
#   This script reads an input segmentation (.seg) file and two CSV files that 
#   contain lists of sample IDs for aneuploid and non-aneuploid samples. It then 
#   creates two new segmentation files: one containing only the rows for samples 
#   in the aneuploid list, and another containing only rows for samples in the 
#   non-aneuploid list.
#
# Usage:
#   Rscript filter_seg_files.R --seg-path <path_to_seg_file> 
#         --aneuploid-samples-path <path_to_aneuploid_csv> 
#         --non-aneuploid-samples-path <path_to_non_aneuploid_csv> 
#         --out <output_folder>
#
# Example:
#   Rscript filter_seg_files.R --seg-path data/TCGA-LAML.seg \
#         --aneuploid-samples-path data/aneuploid.csv \
#         --non-aneuploid-samples-path data/non_aneuploid.csv \
#         --out results/
#
# Output:
#   Two files are created in the output folder:
#     - aneuploid.seg: Contains only rows for samples in the aneuploid sample list.
#     - non_aneuploid.seg: Contains only rows for samples in the non-aneuploid sample list.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  set_here(path = normalizePath('..'))
  library(optparse)
  library(dplyr)
  library(stringr)
})

# Define command-line options
option_list <- list(
  make_option(c("--seg-path"),
              type = "character",
              default = NULL,
              help = "Path to the input segmentation (.seg) file.",
              metavar = "character"),
  make_option(c("--aneuploid-samples-path"),
              type = "character",
              default = NULL,
              help = "Path to CSV file containing the list of aneuploid sample IDs.",
              metavar = "character"),
  make_option(c("--non-aneuploid-samples-path"),
              type = "character",
              default = NULL,
              help = "Path to CSV file containing the list of non-aneuploid sample IDs.",
              metavar = "character"),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "Output folder path for the filtered segmentation files.",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that all mandatory options are provided
if (any(unlist(lapply(option_list, function(option) { is.null(opt[[option@dest]]) })))) {
  print_help(opt_parser)
  stop("Error: All options are mandatory", call. = FALSE)
}

# Read the segmentation file
message("Reading segmentation file from: ", opt$`seg-path`)
seg_data <- read.delim(opt$`seg-path`, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Read the sample lists from the provided CSV files.
# Assumes the CSV files have a single column with header 'x'
message("Reading aneuploid sample list from: ", opt$`aneuploid-samples-path`)
aneuploid_samples <- read.csv(opt$`aneuploid-samples-path`, row.names = 1)$x

message("Reading non-aneuploid sample list from: ", opt$`non-aneuploid-samples-path`)
non_aneuploid_samples <- read.csv(opt$`non-aneuploid-samples-path`, row.names = 1)$x

# Filter the segmentation data based on the sample lists
message("Filtering segmentation data for aneuploid samples...")
seg_aneuploid <- seg_data %>% filter(str_sub(Sample,1,15) %in% aneuploid_samples)

message("Filtering segmentation data for non-aneuploid samples...")
seg_non_aneuploid <- seg_data %>% filter(str_sub(Sample,1,15) %in% non_aneuploid_samples)

# Create the output directory if it doesn't exist
if (!dir.exists(opt$out)) {
  message("Creating output directory: ", opt$out)
  dir.create(opt$out, recursive = TRUE)
}

# Define output file paths
aneuploid_out <- here::here(opt$out, "aneuploid.seg")
non_aneuploid_out <- here::here(opt$out, "non_aneuploid.seg")

# Write the filtered segmentation files
message("Writing aneuploid segmentation file to: ", aneuploid_out)
write.table(seg_aneuploid, file = aneuploid_out, sep = "\t", row.names = FALSE, quote = FALSE)

message("Writing non-aneuploid segmentation file to: ", non_aneuploid_out)
write.table(seg_non_aneuploid, file = non_aneuploid_out, sep = "\t", row.names = FALSE, quote = FALSE)

message("Filtering complete. Files written to: ", opt$out)
