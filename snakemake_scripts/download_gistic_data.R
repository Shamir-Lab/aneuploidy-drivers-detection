#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: download_gistic_data.R
# Description:
#   Downloads GISTIC results from Firebrowse for a given cancer type and extracts 
#   the content into a specified output folder. The downloaded tar.gz file contains 
#   a single folder. Its contents (files and subfolders) are copied into the output 
#   directory (without the container folder).
#
# Usage:
#   Rscript download_gistic_data.R --cancer-type BRCA --out /path/to/output
#
# Example:
#   Rscript download_gistic_data.R --cancer-type BRCA --out results/
#
# Note:
#   For example, for BRCA the expected URL is:
#   http://gdac.broadinstitute.org/runs/analyses__2016_01_28/data/BRCA-TP/20160128/
#       gdac.broadinstitute.org_BRCA-TP.CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz
#
#   In this example, the cancer type is appended with "-TP" to form the ID used in the URL.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
})

# Define command-line options
option_list <- list(
  make_option(c("-c", "--cancer-type"),
              type = "character",
              default = NULL,
              help = "Cancer type (e.g., BRCA)",
              metavar = "character"),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "Output folder path",
              metavar = "character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check mandatory arguments
if (is.null(opt$`cancer-type`) || is.null(opt$out)) {
  print_help(opt_parser)
  stop("Both --cancer-type and --out must be provided", call. = FALSE)
}

# Construct the GISTIC data URL based on the cancer type
# (Assuming the pattern appends "-TP" to the cancer type)
cancer_id <- paste0(opt$`cancer-type`, "-TP")
url <- sprintf("http://gdac.broadinstitute.org/runs/analyses__2016_01_28/data/%s/20160128/gdac.broadinstitute.org_%s.CopyNumber_Gistic2.Level_4.2016012800.0.0.tar.gz",
               cancer_id, cancer_id)

# Inform the user
message("Downloading GISTIC data from: ", url)

# Download the tar.gz file to a temporary file
temp_tar <- tempfile(fileext = ".tar.gz")
download.file(url, temp_tar, mode = "wb")

# Create a temporary directory for extraction
temp_extract_dir <- file.path(tempdir(), "gistic_extract")
dir.create(temp_extract_dir, showWarnings = FALSE)

# Extract the tar.gz file into the temporary directory
untar(temp_tar, exdir = temp_extract_dir)

# The tar is assumed to contain a single folder.
# List the directories in the extraction folder (non-recursively)
extracted_dirs <- list.dirs(temp_extract_dir, full.names = TRUE, recursive = FALSE)
if (length(extracted_dirs) < 1) {
  stop("No directory found in the extracted tar file!")
}
container_folder <- extracted_dirs[1]

# List all files and subdirectories in the container folder
files_to_copy <- list.files(container_folder, full.names = TRUE)

# Create the output folder if it doesn't exist
if (!dir.exists(opt$out)) {
  dir.create(opt$out, recursive = TRUE)
}

# Copy all contents from the container folder into the output folder
success <- file.copy(from = files_to_copy, to = opt$out, recursive = TRUE)
if (!all(success)) {
  warning("Some files may not have been copied successfully.")
} else {
  message("GISTIC data successfully extracted to: ", opt$out)
}

# Clean up: remove the temporary tar file (and optionally the extraction folder)
unlink(temp_tar)
unlink(temp_extract_dir, recursive = TRUE)
