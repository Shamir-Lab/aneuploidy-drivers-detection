#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: filter_mutations_data.R
# Description:
#   Connects to an SQLite database containing mutation data, retrieves the
#   "mutations" table, filters the data based on a list of sample IDs provided
#   in a CSV file, and writes the filtered data to a TSV file.
#
# Usage:
#   Rscript filter_mutations_data.R --db-file <path_to_sqlite_db> 
#         --aneuploid-samples-path <path_to_aneuploid_csv> 
#         --non-aneuploid-samples-path <path_to_non_aneuploid_csv> 
#         --out <output_folder>
#
#   The samples CSV should contain a column named "x" with sample IDs.
#
# Output:
#   A TSV file named "filtered_mutations.tsv" will be written to the specified 
#   output folder.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(optparse)
  library(dplyr)
  library(DBI)
  library(RSQLite)
  library(stringr)
})

# Define command-line options
option_list <- list(
  make_option(c("-d", "--db-file"),
              type = "character",
              default = NULL,
              help = "Path to the SQLite database file containing mutations data.",
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
              help = "Output folder path for the filtered TSV file.",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that all mandatory options are provided
if (any(unlist(lapply(option_list, function(option) { is.null(opt[[option@dest]]) })))) {
  print_help(opt_parser)
  stop("Error: All options are mandatory", call. = FALSE)
}

# Connect to the SQLite database
message("Connecting to SQLite database: ", opt$`db-file`)
db.connection <- dbConnect(RSQLite::SQLite(), opt$`db-file`)

# Reference the mutations table
message("Reading mutations table...")
mutations.table <- tbl(db.connection, "mutations")

# Read the sample lists from the provided CSV files.
# Assumes the CSV files have a single column with header 'x'
message("Reading aneuploid sample list from: ", opt$`aneuploid-samples-path`)
aneuploid_samples <- read.csv(opt$`aneuploid-samples-path`, row.names = 1)$x

message("Reading non-aneuploid sample list from: ", opt$`non-aneuploid-samples-path`)
non_aneuploid_samples <- read.csv(opt$`non-aneuploid-samples-path`, row.names = 1)$x

# Filter the mutations table by sample IDs (assumed column name: Tumor_Sample_Barcode)
message("Filtering mutations for specified sample IDs...")
filtered_mutations_aneuploid <- mutations.table %>%
  filter(str_sub(Tumor_Sample_Barcode,1,15) %in% aneuploid_samples) %>%
  collect()

filtered_mutations_non_aneuploid <- mutations.table %>%
  filter(str_sub(Tumor_Sample_Barcode,1,15) %in% non_aneuploid_samples) %>%
  collect()

# Disconnect from the database
dbDisconnect(db.connection)
message("Disconnected from database.")

# Ensure the output directory exists
if (!dir.exists(opt$out)) {
  message("Creating output directory: ", opt$out)
  dir.create(opt$out, recursive = TRUE)
}

# Define output file paths
aneuploid_out <- here::here(opt$out, "aneuploid_mutations.tsv")
non_aneuploid_out <- here::here(opt$out, "non_aneuploid_mutations.tsv")

# Write the filtered mutations files
message("Writing aneuploid muations file to: ", aneuploid_out)
write.table(filtered_mutations_aneuploid, file = aneuploid_out, sep = "\t", row.names = FALSE, quote = FALSE)

message("Writing non-aneuploid segmentation file to: ", non_aneuploid_out)
write.table(filtered_mutations_non_aneuploid, file = non_aneuploid_out, sep = "\t", row.names = FALSE, quote = FALSE)

message("Filtering complete. Files written to: ", opt$out)
