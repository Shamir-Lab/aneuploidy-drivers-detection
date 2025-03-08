#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: generate_mutation_db.R
# Description: Downloads TCGA mutation data and saves it to a SQLite database.
# Usage: Rscript generate_mutation_db.R --out <output_file_path>
# ------------------------------------------------------------------------------

# Load required packages while suppressing startup messages -----------------
suppressPackageStartupMessages({
  library(optparse)
  library(TCGAbiolinks)
  library(dplyr)
  library(RSQLite)
})

#' Download TCGA Mutation Data and Save as SQLite Database
#'
#' This function downloads mutation data from TCGA using the TCGAbiolinks
#' package and writes the data into a SQLite database file.
#'
#' @param db_file Character. Full path to the SQLite database file.
#' @return None. The SQLite database is created at the specified location.
#' @examples
#' \dontrun{
#'   generate_mutation_db("output/TCGA_mutations.sqlite")
#' }
generate_mutation_db <- function(db_file) {
  # Download mutation data from TCGA
  message("Downloading TCGA mutation data...")
  mutation.df <- getMC3MAF()
  mutation.df <- dplyr::select(mutation.df, -Strand)
  # Connect to the SQLite database (creates file if it doesn't exist)
  message("Connecting to SQLite database...")
  db.connection <- dbConnect(RSQLite::SQLite(), db_file)
  
  # Write the mutation data into a table named "mutations"
  message("Writing mutation data to the database...")
  dbWriteTable(db.connection, "mutations", mutation.df, overwrite = TRUE)
  
  # Disconnect from the database
  dbDisconnect(db.connection)
  message("Database generation complete!")
  
  #Remove maf and maf.gz files
  message("Deleting maf and maf.gz files")
  file.remove('mc3.v0.2.8.PUBLIC.maf')
  file.remove('mc3.v0.2.8.PUBLIC.maf.gz')
}

# Define command-line options -----------------------------------------------
option_list <- list(
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "Output file path where the SQLite database file will be saved",
              metavar = "character")
)

# Parse command-line options -------------------------------------------------
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate that the mandatory 'out' parameter is provided ------------------
if (is.null(opt$out)) {
  print_help(opt_parser)
  stop("Error: The output file path (--out) is mandatory.", call. = FALSE)
}

# Ensure the output file exists (create if necessary) ---------------------
if (!dir.exists(dirname(opt$out))) {
  message("Creating output folder: ", dirname(opt$out))
  dir.create(dirname(opt$out), recursive = TRUE)
}

# Generate the mutation database -------------------------
generate_mutation_db(opt$out)