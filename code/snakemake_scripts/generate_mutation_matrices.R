#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: generate_mutation_matrix.R
# Description: Generates a mutation matrix for a given cancer type and chromosome
#              arm using mutation data stored in a SQLite database.
# Usage: Rscript generate_mutation_matrix.R --cancer-type <cancer_type> --arm <arm>
#              --mutation-db-path <path_to_sqlite_db> --out <output_folder>
# ------------------------------------------------------------------------------

# Load required libraries ----------------------------------------------------
suppressPackageStartupMessages({
  library(here)
  source(here::here('code/genome_annotation_utils.R'))
  library(optparse)
  library(dplyr)
  library(RSQLite)
  library(foreach)
})

#' Generate mutation matrix for a given cancer type and arm
#'
#' This function retrieves mutation data for a specified cancer type and
#' chromosome arm from a SQLite database and generates a mutation matrix.
#' Rows represent genes and columns represent tumor sample barcodes.
#'
#' @param cancer.type Character. The cancer type (e.g., "BRCA").
#' @param arm Character. The chromosome arm name (e.g., "X3p").
#' @param db_file Character. Full path to the SQLite database file containing mutation data.
#'
#' @return A data.frame representing the mutation matrix.
#' @examples
#' \dontrun{
#'   generate_mutation_matrix("BRCA", "X3p", "path/to/MC3MAF.sqlite")
#' }
generate_mutation_matrix <- function(cancer.type, arm, db_file) {
  # Create project identifier based on cancer type
  project_id_ <- paste0("TCGA-", cancer.type)
  
  # Retrieve protein-coding genes on the specified chromosome arm
  arm.genes <- filter(get_arm_genes(arm, GRCh.version = "37"),
                      transcript_biotype == "protein_coding")
  genes.df <- data.frame(Gene = unique(arm.genes$ensembl_gene_id))
  
  # Connect to the SQLite database using the provided db_file path
  message("Connecting to the SQLite database...")
  db.connection <- dbConnect(RSQLite::SQLite(), db_file)
  
  # Reference the mutations table
  mutations.table <- tbl(db.connection, "mutations")
  
  # Filter mutations for the specified genes and cancer type
  arm.mutations <- dplyr::filter(mutations.table,
                                 Gene %in% !!genes.df$Gene,
                                 project_id == project_id_) %>% collect()
  
  # Retrieve distinct tumor sample barcodes for the project
  project.samples <- dplyr::filter(mutations.table, project_id == project_id_) %>%
    dplyr::select("Tumor_Sample_Barcode") %>%
    distinct() %>%
    collect() %>%
    .$Tumor_Sample_Barcode
  
  # Disconnect from the database
  dbDisconnect(db.connection)
  
  # Filter out mutations with low impact
  arm.mutations.non.low <- dplyr::filter(arm.mutations, IMPACT %in% c("HIGH", "MODERATE"))
  
  # Build the mutation matrix: rows are genes, columns are tumor samples
  mut.events <- foreach(sample = project.samples,
                        .init = data.frame(matrix(nrow = nrow(genes.df), ncol = 0)),
                        .combine = cbind) %do% {
                          mut.count <- arm.mutations.non.low %>%
                            dplyr::filter(Tumor_Sample_Barcode == sample) %>%
                            dplyr::count(Gene) %>%
                            data.frame()
                          
                          # Merge to ensure all genes are represented; fill missing with 0
                          x <- merge(genes.df, mut.count, by = "Gene", all.x = TRUE)
                          x[is.na(x)] <- 0
                          x[order(x$Gene),]$n
                        }
  
  colnames(mut.events) <- project.samples
  rownames(mut.events) <- sort(genes.df$Gene)
  
  return(mut.events)
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
  make_option(c("-d", "--mutation-db-path"),
              type = "character",
              default = NULL,
              help = "Path to the SQLite database file containing mutation data",
              metavar = "character"),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "Output folder path to save the mutation matrix CSV",
              metavar = "character")
)

# Parse command-line options -------------------------------------------------
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Validate that all mandatory parameters are provided -----------------------
if (any(sapply(option_list, function(option) is.null(opt[[option@dest]])))) {
  print_help(opt_parser)
  stop("Error: All options (--cancer-type, --arm, --mutation-db-path, --out) are mandatory.", 
       call. = FALSE)
}

# Generate the mutation matrix with error handling --------------------------
res <- generate_mutation_matrix(opt$`cancer-type`, opt$arm, opt$`mutation-db-path`)

# Ensure the output folder exists (create if necessary) ---------------------
if (!dir.exists(opt$out)) {
  message("Creating output folder: ", opt$out)
  dir.create(opt$out, recursive = TRUE)
}

# Write the mutation matrix to a CSV file -----------------------------------
output_file <- here::here(opt$out, 'mutation_matrix.csv')
write.csv(res, output_file, row.names = TRUE)
message("Mutation matrix successfully written to: ", output_file)
