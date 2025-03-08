#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: generate_mutation_data.R
# Description:
#   This script reads a mutation matrix and sample lists (for broad deletion and
#   non-broad deletion samples) from CSV files, computes Fisher's exact test
#   statistics to assess differences in mutation patterns between these groups,
#   and outputs the results as a CSV file.
#
# Usage:
#   Rscript generate_mutation_data.R --mutation-matrix-path <path_to_mutation_matrix_csv> 
#        --broad-del-samples-path <path_to_broad_del_samples_csv> 
#        --no-broad-del-samples-path <path_to_no_broad_del_samples_csv> 
#        --out <output_folder>
# ------------------------------------------------------------------------------

# Load required libraries and set working directory -------------------------
suppressPackageStartupMessages({
  library(here)
  set_here(path=normalizePath('..'))
  source(here::here('genome_annotation_utils.R'))
  
  # Load necessary libraries
  library(optparse)
  library(dplyr)
  library(purrr)
  library(foreach)
})

#' Generate Mutation Data Using Fisher's Exact Test
#'
#' This function processes a mutation matrix and computes Fisher's exact test
#' for each gene to evaluate whether differences in mutation frequencies between
#' broad deletion samples and non-broad deletion samples are statistically significant.
#'
#' @param mut.events A data.frame of mutation events (genes as rows, samples as columns).
#' @param del.samples A character vector of sample identifiers with broad deletions.
#' @param no.del.samples A character vector of sample identifiers without broad deletions.
#'
#' @return A data.frame containing, for each gene, Fisher test statistics, counts, 
#'         mutation probabilities, and the external gene name.
#'
#' @examples
#' \dontrun{
#'   res <- generate_mut_data(mutation_matrix, broad_del_samples, no_broad_del_samples)
#' }
generate_mut_data <- function(mut.events, del.samples, no.del.samples) {
  # Clean column names: replace dots with hyphens and truncate to 15 characters
  colnames(mut.events) <- gsub('\\.', '-', str_sub(colnames(mut.events), 1, 15))
  
  # Select columns corresponding to deletion and non-deletion samples
  mut.events.del <- dplyr::select(mut.events, any_of(del.samples))
  mut.events.no.del <- dplyr::select(mut.events, any_of(no.del.samples))
  
  # Compute Fisher's exact test for each gene to assess mutation patterns
  mut.fisher.df <- data.frame(
    foreach(i = 1:nrow(mut.events), .init = c(), .combine = rbind) %do% {
      # Construct contingency table for the current gene
      tbl <- data.frame(
        no.broad = c(mut = length(which(mut.events.no.del[i, ] > 0)),
                     no.mut = length(which(mut.events.no.del[i, ] == 0))),
        broad = c(mut = length(which(mut.events.del[i, ] > 0)),
                  no.mut = length(which(mut.events.del[i, ] == 0)))
      )
      
      # Perform Fisher's exact test with a pseudocount of 1
      res <- fisher.test(tbl + 1)
      mut.depletion <- fisher.test(tbl + 1, alternative = 'greater')
      
      # Return computed statistics and counts
      c(
        pval = res$p.value,
        odds.ratio = unname(res$estimate),
        mutation.and.no.broad = tbl['mut', 'no.broad'],
        no.mutation.and.no.broad = tbl['no.mut', 'no.broad'],
        mutation.and.broad = tbl['mut', 'broad'],
        no.mutation.and.broad = tbl['no.mut', 'broad'],
        mutation.depletion.pvalue = mut.depletion$p.value
      )
    }
  )
  
  # Add gene identifiers to the results
  mut.fisher.df$gene <- rownames(mut.events)
  
  # Calculate mutation probabilities for both sample groups
  mut.fisher.df <- mut.fisher.df %>% 
    dplyr::mutate(
      mut.given.no.brod.proba = mutation.and.no.broad / (mutation.and.no.broad + no.mutation.and.no.broad),
      mut.given.brod.proba = mutation.and.broad / (mutation.and.broad + no.mutation.and.broad)
    )
  
  # Map internal gene identifiers to external gene names
  mut.fisher.df$name <- map_gene_identifiers(mut.fisher.df$gene, 'ensembl_gene_id', 'external_gene_name')
  
  return(mut.fisher.df)
}

# Define command-line options -----------------------------------------------
option_list <- list(
  make_option(c("-m", "--mutation-matrix-path"),
              type = "character",
              default = NULL,
              help = "Path to CSV file containing the mutation matrix",
              metavar = "character"),
  make_option(c("-b", "--broad-del-samples-path"),
              type = "character",
              default = NULL,
              help = "Path to CSV file with broad deletion sample names",
              metavar = "character"),
  make_option(c("-n", "--no-broad-del-samples-path"),
              type = "character",
              default = NULL,
              help = "Path to CSV file with non-broad deletion sample names",
              metavar = "character"),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "Output folder path for saving the processed mutation data CSV",
              metavar = "character")
)

# Parse command-line options -----------------------------------------------
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that all mandatory options have been provided
if (any(unlist(lapply(option_list, function(option) { is.null(opt[[option@dest]]) })))) {
  print_help(opt_parser)
  stop("Error: All options are mandatory", call. = FALSE)
}

# Read input CSV files --------------------------------------------------------
del.samples <- read.csv(opt$`broad-del-samples-path`, row.names = 1)$x
no.del.samples <- read.csv(opt$`no-broad-del-samples-path`, row.names = 1)$x
mutation.matrix <- read.csv(opt$`mutation-matrix-path`, row.names = 1)

# Generate mutation data using the defined function -------------------------
res <- generate_mut_data(mutation.matrix, del.samples, no.del.samples)

# Ensure the output directory exists ----------------------------------------
if (!dir.exists(opt$out)) {
  message("Creating output directory: ", opt$out)
  dir.create(opt$out, recursive = TRUE)
}

# Write the results to a CSV file in the output directory --------------------
output_file <- here::here(opt$out, 'mutation_data.csv')
write.csv(res, output_file, row.names = FALSE)
message("Mutation data successfully written to: ", output_file)
