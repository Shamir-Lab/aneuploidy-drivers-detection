#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: calculate_anp_effect.R
# Description:
#   Calculates the effect of arm-level deletions on gene expression using raw 
#   count data. The analysis uses inverse probability of treatment weighting (IPTW)
#   to adjust for differences in aneuploidy scores. The input raw counts are assumed 
#   to be stored in an RDS file containing a SummarizedExperiment object.
#
#   The script loads additional ANP data, filters for protein-coding 
#   genes in primary tumors, and then computes a weighted t-test for each gene to 
#   estimate effect sizes. The results (including p-values, effect sizes, and adjusted 
#   q-values) are written to a CSV file.
#
# Usage:
#   Rscript calculate_anp_effect.R --cancer-type <cancer_type> --arm <arm> 
#         --raw-counts-path <path_to_raw_counts_RDS>
#         --out <output_folder>
#
# Example:
#   Rscript calculate_anp_effect.R --cancer-type BRCA --arm X3p \
#         --raw-counts-path data/BRCA_raw_counts.rds \
#         --out results/
#
# Output:
#   A CSV file named 'IPTW_anp_effect.csv' is written to the specified output folder.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  source(here::here('code/utils.R'))
  source(here::here('code/genome_annotation_utils.R'))
  
  library(optparse)
  library(tidyr)
  library(dplyr)
  library(foreach)
  library(effsize)
  library(SummarizedExperiment)
  library(weights)
  library(sjstats)
})

#' Calculate ANP Effect on Gene Expression
#'
#' This function calculates the effect of arm-level deletion on gene expression 
#' in a specified cancer type using raw count data from a SummarizedExperiment object.
#' It performs the following steps:
#'   1. Filters the raw counts to include only protein-coding genes in primary tumors.
#'   2. Computes FPKM values and applies a log10 transformation.
#'   3. Loads additional ANP data and filters it for the specified 
#'      cancer type and arm (excluding samples where the arm value is 1).
#'   4. Fits a logistic regression model (propensity model) and computes stabilized weights.
#'   5. Merges the expression data with ANP data and calculates effect sizes (Cohen's d)
#'      and p-values using a weighted t-test for each gene.
#'
#' @param cancer.type Character. The cancer type (e.g., "BRCA").
#' @param arm Character. The chromosome arm name (e.g., "X3p").
#' @param raw.counts A SummarizedExperiment object containing raw counts and FPKM values.
#'
#' @return A data frame with columns:
#'         - name: Gene identifier (external gene name).
#'         - pval: p-value from the weighted t-test.
#'         - effect.size: Cohen's d effect size.
#'         - qval: Adjusted p-value (FDR).
#'
#' @examples
#' \dontrun{
#'   se <- readRDS("data/BRCA_raw_counts.rds")
#'   effect_df <- get_anp_effect("BRCA", "X3p", se)
#' }
get_anp_effect <- function(cancer.type, arm, raw.counts) {
  
  message("Filtering raw counts for protein-coding genes in primary tumors...")
  raw.counts.tumors <- raw.counts[
    raw.counts@rowRanges$gene_type == 'protein_coding',
    raw.counts@colData$sample_type == 'Primary Tumor'
  ]
  
  message("Extracting FPKM values and computing log10-transformed expression...")
  fpkm <- assay(raw.counts.tumors, "fpkm_unstrand")
  fpkm <- fpkm[! (rownames(fpkm) %>% str_sub(1, 15) %>% duplicated()), ]
  rownames(fpkm) <- rownames(fpkm) %>% str_sub(1, 15)
  log10_fpkm <- log10(fpkm + 1)
  
  message("Loading ANP data")
  anp.data <- load_ANP_data()$tcga.anp.data
  anp.data <- anp.data %>% 
    dplyr::filter(Type == cancer.type, 1 != !!as.name(arm)) %>% 
    dplyr::mutate(arm.del = (-1 == !!as.name(arm))) %>%
    dplyr::select(sample, AneuploidyScore, arm.del)
  
  message("Fitting propensity model to compute stabilized weights...")
  propensity_model <- glm(arm.del ~ AneuploidyScore, data = anp.data, family = "binomial")
  propensity_scores <- predict(propensity_model, type = "response")
  iptw <- ifelse(anp.data$arm.del == 1, 1 / propensity_scores, 1 / (1 - propensity_scores))
  stabilized_weights <- ifelse(anp.data$arm.del, iptw * mean(anp.data$arm.del), iptw * mean(anp.data$arm.del == F))
  anp.data$w <- stabilized_weights
  
  message("Mapping gene identifiers to external gene names...")
  gns <- map_gene_identifiers(rownames(log10_fpkm), 'ensembl_gene_id', 'external_gene_name')
  log10_fpkm <- log10_fpkm[!duplicated(gns), ]
  rownames(log10_fpkm) <- gns[!duplicated(gns)]
  
  message("Merging expression data with ANP data...")
  df <- data.frame(t(log10_fpkm))
  df$sample <- gsub('\\.', '-', rownames(df)) %>% str_sub(1, 15)
  df <- merge(anp.data, df, by = 'sample', all = FALSE)
  
  message("Calculating effect sizes for each gene using weighted t-tests...")
  arm.del.effect <- foreach(gene = colnames(df)[-c(1, 2, 3, 4)], .init = data.frame(), .combine = rbind) %do% {
    tryCatch({
      x <- df[df$arm.del == TRUE, gene]
      y <- df[df$arm.del == FALSE, gene]
      w_x <- df$w[df$arm.del == TRUE]
      w_y <- df$w[df$arm.del == FALSE]
      p.val <- unname(wtd.t.test(x, y, w_x, w_y)$coefficients['p.value'])
      w_sd_x <- weighted_sd(x, w_x)
      w_sd_y <- weighted_sd(y, w_y)
      n_x <- length(x)
      n_y <- length(y)
      c.d <- (wtd.mean(x, w_x) - wtd.mean(y, w_y)) / ((((n_x - 1) * w_sd_x^2 + (n_y - 1) * w_sd_y^2) / (n_x + n_y - 2))^(1/2))
      data.frame(gene = gene, p.val = p.val, c.d = c.d)
    }, error = function(cond) {
      message(paste0('Error working on ', gene, ': ', cond))
      return(c())
    })
  }
  
  colnames(arm.del.effect) <- c('name', 'pval', 'effect.size')
  arm.del.effect$pval <- as.numeric(arm.del.effect$pval)
  arm.del.effect$effect.size <- as.numeric(arm.del.effect$effect.size)
  arm.del.effect$qval <- p.adjust(arm.del.effect$pval, method = 'BH')
  
  message("ANP effect calculation complete. Returning results.")
  return(arm.del.effect)
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
  make_option(c("-r", "--raw-counts-path"),
              type = "character",
              default = NULL,
              help = "Path to the raw counts file (RDS file with a SummarizedExperiment object)",
              metavar = "character"),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "Output folder path for the results CSV file",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check that all mandatory options are provided
if (any(unlist(lapply(option_list, function(option) { is.null(opt[[option@dest]]) })))) {
  print_help(opt_parser)
  stop("Error: All options are mandatory", call. = FALSE)
}

message("Loading raw counts data from: ", opt$`raw-counts-path`)
raw_counts <- readRDS(opt$`raw-counts-path`)

message("Calculating ANP effect for cancer type: ", opt$`cancer-type`, " and arm: ", opt$`arm`)
res <- get_anp_effect(opt$`cancer-type`, opt$`arm`, raw_counts)

# Ensure the output directory exists
if (!dir.exists(opt$out)) {
  message("Creating output directory: ", opt$out)
  dir.create(opt$out, recursive = TRUE)
}

output_file <- here::here(opt$out, 'IPTW_anp_effect.csv')
message("Writing results to: ", output_file)
write.csv(res, output_file, row.names = FALSE)

message("Script completed successfully.")
