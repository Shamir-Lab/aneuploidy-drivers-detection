#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: generate_cnv_matrices.R
# Description:
#   Generates copy-number variation (CNV) matrices for a specified cancer type 
#   and chromosome arm using segmentation and sample cutoff data provided by the user.
#
# Usage:
#   Rscript generate_cnv_matrices.R --cancer-type <cancer_type> --arm <arm> 
#         --seg-file-path <path_to_segmentation_file> --sample-cutoffs-path <path_to_sample_cutoffs_file> --out <output_folder>
#
# Output:
#   Several CSV files with CNV matrices and derived probabilities will be saved
#   in the specified output folder.
# ------------------------------------------------------------------------------

# Load required libraries and set working directory -------------------------
suppressPackageStartupMessages({
  library(here)
  set_here(path=normalizePath('..'))
  source(here::here('genome_annotation_utils.R'))
  
  # Load required packages
  library(optparse)
  library(dplyr)
  library(GenomicRanges)
  library(plyranges)
  library(foreach)
})

#' Generate CNV Matrices for a Given Cancer Type and Chromosome Arm
#'
#' This function reads segmentation data and sample cutoff information provided 
#' by the user to generate CNV matrices. It computes matrices for both focal and 
#' broad CNVs, as well as associated loss-of-heterozygosity (LOH) events, and 
#' performs a Fisher's exact test to assess enrichment of focal events given broad events.
#'
#' @param cancer.type Character. Cancer type identifier (e.g., "BRCA").
#' @param arm Character. Chromosome arm name (e.g., "X3p").
#' @param seg_file_path Character. Path to the segmentation file.
#' @param sample_cutoffs_path Character. Path to the sample cutoffs file.
#' @param focal.threshold.prec Numeric. Proportion of arm length to define focal events.
#'        Default is 0.7.
#'
#' @return A list containing:
#'         - short.mat: focal CNV matrix.
#'         - long.mat: broad CNV matrix.
#'         - short.loh.mat: focal LOH CNV matrix.
#'         - long.loh.mat: broad LOH CNV matrix.
#'         - short.given.long.cnv.probabilities: probability matrix for focal CNVs
#'           given broad CNVs.
#'         - short.given.no.long.cnv.probabilities: probability matrix for focal CNVs
#'           given absence of broad CNVs.
#'         - del.fisher.df: Data frame with Fisher test statistics for each gene.
#'         - short.and.long: Data frame with counts of genes having both focal and broad events.
#'
#' @examples
#' \dontrun{
#'   res <- generate_cnv_matrices("BRCA", "X3p", "/path/to/seg_file.txt", "/path/to/sample_cutoffs.txt")
#' }
generate_cnv_matrices <- function(cancer.type, arm, seg_file_path, sample_cutoffs_path, focal.threshold.prec = 0.7) {
  
  message("Starting CNV matrix generation for cancer type: ", cancer.type, " and arm: ", arm)
  
  # Get protein-coding genes on the specified arm
  message("Retrieving protein-coding genes for arm: ", arm)
  arm.genes <- dplyr::filter(get_arm_genes(arm, GRCh.version = "37"),
                             transcript_biotype == "protein_coding")
  
  # Read segmentation data from the provided segmentation file
  message("Reading segmentation file from: ", seg_file_path)
  seg.data <- read.table(seg_file_path, header = TRUE, sep = "\t")
  
  # Read sample cutoffs from the provided sample cutoff file
  message("Reading sample cutoff file from: ", sample_cutoffs_path)
  sample.cutoffs <- read.csv(sample_cutoffs_path, header = TRUE, sep = "\t", skip = 1)
  
  # Truncate sample names to 15 characters
  message("Truncating sample names in segmentation data.")
  seg.data$Sample <- substr(seg.data$Sample, 1, 15)
  g.ranges <- GRanges(seg.data)
  
  # Determine chromosome and arm boundaries based on arm name
  message("Determining chromosome and arm boundaries for arm: ", arm)
  chromosome <- arm_to_chromosome(arm)
  chromosome.info <- get_chromosome_info_hg19() %>% dplyr::filter(chrom == chromosome)
  if (substr(arm, nchar(arm), nchar(arm)) == "p") {
    arm.start <- 1
    arm.end <- chromosome.info$centromerStart
  } else {
    arm.start <- chromosome.info$centromerEnd
    arm.end <- chromosome.info$length
  }
  
  # Create genomic ranges for genes on the specified arm
  message("Creating genomic ranges for genes on the specified arm.")
  gene.ranges <- GRanges(arm.genes %>% dplyr::rename(Start = start_position, End = end_position))
  
  # Process each patient sample
  message("Processing patient samples...")
  dels <- foreach(patient = unique(substr(unique(g.ranges$Sample), 1, 12)), .init = c(), .combine = rbind) %do% {
    tryCatch({
      sample <- paste0(patient, "-01")
      message("Processing sample: ", sample)
      dell.cutoff <- -0.1
      
      # Retrieve cutoff for current sample; if not found, use a default value
      loh.cutoff <- dplyr::filter(sample.cutoffs, substr(Sample, 1, 15) == sample)$Low
      if (!sample %in% substr(sample.cutoffs$Sample, 1, 15) || loh.cutoff == -0.1) {
        message("No valid cutoff found for sample: ", sample, "; setting LOH cutoff to -10")
        loh.cutoff <- -10
      }
      
      res <- c()
      # Calculate LOH regions
      loh <- GenomicRanges::reduce(dplyr::filter(g.ranges[seqnames(g.ranges) == chromosome],
                                                 Sample == sample,
                                                 Segment_Mean < loh.cutoff), min.gapwidth = 400000)
      # Calculate deletion regions
      del <- GenomicRanges::reduce(dplyr::filter(g.ranges[seqnames(g.ranges) == chromosome],
                                                 Sample == sample,
                                                 Segment_Mean < dell.cutoff), min.gapwidth = 400000)
      
      # Limit regions to the boundaries of the arm
      del <- dplyr::filter(del, (start >= arm.start & start <= arm.end) | (end >= arm.start & end <= arm.end)) %>%
        dplyr::mutate(start = pmax(start, arm.start), end = pmin(end, arm.end))
      loh <- dplyr::filter(loh, start >= arm.start, start <= arm.end) %>%
        dplyr::mutate(start = pmax(start, arm.start), end = pmin(end, arm.end))
      
      # Join overlapping regions with gene ranges and annotate type
      if (length(del) > 0) {
        res <- rbind(res,
                     plyranges::join_overlap_inner(del, gene.ranges) %>%
                       data.frame() %>%
                       dplyr::select(ensembl_gene_id, start, end, width) %>%
                       unique() %>%
                       dplyr::mutate(sample = sample, type = "del"))
      }
      if (length(loh) > 0) {
        res <- rbind(res,
                     plyranges::join_overlap_inner(loh, gene.ranges) %>%
                       data.frame() %>%
                       dplyr::select(ensembl_gene_id, start, end, width) %>%
                       unique() %>%
                       dplyr::mutate(sample = sample, type = "loh"))
      }
      res
    }, error = function(cond) {
      message("Error processing sample: ", sample)
      print(cond)
      c()
    })
  }
  
  # Define the set of all samples based on the segmentation data
  all.samples <- paste0(unique(substr(unique(g.ranges$Sample), 1, 12)), "-01")
  message("Total number of samples processed: ", length(all.samples))
  
  # Define focal threshold as a proportion of the arm length
  focal.threshold <- focal.threshold.prec * (arm.end - arm.start)
  message("Focal threshold set to: ", focal.threshold)
  
  # Initialize matrices for focal and broad CNVs and LOH events
  unique.genes <- unique(dels$ensembl_gene_id)
  short.mat <- long.mat <- short.loh.mat <- long.loh.mat <-
    matrix(0, nrow = length(unique.genes), ncol = length(all.samples))
  rownames(short.mat) <- rownames(long.mat) <- rownames(short.loh.mat) <- rownames(long.loh.mat) <- unique.genes
  colnames(short.mat) <- colnames(long.mat) <- colnames(short.loh.mat) <- colnames(long.loh.mat) <- all.samples
  
  message("Filling CNV matrices for each sample...")
  # Fill matrices based on event lengths compared to focal threshold
  for (smp in all.samples) {
    total.del.len <- sum(dplyr::select(dplyr::filter(dels, type == "del", sample == smp), start, end, width) %>%
                           unique() %>% .$width)
    total.loh.len <- sum(dplyr::select(dplyr::filter(dels, type == "loh", sample == smp), start, end, width) %>%
                           unique() %>% .$width)
    
    if (total.del.len < focal.threshold) {
      short.mat[, smp] <- ifelse(unique.genes %in%
                                   (dplyr::filter(dels, type == "del", sample == smp) %>% .$ensembl_gene_id),
                                 1, 0)
    } else {
      long.mat[, smp] <- ifelse(unique.genes %in%
                                  (dplyr::filter(dels, type == "del", sample == smp) %>% .$ensembl_gene_id),
                                1, 0)
    }
    
    if (total.loh.len < focal.threshold) {
      short.loh.mat[, smp] <- ifelse(unique.genes %in%
                                       (dplyr::filter(dels, type == "loh", sample == smp) %>% .$ensembl_gene_id),
                                     1, 0)
    } else {
      long.loh.mat[, smp] <- ifelse(unique.genes %in%
                                      (dplyr::filter(dels, type == "loh", sample == smp) %>% .$ensembl_gene_id),
                                    1, 0)
    }
  }
  
  # Combine focal and broad matrices
  all.short.mat <- short.mat + short.loh.mat
  all.long.mat <- long.mat + long.loh.mat
  
  # Compute probabilities for broad events and conditional focal events
  long.cnv.probabilities <- rowMeans(all.long.mat)
  no.long.cnv.probabilities <- 1 - long.cnv.probabilities
  short.given.long.cnv.probabilities <- rowMeans((all.short.mat > 0) & (all.long.mat > 0)) / long.cnv.probabilities
  short.given.no.long.cnv.probabilities <- rowMeans((all.short.mat > 0) & (all.long.mat == 0)) / no.long.cnv.probabilities
  
  # Compute Fisher's exact test statistics for co-deletion enrichment per gene
  message("Computing Fisher's exact test statistics for co-deletion enrichment...")
  del.fisher.df <- data.frame(
    foreach(i = 1:nrow(all.short.mat), .init = c(), .combine = rbind) %do% {
      long.samples <- colnames(all.long.mat)[all.long.mat[i, ] > 0]
      tbl <- data.frame(
        gene = c(focal = length(which(all.short.mat[i, long.samples] > 0)),
                 no.focal = length(which(all.short.mat[i, long.samples] == 0))),
        background = c(focal = length(which(all.short.mat[-i, long.samples] > 0)),
                       no.mut = length(which(all.short.mat[-i, long.samples] == 0)))
      )
      res <- fisher.test(tbl + 1, alternative = "greater")
      c(name = rownames(all.short.mat)[i], pval = res$p.value, odds.ratio = unname(res$estimate))
    }
  )
  del.fisher.df$qval <- p.adjust(del.fisher.df$pval, method = "BH")
  
  message("CNV matrices successfully computed. Preparing to return results.")
  
  # Return results as a list of matrices and data frames
  return(list(
    short.mat = short.mat,
    long.mat = long.mat,
    short.loh.mat = short.loh.mat,
    long.loh.mat = long.loh.mat,
    short.given.long.cnv.probabilities = data.frame(prob = short.given.long.cnv.probabilities),
    short.given.no.long.cnv.probabilities = data.frame(prob = short.given.no.long.cnv.probabilities),
    del.fisher.df = del.fisher.df,
    short.and.long = data.frame(count = rowSums((all.short.mat > 0) & (all.long.mat > 0)))
  ))
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
  make_option(c("-s", "--seg-file-path"),
              type = "character",
              default = NULL,
              help = "Path to the segmentation file",
              metavar = "character"),
  make_option(c("-p", "--sample-cutoffs-path"),
              type = "character",
              default = NULL,
              help = "Path to the sample cutoffs file",
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

# Check that all mandatory options are provided -----------------------------
if (any(unlist(lapply(option_list, function(option) {
  is.null(opt[[option@dest]])
})))) {
  print_help(opt_parser)
  stop("Error: All options are mandatory", call. = FALSE)
}

message("Starting CNV matrix generation process...")
# Generate CNV matrices using the provided parameters -------------------------
res <- generate_cnv_matrices(opt$`cancer-type`, opt$arm, opt$`seg-file-path`, opt$`sample-cutoffs-path`)

# Ensure the output directory exists -----------------------------------------
if (!dir.exists(opt$out)) {
  message("Creating output directory: ", opt$out)
  dir.create(opt$out, recursive = TRUE)
}

# Write results to CSV files in the output directory -------------------------
message("Writing focal CNVs matrix to CSV...")
write.csv(res$short.mat, here::here(opt$out, "focal_cnvs_matrix.csv"), row.names = TRUE)
message("Writing broad CNVs matrix to CSV...")
write.csv(res$long.mat, here::here(opt$out, "broad_cnvs_matrix.csv"), row.names = TRUE)
message("Writing focal LOH matrix to CSV...")
write.csv(res$short.loh.mat, here::here(opt$out, "focal_loh_cnvs_matrix.csv"), row.names = TRUE)
message("Writing broad LOH matrix to CSV...")
write.csv(res$long.loh.mat, here::here(opt$out, "broad_loh_cnvs_matrix.csv"), row.names = TRUE)
message("Writing focal given broad probabilities to CSV...")
write.csv(res$short.given.long.cnv.probabilities, here::here(opt$out, "focal_given_broad_probabilities.csv"), row.names = TRUE)
message("Writing focal given no broad probabilities to CSV...")
write.csv(res$short.given.no.long.cnv.probabilities, here::here(opt$out, "focal_given_no_broad_probabilities.csv"), row.names = TRUE)
message("Writing co-deletion Fisher results to CSV...")
write.csv(res$del.fisher.df, here::here(opt$out, "co_deletion_fisher.csv"), row.names = FALSE)
message("Writing focal and broad counts to CSV...")
write.csv(res$short.and.long, here::here(opt$out, "focal_and_broad.csv"), row.names = TRUE)

message("CNV matrices and analysis results successfully written to folder: ", opt$out)
