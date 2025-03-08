#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# Script: analyse_region.R
# Description:
#   This script performs region-specific analysis by integrating GISTIC results,
#   HiSeq expression data, and segmentation data. It:
#     1. Extracts differential expression (DEG) results for a specified genomic region
#        by comparing samples with deletions versus those without using DESeq2.
#     2. Uses the filtered expression matrix, together with pathway information,
#        to run the PRODIGY algorithm for pathway analysis.
#
# Usage:
#   Rscript analyse_region.R --cancer-type <cancer_type> --arm <arm> 
#         --region-start <start_coordinate> --region-end <end_coordinate> 
#         --hiSeq-data-path <path_to_HiSeq_RDS> --seg-path <path_to_segmentation_file>
#         --pathway-list-path <path_to_pathway_list_RDS> --out <output_folder>
#
# Example:
#   Rscript analyse_region.R --cancer-type BRCA --arm X3p --region-start 1000000 --region-end 5000000 \
#         --hiSeq-data-path data/BRCA_HiSeq.rds --seg-path data/BRCA.seg \
#         --pathway-list-path data/pathway_list.rds --out results/
#
# Output:
#   A CSV file containing the PRODIGY results for the specified genomic region will be
#   written to the output folder.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(here)
  set_here(path = normalizePath('..'))
  source(here::here('utils.R'))
  source(here::here('genome_annotation_utils.R'))
  
  library(PRODIGY)
  library(ff)
  library(optparse)
  library(dplyr)
  library(stringr)
  library(GenomicRanges)
  library(DESeq2)
})

# ------------------------------------------------------------------------------
# Function: get_deletion_DEG
# Purpose: Given HiSeq expression data and segmentation data along with a 
#          genomic region (as a GRanges object), filter the HiSeq data to only 
#          include primary tumors and compute differential expression (DEG) between 
#          samples with and without a deletion in that region using DESeq2.
#
# Parameters:
#   - HiSeq.data: A SummarizedExperiment object containing expression data.
#   - seg.data: A data frame with segmentation data.
#   - region: A GRanges object defining the genomic region of interest.
#
# Returns:
#   A DESeq2 results object (data frame) containing differential expression statistics.
# ------------------------------------------------------------------------------
get_deletion_DEG <- function(HiSeq.data, seg.data, region) {
  # Standardize sample names to 15 characters
  HiSeq.data@colData$sample <- str_sub(HiSeq.data@colData$sample, 1, 15)
  
  # Select primary tumor samples
  cancer.sample.ids <- subset(HiSeq.data, select = colData(HiSeq.data)$sample_type == 'Primary Tumor')$sample
  
  # Convert segmentation data to GRanges
  seg.data <- GRanges(seg.data)
  
  # Identify samples with deletion in the region (cutoff: Segment_Mean < -0.1)
  samples_with_deletion <- mergeByOverlaps(region, seg.data) %>% 
    data.frame() %>% 
    dplyr::filter(Segment_Mean < -0.1) %>% 
    .$seg.data.Sample %>% 
    unique() %>% 
    str_sub(1, 15)
  
  samples_with_deletion <- intersect(samples_with_deletion, cancer.sample.ids)
  samples_without_deletion <- setdiff(cancer.sample.ids, samples_with_deletion)
  samples_without_deletion <- intersect(samples_without_deletion, str_sub(seg.data$Sample, 1, 15))
  
  # Filter HiSeq data for relevant samples
  HiSeq.data.filtered <- subset(HiSeq.data, select = (colData(HiSeq.data)$sample %in% c(samples_with_deletion, samples_without_deletion)))
  expression.matrix <- assay(HiSeq.data.filtered, "raw_count")
  
  # Ensure unique gene identifiers
  expression.matrix <- expression.matrix[!duplicated(HiSeq.data@rowRanges$ensembl_gene_id), ]
  rownames(expression.matrix) <- HiSeq.data@rowRanges$ensembl_gene_id[!duplicated(HiSeq.data@rowRanges$ensembl_gene_id)]
  
  # Annotate samples: group2 for samples with deletion, group1 for others
  sample.annotation.vector <- ifelse(colData(HiSeq.data.filtered)$sample %in% samples_with_deletion, 'group2', 'group1')
  group <- factor(sample.annotation.vector)
  
  ddsMat <- DESeqDataSetFromMatrix(countData = round(expression.matrix),
                                   colData = DataFrame(group),
                                   design = ~ group)
  
  deg_result <- results(DESeq(ddsMat))
  return(deg_result)
}

# ------------------------------------------------------------------------------
# Function: analyse_region
# Purpose: Analyze a genomic region defined by user parameters (cancer type, arm,
#          region start, and region end) by:
#          1. Constructing a genomic region as a GRanges object.
#          2. Calling get_deletion_DEG() to compute differential expression (DEG).
#          3. Creating an expression matrix from the HiSeq data.
#          4. Filtering DEGs based on significance and fold change.
#          5. Running the PRODIGY algorithm for pathway analysis.
#
# Parameters:
#   - cancer.type: Character. The cancer type (e.g., "BRCA").
#   - arm: Character. The chromosome arm (e.g., "X3p").
#   - region_start: Numeric. Start coordinate of the region.
#   - region_end: Numeric. End coordinate of the region.
#   - HiSeq.data: A SummarizedExperiment object with HiSeq expression data.
#   - seg.data: A data frame with segmentation data.
#   - pathway.list: A list or object defining pathways for PRODIGY analysis.
#
# Returns:
#   The output of the PRODIGY algorithm (typically a list or data frame with pathway analysis results).
# ------------------------------------------------------------------------------
analyse_region <- function(cancer.type, arm, region_start, region_end, HiSeq.data, seg.data, pathway.list) {
  message("Constructing genomic region for analysis...")
  # Construct the genomic region from provided parameters
  region <- GRanges(seqnames = arm_to_chromosome(arm), 
                    ranges = IRanges(start = region_start, end = region_end))
  
  message("Calculating differential expression for the specified region...")
  # Compute differential expression for the region
  DEGs <- data.frame(get_deletion_DEG(HiSeq.data, seg.data, region))
  
  message("Preparing expression matrix...")
  # Create an expression matrix from HiSeq data (using all primary tumor samples)
  
  expression.matrix <- assay(HiSeq.data, "raw_count")
  genes <- HiSeq.data@rowRanges$gene_id
  expression.matrix <- data.frame(expression.matrix[!duplicated(genes), ])
  rownames(expression.matrix) <- genes[!duplicated(genes)]
  
  message("Filtering DEG results for significance...")
  # Filter DEG results: p-value < 0.05 and absolute log2FoldChange > 2
  DEGs <- dplyr::filter(DEGs, pvalue < 0.05, abs(log2FoldChange) > 2)
  diff_genes <- abs(DEGs$log2FoldChange)
  names(diff_genes) <- map_gene_identifiers(rownames(DEGs), 'ensembl_gene_id', 'external_gene_name')
  
  message("Running PRODIGY pathway analysis...")
  # Use a global STRING network (assumed to be defined in utils.R)
  network <- STRING_network
  
  # Get candidate genes in the region (function defined in genome_annotation_utils.R)
  candidate_genes <- get_region_genes(region_start, region_end, arm_to_chromosome(arm))
  
  res <- PRODIGY(mutated_genes = candidate_genes,
                 expression_matrix = expression.matrix, network = network, sample = NULL,
                 diff_genes = diff_genes, alpha = 0.05, pathway_list  = pathway.list,
                 write_results = FALSE, beta = 2, gamma = 0.05, delta = 0.05,
                 num_of_cores = 1)
  
  return(res)
}

# ------------------------------------------------------------------------------
# Command-Line Interface and Script Execution
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("-c", "--cancer-type"),
              type = "character",
              default = NULL,
              help = "Cancer type (e.g., BRCA)",
              metavar = "character"),
  make_option(c("-a", "--arm"),
              type = "character",
              default = NULL,
              help = "Chromosome arm (e.g., X3p)",
              metavar = "character"),
  make_option(c("--region-start"),
              type = "numeric",
              default = NULL,
              help = "Start coordinate of the region",
              metavar = "number"),
  make_option(c("--region-end"),
              type = "numeric",
              default = NULL,
              help = "End coordinate of the region",
              metavar = "number"),
  make_option(c("--hiSeq-data-path"),
              type = "character",
              default = NULL,
              help = "Path to the HiSeq data RDS file (SummarizedExperiment object)",
              metavar = "character"),
  make_option(c("--seg-path"),
              type = "character",
              default = NULL,
              help = "Path to the segmentation data file (text file with headers)",
              metavar = "character"),
  make_option(c("--pathway-list-path"),
              type = "character",
              default = NULL,
              help = "Path to the pathway list file (RDS)",
              metavar = "character"),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "Output file path for saving the PRODIGY results CSV file",
              metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Verify that all mandatory options are provided
if (any(unlist(lapply(option_list, function(option) { is.null(opt[[option@dest]]) })))) {
  print_help(opt_parser)
  stop("Error: All options are mandatory", call. = FALSE)
}

message("Loading HiSeq data from: ", opt$`hiSeq-data-path`)
HiSeq.data <- readRDS(opt$`hiSeq-data-path`)

message("Loading segmentation data from: ", opt$`seg-path`)
seg.data <- read.table(opt$`seg-path`, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

message("Loading pathway list from: ", opt$`pathway-list-path`)
pathway.list <- readRDS(opt$`pathway-list-path`)

# Convert the provided index to numeric
index <- as.numeric(opt$index)

# Run analysis on the specified region
message("Analyzing region with index: ", index)
res <- analyse_region(opt$`cancer-type`, opt$arm, opt$`region-start`, opt$`region-end`, HiSeq.data, seg.data, pathway.list)

# Create the output directory if it does not exist
if (!dir.exists(dirname(opt$out))) {
  message("Creating output directory: ", dirname(opt$out))
  dir.create(dirname(opt$out), recursive = TRUE)
}

# Write the PRODIGY results to a CSV file
message("Writing analysis results to: ", opt$out)
write.csv(as.data.frame(res), opt$out, row.names = TRUE)

message("Analysis complete. Results saved to: ", opt$out)
