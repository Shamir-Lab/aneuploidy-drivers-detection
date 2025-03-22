suppressPackageStartupMessages({
  library(here)
  library(optparse)
  library(TCGAbiolinks)
})

download_gene_expression <- function(cancer.type){
  project <- paste0('TCGA-',cancer.type)
  query.exp.hg38 <- GDCquery(project = project, 
                             data.category = "Transcriptome Profiling", 
                             data.type = "Gene Expression Quantification", 
                             workflow.type = "STAR - Counts")
  GDCdownload(query.exp.hg38)
  raw.counts <- GDCprepare(query = query.exp.hg38, save = FALSE)
  return(raw.counts)
}


option_list = list(
  make_option(c("-c", "--cancer-type"),
              type = "character",
              default = NULL,
              help = "cancer type",
              metavar = "character"),
  make_option(c("-o", "--out"),
              type = "character",
              default = NULL,
              help = "output folder path",
              metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
if (any(unlist(lapply( option_list, function(option){is.null(opt[[option@dest]])})))){
  print_help(opt_parser)
  stop("All options are mandatory", call. = FALSE)
}

res <- download_gene_expression(opt$`cancer-type`)

if(!dir.exists(opt$out)){
  dir.create(opt$out, recursive = TRUE)
}

saveRDS(res,here::here(opt$out,paste0(opt$`cancer-type`,'_gene_expression.rds')))