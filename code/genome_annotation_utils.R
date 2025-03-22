library(stringr)
library(dplyr)
library(biomaRt)
library(rCGH)
library(GenomicRanges)
library(plyranges)

library(here)
source(here::here('code/utils.R'))


get_genes_form_biomaRt <- function(GRCh.version = '37'){
  if((GRCh.version != '37') & (GRCh.version != '38')){
    warning('GRCh version not supported')
    return(NULL)
  }
  
  GRCh.version.name <- paste0('grch',GRCh.version)
  genes.rds.path <- here::here(load_configuration()$data_folder,paste0(GRCh.version.name,'.rds'))
  if(file.exists(genes.rds.path)){
    return(readRDS(genes.rds.path))
  }
  else{
    if(GRCh.version == '37'){
      mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",host="https://grch37.ensembl.org",path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
    } else{
      mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",path="/biomart/martservice",dataset="hsapiens_gene_ensembl")
    }
    genes <- biomaRt::getBM(attributes = c("entrezgene_id","ensembl_gene_id","external_gene_name", "chromosome_name","transcript_biotype","start_position","end_position","band"), 
                            mart = mart,useCache = F)
    saveRDS(genes,genes.rds.path)
    return(genes)
  }
}

map_gene_identifiers <- function(ids, from, to, GRCh.version = '37', only_protein_coding = F){
  df <- data.frame(col1 = ids,idx=seq(length(ids)))
  colnames(df)[1] <- from
  genes <- get_genes_form_biomaRt(GRCh.version)
  if(only_protein_coding){
    genes <- filter(genes, transcript_biotype=='protein_coding')
  }
  genes <- filter(genes,!is.na(.data[[to]]))
  res <- (merge(df,genes,by=from,all.x=T) %>% filter(!duplicated(.data[[from]])) %>% arrange(idx))[,to]
  res[is.na(res)] <- ids[is.na(res)]
  return(res)
}

get_arm_genes <- function(arm, GRCh.version='37'){
  chr <- substr(arm ,2,nchar(arm)-1)
  arm <- substr(arm ,nchar(arm),nchar(arm))
  genes <- get_genes_form_biomaRt(GRCh.version)
  return(filter(genes, chromosome_name == chr, substr(band,1,1) == arm))
}

arm_to_chromosome <- function(arm){
  return(str_match(arm,'^X(\\d+|X|Y)(p|q)')[,2])
}

get_chromosome_info_hg19 <- function(){
  chromosomes.rds.path <- here::here(load_configuration()$data_folder,'grch19_chromosomes.rds')
  if(file.exists(chromosomes.rds.path)){
    return(readRDS(chromosomes.rds.path))
  }
  chromosome.info.df <- rCGH::hg19
  saveRDS(chromosome.info.df, chromosomes.rds.path)
  return(chromosome.info.df)
}

cytoband_to_arm <- function(cytobands){
  return(str_match(cytobands,'^(\\d+|X|Y)(p|q)')[,1])
}

get_region_genes <- function(start, end, chromosome, GRCh.version='37'){
  gene.info <- GRanges(get_genes_form_biomaRt(GRCh.version) %>% 
                         filter(chromosome_name %in% seq(1:22)) %>%
                         dplyr::select(Start=start_position,End=end_position,
                                       chromosome_name,external_gene_name, band) %>%
                         unique())
  gene.info <- gene.info[!duplicated(gene.info$external_gene_name),]
  ranges <- data.frame(start, end, chromosome) %>%
    GRanges %>% 
    plyranges::join_overlap_inner(gene.info)
  return(unique(ranges$external_gene_name))
}