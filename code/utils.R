library(yaml)
library(here)

load_configuration <- function(){
  return(read_yaml(here::here('config.yaml'))) 
}

load_ANP_data <- function() {
  csv.path <- here::here(load_configuration()$data_folder,'Taylor_et_al._Arm-Level_WGD_TCGA_data.csv')
  tcga.anp.data <- read.csv(csv.path)
  
  #rename 'i..Sample' column to 'sample' and add 'patient' column for joining with subtypes
  names(tcga.anp.data)[1] <-"sample"
  tcga.anp.data["patient"] <- tcga.anp.data['sample']
  
  #rename 'AneuploidyScore.AS.' column
  names(tcga.anp.data)[names(tcga.anp.data) == 'AneuploidyScore.AS.'] <- 'AneuploidyScore'
  
  
  #rename 'patient' column removing '-XX' suffix 
  tcga.anp.data$patient<-substr(tcga.anp.data$patient,start=0,stop=12)
  
  #get column names for chromosome arms
  fisrt.column.idx<- which(names(tcga.anp.data)=='X1p')
  last.column.idx<- which(names(tcga.anp.data)=='X22q')
  chromosome.column.names <-  names(tcga.anp.data)[fisrt.column.idx:last.column.idx]
  
  tcga.anp.data[,chromosome.column.names][is.na(tcga.anp.data[,chromosome.column.names])]<-0
  return(list(tcga.anp.data=tcga.anp.data,chromosome.column.names=chromosome.column.names))
}

load_arm_cancer_type_pairs <- function(){
  csv.path <- here::here(load_configuration()$data_folder,'arm_cancer_type_paires.csv')
  return(read.csv(csv.path))
}

load_TSGs <- function(){
  all.tsgs <- read.csv(here::here(load_configuration()$data_folder,'Human_TSGs.txt'),header = T,sep = '\t')
  cancer.census.genes <- read.csv(here::here(load_configuration()$data_folder,'Census_allSat Oct 15 12_04_35 2022.csv') ,header = T)
  cancer.census.genes.tsgs <- cancer.census.genes[grep('TSG|^$',cancer.census.genes$Role.in.Cancer),]
  return(c(all.tsgs$GeneSymbol, cancer.census.genes.tsgs$Gene.Symbol) %>% unique())
}

load_oncogenes <- function(){
  oncokb_data <- read.csv(here::here(load_configuration()$data_folder,'cancerGeneList_oncokb.tsv'), sep='\t')
  return(dplyr::filter(oncokb_data, `Is.Oncogene`=='Yes',`Is.Tumor.Suppressor.Gene`!='Yes') %>%
           .$`Hugo.Symbol`)
}