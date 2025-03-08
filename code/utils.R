load_ANP_data <- function(csv.path)
{
  tcga.anp.data<-read.csv(csv.path)
  
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