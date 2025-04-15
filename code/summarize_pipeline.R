library(here)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)
library(foreach)
library(PRODIGY)
library(GenomicRanges)

source(here::here('code/utils.R'))
source(here::here('code/genome_annotation_utils.R'))

load_snakemake_results <- function(){
  save.to = here::here(load_configuration()$out_folder, 'snakemake_results.csv')
  if(file.exists(save.to)){
    return(read.csv(save.to))
  }
  arm_ct_pairs = load_arm_cancer_type_pairs()
  all.genes <- 
    foreach(arm_ = unique(arm_ct_pairs$arm), .combine=rbind.fill, .init= data.frame()) %do%{
      foreach(cancer.type =(unique(dplyr::filter(arm_ct_pairs,arm == arm_) %>% .$cancer.type)), .init=data.frame(),.combine = rbind.fill) %do% {
          genes <- read.csv(here::here(load_configuration()$out_folder,arm_,cancer.type,'gene_features.csv'))
          genes$cancer.type = cancer.type
          genes$arm = arm_
          genes$mutation.depletion.qvalue <- p.adjust(genes$mutation.depletion.pvalue, method = 'BH')
          genes$mutation.fisher.qval <- p.adjust(genes$mutation.fisher.pval, method = 'BH')
          genes <- genes %>% replace(is.na(.), 0)
          genes
      }
    }
  write.csv(all.genes,save.to)
  return(all.genes)
}

load_single_gistic_results <- function(res.file){
  res <- read.table(res.file,header = FALSE,sep = '\t')
  res <- data.frame(t(res[1:4,-1]))
  rownames(res) <- NULL
  colnames(res) <- c('cytoband','qvalue','residual_qvalue','wide_peak_boundaries')
  res = res %>% dplyr::filter(! is.na(cytoband)) %>% 
    dplyr::mutate(qvalue = as.numeric(qvalue),
                  residual_qvalue = as.numeric(residual_qvalue),
                  arm = paste0('X',cytoband_to_arm(cytoband)),
                  start = as.numeric(str_extract(wide_peak_boundaries,'([\\w|\\d]+):([\\w|\\d]+)-([\\w|\\d]+)',group = 2)),
                  end = as.numeric(str_extract(wide_peak_boundaries,'([\\w|\\d]+):([\\w|\\d]+)-([\\w|\\d]+)',group = 3)))
  return(res)
}

load_all_gistic_results <- function(){
  result_foldr <- load_configuration()$out_folder
  given.arm.loss.path <- here::here(result_foldr, 'gistic_regions_given_arm_loss.csv')
  given.no.arm.loss.path <- here::here(result_foldr, 'gistic_regions_given_no_arm_loss.csv')
  if(file.exists(given.arm.loss.path) & file.exists(given.no.arm.loss.path)){
    given.arm.loss.resutls <- read.csv(given.arm.loss.path)
    given.no.arm.loss.results <- read.csv(given.no.arm.loss.path)
  } else {
    arm_ct_pairs = load_arm_cancer_type_pairs()
    
    given.arm.loss.resutls <- data.frame()
    given.no.arm.loss.results <- data.frame()
    for(arm_ in unique(arm_ct_pairs$arm)){
      for(cancer.type_ in unique((dplyr::filter(arm_ct_pairs, arm == arm_))$cancer.type)){
        gistic.res.path <- here::here(result_foldr, arm_, cancer.type_, 'GISTIC','aneuploid','del_genes.conf_99.txt')
        if(! file.exists(gistic.res.path)){
         print(paste0('File ',gistic.res.path,' does not exist. Skipping')) 
        } else {
          gistic.res <- load_single_gistic_results(gistic.res.path)
          gistic.res <- gistic.res %>% dplyr::filter(arm == arm_)
          given.arm.loss.resutls <- rbind(given.arm.loss.resutls, gistic.res %>% dplyr::mutate(cancer.type = cancer.type_))
        }
        
        gistic.res.path <- here::here(result_foldr, arm_, cancer.type_, 'GISTIC','non_aneuploid','del_genes.conf_99.txt')
        if(! file.exists(gistic.res.path)){
          print(paste0('File ',gistic.res.path,' does not exist. Skipping')) 
        } else{
          gistic.res <- load_single_gistic_results(gistic.res.path)
          gistic.res <- gistic.res %>% dplyr::filter(arm == arm_)
          given.no.arm.loss.results <- rbind(given.no.arm.loss.results, gistic.res %>% dplyr::mutate(cancer.type = cancer.type_))
        }
      }
    }
    write.csv(given.arm.loss.resutls, given.arm.loss.path, row.names = F)
    write.csv(given.no.arm.loss.results, given.no.arm.loss.path, row.names = F)
  }

  return(list(given.arm.loss = given.arm.loss.resutls,
              given.no.arm.loss = given.no.arm.loss.results))
}

filter_gistic_regions <- function(gistic.regions, regions.to.exclude) {
  gistic.regions$index <- 1:nrow(gistic.regions)
  
  chromosome_info <- get_chromosome_info_hg19()
  arm_sizes <- data.frame()
  for(i in seq(22)){
    arm_sizes <- rbind(arm_sizes,
                       data.frame(chrom = i,arm = paste0('X',i,'p'), size = chromosome_info$centromerStart[i]))
    arm_sizes <- rbind(arm_sizes,
                       data.frame(chrom = i,arm = paste0('X',i,'q'), size = (chromosome_info$length[i] - chromosome_info$centromerEnd[i])))
  }
  
  gistic.regions.fillterd <- gistic.regions %>% dplyr::mutate(chrom =  arm_to_chromosome(arm)) %>%
    merge(chromosome_info, by='chrom') %>%
    dplyr::filter(!((start < centromerStart) & (end > centromerStart)))
  
  gistic.regions.fillterd <- gistic.regions.fillterd %>% merge(arm_sizes, by='arm') %>% dplyr::mutate(del_size = end - start) %>%
    dplyr::mutate(fraction = del_size/size) %>% dplyr::filter(fraction < 0.9)
  
  
  regions = makeGRangesFromDataFrame(dplyr::select(gistic.regions.fillterd,seqnames = arm, start,end,cancer.type),keep.extra.columns=T)
  
  regions.to.exclude = makeGRangesFromDataFrame(dplyr::select(regions.to.exclude,seqnames = arm, start,end,cancer.type),keep.extra.columns=T)
  
  all_overlaps <- data.frame()
  for(arm_ in unique(gistic.regions.fillterd$arm)){
    for(cancer.type_ in unique(gistic.regions.fillterd %>% dplyr::filter(arm == arm_) %>% .$cancer.type)){
      regions_ = regions %>% dplyr::filter(seqnames == arm_,cancer.type == cancer.type_)
      regions.to.exclude_ = regions.to.exclude %>% dplyr::filter(seqnames == arm_,cancer.type == cancer.type_)
      overlaps = data.frame(findOverlaps(query = regions_, subject = regions.to.exclude_))
      if(nrow(overlaps) == 0){
        next
      }
      for(i in seq(nrow(overlaps))){
        all_overlaps <- rbind(all_overlaps,
                                data.frame(regions_[overlaps$queryHits[i],]) %>% 
                                  dplyr::select(arm = seqnames, cancer.type,
                                                start, end, width) %>% 
                                cbind(data.frame(overlap_size =
                                                   width(intersect(regions_[overlaps$queryHits[i],],regions.to.exclude_[overlaps$subjectHits[i],])))))
      }
    }
  }
  
  to_filter = all_overlaps %>%
    dplyr::mutate(r = overlap_size / width,
                  name = paste0(cancer.type,',',arm,':',start,':',end)) %>% 
    dplyr::filter(r > 0.5) %>% .$name
  
  gistic.regions.fillterd <- gistic.regions.fillterd %>% dplyr::mutate(name = paste0(cancer.type,',',arm,':',start,':',end)) %>%  dplyr::filter(! name %in% to_filter)
  return(gistic.regions.fillterd)
}

load_mono_cnv_genes <- function(){
  TSGs <- load_TSGs()
  oncos <- load_oncogenes()
  
  gistic.resutls <- load_all_gistic_results()
  regions.fillterd <- filter_gistic_regions(gistic.resutls$given.no.arm.loss, 
                                            gistic.resutls$given.arm.loss)
  
  all.prodigy.results <- data.frame()
  for(i in regions.fillterd$index){
    region <- regions.fillterd %>% dplyr::filter(index == i)
    file.path <- here::here(load_configuration()$out_folder,"PRODIGY",
                            paste0(c(i,region$cancer.type,region$arm,
                                     region$start,region$end
                                     ,"res.csv"), collapse = '_'))
    if((file.size(file.path) < 10 ) && (readLines(file.path) == '""')){
      next
    }
    influence.matrix <- read.csv(file.path)
    region.genes <- get_region_genes(region$start,
                                     region$end,
                                     arm_to_chromosome(region$arm))
    
    if(length(region.genes) == 0){
      next
    }
    
    region.influence.matrix <- dplyr::select(influence.matrix,any_of(region.genes))
    rownames(region.influence.matrix) <- seq(nrow(region.influence.matrix))
    pathways.to.take <- names(which(apply(region.influence.matrix, 
                                          1, function(x) length(which(x > 0))) < ncol(region.influence.matrix)/2))
    if((length(pathways.to.take) == 0) & (ncol(region.influence.matrix) >= 4)){
      next
    } else {
      for(j in seq(10)){
        anp.PRODIGY.res <- analyze_PRODIGY_results(list(as.matrix(region.influence.matrix)))
        if(length(anp.PRODIGY.res) == 0){
          next
        } else {
          if(length(anp.PRODIGY.res[[1]]) == 0){
            next
          }
          anp.PRODIGY.res <- anp.PRODIGY.res[[1]]
          rg <- data.frame(cytoband = region$cytoband,
                           arm = region$arm,
                           cancer.type = region$cancer.type,
                           PRODIGY.res = anp.PRODIGY.res,
                           rank = seq(length(anp.PRODIGY.res)),
                           idx = j
          )
          all.prodigy.results <- rbind(all.prodigy.results, rg)
        }
      }
    }
  }
  
  all.prodigy.results.summ <- all.prodigy.results %>%
    dplyr::group_by(cytoband,arm,cancer.type,PRODIGY.res) %>%
    dplyr::summarise(max_rank = max(rank)) %>% 
    dplyr::rename(name = PRODIGY.res)
  
  all.prodigy.results.summ <- all.prodigy.results.summ %>% 
    mutate(isTSG = (name %in% TSGs),
           isOncoKB = (name %in% oncos))
  all.prodigy.results.summ <- dplyr::filter(all.prodigy.results.summ,max_rank <= 5 ,isTSG)
  return(all.prodigy.results.summ)
}

load_co_cnv_genes <- function(){
  gistic.region.given.del <- load_all_gistic_results()$given.arm.loss
  
  for(i in seq(nrow(gistic.region.given.del))){
    region.genes <- get_region_genes(gistic.region.given.del[i,'start'],
                                     gistic.region.given.del[i,'end'],
                                     arm_to_chromosome(gistic.region.given.del[i,'arm']))
    gistic.region.given.del[i,'genes'] <- paste(region.genes, collapse = ',')
    gistic.region.given.del[i,'size'] <- length(region.genes)
  }
  
  cnv.given.del.stats <- load_snakemake_results() %>% 
    dplyr::select(name, cancer.type, focal.and.broad, 
                  focal.given.broad.probabilities)
  
  
  genes <- gistic.region.given.del %>% 
    separate_rows(genes,sep=',') %>% 
    dplyr::rename(name = genes) %>% 
    merge(cnv.given.del.stats, by = c('name','cancer.type')) %>% 
    dplyr::filter(focal.and.broad > 5, 
                  focal.given.broad.probabilities > 0.05) %>%
    dplyr::group_by(cancer.type,arm,start,end) %>%
    dplyr::summarise(name = paste(name, collapse = ','),
                     size = n()) %>% 
    separate_rows(name,sep=',')
  
  genes.filterd <- merge(genes,cnv.given.del.stats, by = c('name','cancer.type')) %>%
    dplyr::group_by(cancer.type,arm,start,end) %>%
    dplyr::filter(focal.given.broad.probabilities == max(focal.given.broad.probabilities)) %>%
    dplyr::filter(size <= 3)
  
  genes.filterd <- genes %>% dplyr::filter(name %in% genes.filterd$name)
  
  return(genes.filterd)
}

load_co_mutation_results <- function(){
  # if(is.null(save.to)){
  #   save.to = here::here(OUTPUT.DIR,'aneuploidy_drivers','co_occurring_mutations.csv')
  # }
  # if(file.exists(save.to)){
  #   return(read.csv(save.to))
  # }
  
  arm_ct_pairs = load_arm_cancer_type_pairs()
  TSGs <- load_TSGs()
  oncos <- load_oncogenes()
  
  all.mut.res <- data.frame()
  for(cancer.type_ in unique(arm_ct_pairs$cancer.type)){
    for(arm_ in unique(dplyr::filter(arm_ct_pairs, cancer.type == cancer.type_)$arm)){
      arm.genes <- get_arm_genes(arm_) %>% .$external_gene_name %>% unique()
      gistic.res.path <- here::here(load_configuration()$out_folder, arm_, cancer.type_, 'MutSig','aneuploid_mutations','sig_genes.txt')
      mutsig.res <- read.table(here::here(load_configuration()$out_folder, arm_, cancer.type_, 'MutSig','aneuploid_mutations','sig_genes.txt')
                               , sep = '\t', header = T) %>%
        dplyr::filter(gene %in% arm.genes) %>% mutate(cancer.type = cancer.type_, arm = arm_)
      
      all.mut.res <- rbind(all.mut.res, mutsig.res)
    }
  }
  all.mut.res.filterd <- dplyr::filter(all.mut.res, q < 0.25) %>% dplyr::select(name = gene,cancer.type,arm,q)
  all.mut.res.filterd <- dplyr::mutate(all.mut.res.filterd, 
                                       isTSG = (name %in% TSGs),
                                       isOncoKB = (name %in% oncos))
  # write.csv(all.mut.res.filterd, save.to, row.names = F)
  return(all.mut.res.filterd)
}
