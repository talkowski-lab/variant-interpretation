#!/usr/bin/env R

#This script runs C-Alpha on a pair of phenotypes
#Author: Alba Sanchis Juan, based on Jack Fu's scripts

#############
##LIBRARIES##
#############
library(data.table)
library(optparse)
library(tidyverse)
library(parallel)

#############
##ARGUMENTS##
#############
option_list = list(
  make_option(c("-s", "--input_snv"), type="character", default=NULL, 
              help="input file name", metavar="character"),
  make_option(c("-p", "--pedigree"), type="character", default=NULL, 
              help="pedigree file", metavar="character"),
  make_option(c("-g", "--genes"), type="character", default=".", 
              help="list of genes to run c-alpha on [default= %default]", metavar="character"),
  make_option(c("-a", "--phenotype1"), type="character", default=".", 
              help="phenotype1 [default= %default]", metavar="character"),
  make_option(c("-b", "--phenotype2"), type="character", default=".", 
              help="phenotype2 [default= %default]", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=".", 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-c", "--column_name"), type="character", default=".", 
              help="column name to run c-alpha on [default= %default]", metavar="character"),
  make_option(c("-u", "--subsample"), type="character", default=".", 
              help="column name to run c-alpha on [default= %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt <- parse_args(opt_parser)

dv_total <- fread(opt$input_snv)
pedigree_final <- fread(opt$pedigree)
genes <- fread(opt$genes, header = F)

outfile <- opt$output
column_name <- opt$column_name

Group1 <- opt$phenotype1
Group2 <- opt$phenotype2

do_subsample <- opt$subsample

#############
##FUNCTIONS##
#############
calpha_method_newer <- function(variant_list_with_phenotypes, p0, singletons) {
  # ASSUMES ALL VARIANTS ARE HETS AND THERE IS A BINARY "PHENO" COLUMN
  
  # Variant counts by gene and phenotype
  gene_table = group_by(variant_list_with_phenotypes, Gene)
  gene_table = summarise(gene_table, case_alleles = sum(Pheno), total_alleles = sum(nNonRef))
  
  if(singletons == 'group') {
    singles = subset(gene_table, total_alleles == 1)
    gene_table = subset(gene_table, total_alleles > 1)
    
    gene_table = rbind(gene_table, data.frame(Gene = 'SINGLETONS', 
                                              case_alleles = sum(singles$case_alleles), 
                                              total_alleles = sum(singles$total_alleles)))
  }
  
  # Test statistic
  Talpha = sum((gene_table$case_alleles - p0*gene_table$total_alleles)^2 - (gene_table$total_alleles * p0 * (1-p0)))
  
  return(Talpha)
}

calpha_function_newer <- function(variant_list_with_phenotypes, perm, singletons) {
  # ASSUMES ALL VARIANTS ARE HETS AND THERE IS A BINARY "PHENO" COLUMN
  
  start.time = Sys.time()
  
  if (!singletons %in% c('keep', 'group', 'drop')) { print ('No special handling for singletons by default') }
  
  # Variant counts by gene
  gene_table = table(variant_list_with_phenotypes$Gene)
  
  # Group singletons?
  if(singletons == 'group') {
    singles = gene_table[gene_table == 1]
    gene_table = gene_table[gene_table > 1]
    
    gene_table = c(gene_table, SINGLETONS = sum(singles))
  }
  
  # Drop singletons?
  if(singletons == 'drop') {
    gene_table = gene_table[gene_table > 1]
    
    # Adjust p0
    variant_list_with_phenotypes = subset(variant_list_with_phenotypes, Gene %in% rownames(gene_table))
  }
  
  # Get number of genes
  num_genes = length(gene_table)
  
  # Calculate p0...
  # How many alleles from cases
  nAff = sum(variant_list_with_phenotypes$Pheno)
  # How many total alleles
  nTot = nrow(variant_list_with_phenotypes)
  # Proportion of alleles from cases
  p0 = nAff / nTot  
  
  # Test statistic 
  calpha.stat = calpha_method_newer(variant_list_with_phenotypes, p0, singletons)
  
  # Variance of Talpha
  Valpha = 0
  for (i in 1:num_genes) {
    for (u in 0:gene_table[i]) {
      Valpha = Valpha + (((u - gene_table[i]*p0)^2 - gene_table[i]*p0*(1-p0))^2)*dbinom(u, gene_table[i], p0)
    }
  }
  names(Valpha) = NULL
  
  # Asymptotic p-value
  if (Valpha==0) { 
    asym.pval = 1 
  } else {
    asym.pval = 1 - pchisq(calpha.stat^2 / Valpha, df=1)
  }
  
  # Permutations
  perm.pval = NA
  
  #run permutations in parallel
  Pheno_original <- variant_list_with_phenotypes$Pheno
  n <- nrow(variant_list_with_phenotypes)
  
  x.perm <- unlist(
    mclapply(1:perm, function(i) {
      shuffled <- Pheno_original[sample.int(n)]
      variant_list_with_phenotypes$Pheno <- shuffled
      calpha_method_newer(variant_list_with_phenotypes, p0, singletons)
    }, mc.cores = 8)  # adjust number of cores
  )
  
  # 
  # x.perm = rep(0, perm)
  # for (i in 1:perm) {
  #   permuted_variant_list = variant_list_with_phenotypes
  #   permuted_variant_list$Pheno = sample(permuted_variant_list$Pheno)
  #   x.perm[i] = calpha_method_newer(permuted_variant_list, p0, singletons) 
  # }
  # p-value 
  perm.pval = sum(x.perm^2 > calpha.stat^2) / perm
  
  ## results
  name = "CALPHA: c-alpha Test"
  arg.spec = c(nAff, nTot-nAff, num_genes, perm)
  names(arg.spec) = c("case.alleles", "control.alleles", "genes", "n.perms")
  res = list(calpha.stat = calpha.stat, 
             asym.pval = asym.pval, 
             perm.pval = perm.pval, 
             args = arg.spec, 
             name = name,
             perms = x.perm)
  class(res) = "assoctest"  
  
  end.time = Sys.time()
  print(end.time - start.time)
  
  print(res)
  
  return(res)
}

runCalpha_1_vs_1_2 <- function(pedigree, group1, group2, genes, group1_level, group2_level, subsample){
  # #test
  # group1 <- "Abnormality of the nervous system"
  # group2 <- "Autism"
  # group2 <- "Abnormality of the eye"
  # group1_level <- "hpo_present_3th_ancest"
  # group2_level <- "hpo_present_name"
  # genes <- constrained_genes
  
  pedigree_aff <- subset(pedigree, affected == 2) ##keep affected probands only for now
  
  all_calls <- dv_total[dv_total$SYMBOL %in% genes,]
  
  if(subsample){
    group1_count <- nrow(subset(pedigree_aff, grepl(group1, pedigree_aff[[group1_level]])))
    group1_2_count <- nrow(subset(pedigree_aff, grepl(group1, pedigree_aff[[group1_level]]) &  grepl(group2, pedigree_aff[[group2_level]])))
    min_downsample <- min(group1_count, group1_2_count)
    
    pedigree_group1 <- subset(pedigree_aff, grepl(group1, pedigree_aff[[group1_level]]) & !grepl(group2, pedigree_aff[[group2_level]]))
    pedigree_group1_2 <- subset(pedigree_aff, grepl(group1, pedigree_aff[[group1_level]]) & grepl(group2, pedigree_aff[[group2_level]]))
    
    if(group1_count > group1_2_count){
      set.seed(1337)
      
      pedigree_group1_downsampled <- pedigree_group1 %>%
        slice_sample(n = min_downsample)
      
      pedigree_downsampled <- rbind(pedigree_group1_downsampled, pedigree_group1_2)
      
      all_calls$group1 <- ifelse(all_calls$snv_pipeline_id_batch %in% subset(pedigree_downsampled, grepl(group1, pedigree_downsampled[[group1_level]]))$snv_pipeline_id_batch, 1, 0)
      all_calls$group2 <- ifelse(all_calls$snv_pipeline_id_batch %in% subset(pedigree_downsampled, grepl(group2, pedigree_downsampled[[group2_level]]))$snv_pipeline_id_batch, 1, 0)
      
    }else{
      
      pedigree_group1_2_downsampled <- pedigree_group1_2 %>%
        slice_sample(n = min_downsample)
      
      pedigree_downsampled <- rbind(pedigree_group1, pedigree_group1_2_downsampled)
      
      all_calls$group1 <- ifelse(all_calls$snv_pipeline_id_batch %in% subset(pedigree_downsampled, grepl(group1, pedigree_downsampled[[group1_level]]))$snv_pipeline_id_batch, 1, 0)
      all_calls$group2 <- ifelse(all_calls$snv_pipeline_id_batch %in% subset(pedigree_downsampled, grepl(group2, pedigree_downsampled[[group2_level]]))$snv_pipeline_id_batch, 1, 0)
    }
    
  }else{
    all_calls$group1 <- ifelse(all_calls$snv_pipeline_id_batch %in% subset(pedigree_aff, grepl(group1, pedigree_aff[[group1_level]]))$snv_pipeline_id_batch, 1, 0)
    all_calls$group2 <- ifelse(all_calls$snv_pipeline_id_batch %in% subset(pedigree_aff, grepl(group2, pedigree_aff[[group2_level]]))$snv_pipeline_id_batch, 1, 0)
  }
  
  all_calls_group1 <- subset(all_calls, group1 ==1) #take only group1 to compare group1-group2 vs group1+group2
  # all_calls_NDD[all_calls_NDD$NDD_ASD == 1,]$NDD <- 0
  
  ###########
  ##C-ALPHA##
  ###########
  
  # Define permutations
  perms = 10000
  # perms = 100000
  
  ##Synonymous
  vars_syn <- subset(all_calls_group1, isSYN)
  vars_syn$Pheno <- vars_syn$group2
  vars_syn$nNonRef <- 1 #Kyle says "nNonRef is just the number of non-ref alleles in a given genotype call (and for the de novos I have been working with, it's always been 1)"
  # vars_syn$nNonRef <- 1 - vars_syn$Pheno
  
  set.seed(1337)
  calpha_syn <- calpha_function_newer(vars_syn, perms, singletons = 'keep')
  
  ##PTV
  vars_ptv <- subset(all_calls_group1, isCoding == TRUE & isPTV)
  vars_ptv$Pheno <- vars_ptv$group2
  vars_ptv$nNonRef <- 1 #Kyle says "nNonRef is just the number of non-ref alleles in a given genotype call (and for the de novos I have been working with, it's always been 1)"
  
  set.seed(1337)
  calpha_ptv <- calpha_function_newer(vars_ptv, perms, singletons = 'keep')
  
  ##MIS2
  vars_mis2 <- subset(all_calls_group1, isCoding == TRUE & (isMIS & MPC_v2 >= 2 & am_pathogenicity >= 0.97))
  vars_mis2$Pheno <- vars_mis2$group2
  vars_mis2$nNonRef <- 1
  
  set.seed(1337)
  calpha_mis2 <- calpha_function_newer(vars_mis2, perms, singletons = 'keep')
  
  ##DMG = PTV+MIS2
  # vars <- dv_total[which(dv_total$Gene %in% qval6$gene[which(qval6$FDR<0.001)] & (dv_total$isPTV | dv_total$isMIS2 | dv_total$isOS)),]
  vars_dmg <- subset(all_calls_group1, isCoding == TRUE & (isPTV | (isMIS & MPC_v2 >= 2 & am_pathogenicity >= 0.97)))
  vars_dmg$Pheno <- vars_dmg$group2
  vars_dmg$nNonRef <- 1
  
  set.seed(1337)
  
  calpha_dmg <- calpha_function_newer(vars_dmg, perms, singletons = 'keep')
  
  calpha_list <- list(
    syn  = calpha_syn,
    mis2 = calpha_mis2,
    ptv  = calpha_ptv,
    dmg  = calpha_dmg
  )
  
  if (subsample) {
    calpha_list$min_downsample <- min_downsample
  }
  
  calpha_list
}

runCalpha_1_vs_1 <- function(pedigree, group1, group2, genes, group1_level, group2_level, subsample){
  # #test
  # group1 <- "Abnormality of the nervous system"
  # group2 <- "Autism"
  # group1_level <- "hpo_present_3th_ancest"
  # group2_level <- "hpo_present_name"
  # genes <- constrained_genes
  
  pedigree_aff <- subset(pedigree, affected == 2) ##keep affected probands only for now
  
  all_calls <- dv_total[dv_total$SYMBOL %in% genes,]
  
  if(subsample){
    group1_count <- nrow(subset(pedigree_aff, grepl(group1, pedigree_aff[[group1_level]]) & ! grepl(group2, pedigree_aff[[group2_level]])))
    group2_count <- nrow(subset(pedigree_aff, !grepl(group1, pedigree_aff[[group1_level]]) &  grepl(group2, pedigree_aff[[group2_level]])))
    min_downsample <- min(group1_count, group2_count)
    
    pedigree_group1 <- subset(pedigree_aff, grepl(group1, pedigree_aff[[group1_level]]) & !grepl(group2, pedigree_aff[[group2_level]]))
    pedigree_group2 <- subset(pedigree_aff, !grepl(group1, pedigree_aff[[group1_level]]) & grepl(group2, pedigree_aff[[group2_level]]))
    
    if(group1_count > group2_count){
      set.seed(1337)

      pedigree_group1_downsampled <- pedigree_group1 %>%
        slice_sample(n = min_downsample)
      
      pedigree_downsampled <- rbind(pedigree_group1_downsampled, pedigree_group2)
      
      all_calls$group1 <- ifelse(all_calls$snv_pipeline_id_batch %in% subset(pedigree_downsampled, grepl(group1, pedigree_downsampled[[group1_level]]))$snv_pipeline_id_batch, 1, 0)
      all_calls$group2 <- ifelse(all_calls$snv_pipeline_id_batch %in% subset(pedigree_downsampled, grepl(group2, pedigree_downsampled[[group2_level]]))$snv_pipeline_id_batch, 1, 0)

    }else{
      set.seed(1337)
      
      pedigree_group2_downsampled <- pedigree_group2 %>%
        slice_sample(n = min_downsample)
      
      pedigree_downsampled <- rbind(pedigree_group2_downsampled, pedigree_group1)
      
      all_calls$group1 <- ifelse(all_calls$snv_pipeline_id_batch %in% subset(pedigree_downsampled, grepl(group1, pedigree_downsampled[[group1_level]]))$snv_pipeline_id_batch, 1, 0)
      all_calls$group2 <- ifelse(all_calls$snv_pipeline_id_batch %in% subset(pedigree_downsampled, grepl(group2, pedigree_downsampled[[group2_level]]))$snv_pipeline_id_batch, 1, 0)
      
    }
  }else{
    
    all_calls$group1 <- ifelse(all_calls$snv_pipeline_id_batch %in% subset(pedigree_aff, affected == 2 & grepl(group1, pedigree_aff[[group1_level]]))$snv_pipeline_id_batch, 1, 0)
    all_calls$group2 <- ifelse(all_calls$snv_pipeline_id_batch %in% subset(pedigree_aff, affected == 2 & grepl(group2, pedigree_aff[[group2_level]]))$snv_pipeline_id_batch, 1, 0)
    
  }
  
  all_calls_group1_2 <- subset(all_calls, (group1 ==1 & group2 == 0) | (group2 ==1 & group1 == 0)) #take only group1 to compare group1 vs group2
  all_calls_group1_2$Pheno <- ifelse(all_calls_group1_2$group1 == 1, 1, 0) #group1 = 1, group2 = 0
  all_calls_group1_2$nNonRef <- 1
  # all_calls_NDD[all_calls_NDD$NDD_ASD == 1,]$NDD <- 0
  
  ###########
  ##C-ALPHA##
  ###########
  
  # Define permutations
  perms = 10000
  # perms = 100000
  
  ##Synonymous
  vars_syn <- subset(all_calls_group1_2, isSYN)
  #vars_syn$Pheno <- vars_syn$group2
  #vars_syn$nNonRef <- 1 #Kyle says "nNonRef is just the number of non-ref alleles in a given genotype call (and for the de novos I have been working with, it's always been 1)"
  #vars_syn$nNonRef <- 1 - vars_syn$Pheno
  
  set.seed(1337)
  calpha_syn <- calpha_function_newer(vars_syn, perms, singletons = 'keep')
  
  ##PTV
  vars_ptv <- subset(all_calls_group1_2, isCoding == TRUE & isPTV)
  #vars_ptv$Pheno <- vars_ptv$group2
  #vars_ptv$nNonRef <- 1 #Kyle says "nNonRef is just the number of non-ref alleles in a given genotype call (and for the de novos I have been working with, it's always been 1)"
  
  set.seed(1337)
  calpha_ptv <- calpha_function_newer(vars_ptv, perms, singletons = 'keep')
  
  ##MIS2
  vars_mis2 <- subset(all_calls_group1_2, isCoding == TRUE & (isMIS & MPC_v2 >= 2 & am_pathogenicity >= 0.97))
  #vars_mis2$Pheno <- vars_mis2$group2
  #vars_mis2$nNonRef <- 1
  
  set.seed(1337)
  calpha_mis2 <- calpha_function_newer(vars_mis2, perms, singletons = 'keep')
  
  ##DMG = PTV+MIS2
  # vars <- dv_total[which(dv_total$Gene %in% qval6$gene[which(qval6$FDR<0.001)] & (dv_total$isPTV | dv_total$isMIS2 | dv_total$isOS)),]
  vars_dmg <- subset(all_calls_group1_2, isCoding == TRUE & (isPTV | (isMIS & MPC_v2 >= 2 & am_pathogenicity >= 0.97)))
  #vars_dmg$Pheno <- vars_dmg$group2
  #vars_dmg$nNonRef <- 1
  
  set.seed(1337)
  
  calpha_dmg <- calpha_function_newer(vars_dmg, perms, singletons = 'keep')
  
  calpha_list <- list(
    syn  = calpha_syn,
    mis2 = calpha_mis2,
    ptv  = calpha_ptv,
    dmg  = calpha_dmg
  )
  
  if (subsample) {
    calpha_list$min_downsample <- min_downsample
  }
  
  calpha_list
}

###################
###Run everything##
###################

constrained_genes <- subset(genes, !is.na(V1))$V1

##Run group1 vs group1+group2
#outfile <- gsub(" ", "_", sprintf(paste0(c_alpha_outdir, "/%s_%s.RData"), gsub("/", "_", row$Group1), gsub("/", "_", row$Group2)))
  
res_1_12 <- runCalpha_1_vs_1_2(
  pedigree = pedigree_final,
  group1 = Group1,
  group2 = Group2,
  genes = constrained_genes,
  group1_level = column_name,
  group2_level = column_name,
  subsample = do_subsample
)
  
##Run group2 vs group1+group2
res_2_12 <- runCalpha_1_vs_1_2(
  pedigree = pedigree_final,
  group1 = Group2,
  group2 = Group1,
  genes = constrained_genes,
  group1_level = column_name,
  group2_level = column_name,
  subsample = do_subsample
)
  
res_1_2 <- runCalpha_1_vs_1(
  pedigree = pedigree_final,
  group1 = Group1,
  group2 = Group2,
  genes = constrained_genes,
  group1_level = column_name,
  group2_level = column_name,
  subsample = do_subsample
)

save(res_1_12, res_2_12, res_1_2, Group1, Group2, file = outfile)
