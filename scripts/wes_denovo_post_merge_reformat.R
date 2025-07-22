#!/usr/bin/env Rscript
#Author: Alba Sanchis-Juan

##Update merged calls and reformat WES de novo callset

##Load libraries
library(dplyr)
library(tidyverse)
library(data.table)
library(optparse)

option_list = list(
  make_option(c("-d", "--denovo"), type="character", default=NULL,
              help="denovo callset", metavar="character"),
  make_option(c("-m", "--merged"), type="character", default=NULL,
              help="merged callset", metavar="character"),
  make_option(c("-f", "--flipbook"), type="character", default=NULL,
              help="flipbook responses", metavar="character"),
  make_option(c("-r", "--release"), type="character", default=NULL,
              help="release date", metavar="character"),
  make_option(c("-i", "--ids"), type="character", default=NULL,
              help="ids_corresp file", metavar="character"),
  make_option(c("-o", "--outputdir"), type="character", default=NULL,
              help="output file", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

#Define input parameteres and read files
denovo <- opt$denovo
denovo_merged <- opt$merged
flipbook <- opt$flipbook
ids_corresp <- opt$ids
release <- opt$release
output_dir <- opt$outputdir

##read files
wes_denovo <- fread(denovo)
wes_denovo_merged <- fread(denovo_merged)
names(wes_denovo_merged) <- c("chrom_type_sample", "start_merged", "end_merged", "name_merged")
wes_denovo_merged$name <- wes_denovo_merged$name_merged

wes_denovo_merged %>%
  separate_rows(name, sep = ",", convert = TRUE) -> wes_denovo_merged_ref

wes_denovo_name <- merge(wes_denovo, unique(subset(wes_denovo_merged_ref, select = c(name_merged, name))), by = "name", all.x = T, all.y = F)

##Add flipbook
if( length(flipbook) >0 & length(ids_corresp) >0){
  flipbook <- fread(flipbook)
  ids_corresp <- fread(ids_corresp)

  wes_denovos_path <- merge(wes_denovo_name, ids_corresp, by = "name", all.x = T, all.y = F)
  wes_denovos_flipbook <- merge(wes_denovos_path, flipbook, by = "Path", all.x = T, all.y = F)
}else{
  wes_denovos_flipbook <- wes_denovo_name
}
# wes_denovos_flipbook[wes_denovos_flipbook$cohort == "ASD",]$`Is de novo` <- "yes"

##Reformat
wes_denovos_flipbook %>%
  group_by(name_merged) %>%
  reframe("chrom" = unique(chr),
          "start" = min(start),
          "end" = max(end),
          "svtype" = unique(svtype),
          "maternal_id" = ifelse("maternal_id" %in% names(wes_denovos_flipbook), unique(maternal_id), "0"),
          "paternal_id" = ifelse("paternal_id" %in% names(wes_denovos_flipbook), unique(paternal_id), "0"),
          # "paternal_id" = unique(paternal_id),
          "sample" = unique(sample),
          "batch" = ifelse("batch" %in% names(wes_denovos_flipbook), unique(batch), "0"),
          "CN" = paste(unique(CN), collapse = ","),
          "NP" = sum(NP),
          "NEx" = ifelse("NEx" %in% names(wes_denovos_flipbook), sum(NEx), "NA"),
          # "NEx" = sum(NEx),
          'QA' = paste(unique(QA), collapse = ","),
          'QS' = paste(unique(QS), collapse = ","),
          'QSS' = paste(unique(QSS), collapse = ","),
          'QSE' = paste(unique(QSE), collapse = ","),
          'PASS_SAMPLE' = paste(unique(PASS_SAMPLE), collapse = ","),
          'PASS_QS' = paste(unique(PASS_QS), collapse = ","),
          'PASS_FREQ' = paste(unique(PASS_FREQ), collapse = ","),
          'HIGH_QUALITY' = paste(unique(HIGH_QUALITY), collapse = ","),
          'cohort' = ifelse("cohort" %in% names(wes_denovos_flipbook), paste(unique(cohort), collapse = ","), "NA"),
          'variant_name' = paste(unique(variant_name), collapse = ","),
          'ID' = paste(unique(ID), collapse = ","),
          'rmsstd' = ifelse("rmsstd" %in% names(wes_denovos_flipbook), paste(unique(rmsstd), collapse = ","), "NA"),
          'sc' = paste(unique(sc), collapse = ","),
          'sf' = paste(unique(sf), collapse = ","),
          'cov_p' = ifelse("cov_p" %in% names(wes_denovos_flipbook), paste(unique(cov_p), collapse = ","), "NA"),
          'cov_m' = ifelse("cov_m" %in% names(wes_denovos_flipbook), paste(unique(cov_m), collapse = ","), "NA"),
          'cov_p_bs' = ifelse("cov_p_bs" %in% names(wes_denovos_flipbook), paste(unique(cov_p_bs), collapse = ","), "NA"),
          'cov_m_bs' = ifelse("cov_m_bs" %in% names(wes_denovos_flipbook), paste(unique(cov_m_bs), collapse = ","), "NA"),
          'inheritance' = ifelse("inheritance" %in% names(wes_denovos_flipbook), paste(unique(inheritance), collapse = ","), "NA"),
          'family_id' = ifelse("family_id" %in% names(wes_denovos_flipbook), paste(unique(family_id), collapse = ","), "NA"),
          'sample_fix' = ifelse("sample_fix" %in% names(wes_denovos_flipbook), paste(unique(sample_fix), collapse = ","), "NA"),
          'GT' = paste(unique(GT), collapse = ","),
          'ploidy' = paste(unique(ploidy), collapse = ","),
          'strand' = paste(unique(strand), collapse = ","),
          'defragmented' = paste(unique(defragmented), collapse = ","),
          'chr_hg19' = ifelse("chr_hg19" %in% names(wes_denovos_flipbook), paste(unique(chr_hg19), collapse = ","), "NA"),
          'start_hg19' = ifelse("start_hg19" %in% names(wes_denovos_flipbook), paste(unique(start_hg19), collapse = ","), "NA"),
          'end_hg19' = ifelse("end_hg19" %in% names(wes_denovos_flipbook), paste(unique(end_hg19), collapse = ","), "NA"),
          'sample_ori' = paste(unique(sample_ori), collapse = ","),
          'chrom_type_sample' = paste(unique(chrom_type_sample), collapse = ","),
          'start_win' = paste(unique(start_win), collapse = ","),
          'end_win' = paste(unique(end_win), collapse = ","),
          'name_merged' = paste(unique(name_merged), collapse = ","),
          'is_de_novo' =  ifelse("is_de_novo" %in% names(wes_denovos_flipbook), if (any(`Is de novo` == "no")) {
            "no"
          } else if (any(`Is de novo` == "yes")) {
            "yes"
          } else {
            "unsure"
          }, "NA"),
          #
          # 'is_de_novo' = if (any(`Is de novo` == "no")) {
          #   "no"
          # } else if (any(`Is de novo` == "yes")) {
          #   "yes"
          # } else {
          #   "unsure"
          # },
          'reason_unsure' = ifelse("reason_unsure" %in% names(wes_denovos_flipbook), paste(unique(`Reason Unsure - Follow up`), collapse = ";"), "NA"),
          'notes' = ifelse("notes" %in% names(wes_denovos_flipbook), paste(unique(Notes), collapse = ";"), "NA"),
          'name' = paste( unique(sample), unique(chrom), unique(start), unique(end), unique(svtype), sep = "_"),
          'path' = ifelse("notes" %in% names(wes_denovos_flipbook), paste(unique(Path), collapse = ";"), "NA"),
          ) -> wes_denovo_final

##Add columns for annotation
wes_denovo_final$SVTYPE <- wes_denovo_final$svtype
wes_denovo_final$CPX_TYPE_VCF <- ""
wes_denovo_final$CHR2 <- ""
wes_denovo_final$END2 <- ""

##Exclude NO calls
wes_denovo_final <- subset(wes_denovo_final, is_de_novo != "no")

##Write output
write.table(wes_denovo_final, paste0(output_dir, "/denovo_wes-", release, ".for_annotation.bed"), sep = "\t", quote = F, row.names = F, col.names = T)