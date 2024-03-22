#!/usr/bin/env Rscript
#Author: Alba Sanchis-Juan
##This script creates the tlocs summary table

#Load libraries
library(data.table)
library(optparse)

##Arguments
option_list = list(
  make_option(c("-s", "--samples"), type="character", default=NULL,
              help="samples file name", metavar="character"),
  make_option(c("-i", "--info"), type="character", default=NULL,
              help="info name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=".",
              help="output file name", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#Read files
samples <- fread("~/tmp/debug3/mg_batch03.chr7.final_cleanup_CTX_chr7_1_samples.bed")
info <- fread("~/tmp/debug3/mg_batch03.chr7.final_cleanup_CTX_chr7_1.info")
out <- opt$out

#Create summary table
summary_df <- do.call(rbind, apply(samples, 1, function(row){
  sample_ID <- row['V6']
  tloc_ID <- row['V5']

  BP1_1_plus <- min(subset(info, V7 == sample_ID & V3 == "+")$V2)
  BP1_2_plus <- max(subset(info, V7 == sample_ID & V3 == "+")$V2)
  BP1_1_minus <- min(subset(info, V7 == sample_ID & V3 == "-")$V2)
  BP1_2_minus <- max(subset(info, V7 == sample_ID & V3 == "-")$V2)
  BP1_int_plus <- BP1_1_plus-BP1_2_plus
  BP1_int_minus <- BP1_1_minus-BP1_2_minus
  BP2_1_plus <- min(subset(info, V7 == sample_ID & V6 == "+")$V5)
  BP2_2_plus <- max(subset(info, V7 == sample_ID & V6 == "+")$V5)
  BP2_1_minus <- min(subset(info, V7 == sample_ID & V6 == "-")$V5)
  BP2_2_minus <- max(subset(info, V7 == sample_ID & V6 == "-")$V5)
  BP2_int_plus <- BP2_1_plus-BP2_2_plus
  BP2_int_minus <- BP2_1_minus-BP2_2_minus

  return(data.frame(BP1_1_plus, BP1_2_plus, BP1_1_minus, BP1_2_minus, BP1_int_plus, BP1_int_minus,
                    BP2_1_plus, BP2_2_plus, BP2_1_minus, BP2_2_minus, BP2_int_plus, BP2_int_minus))

}))

#Append to main samples table
samples_summary <- cbind(samples, summary_df)

#Write output
write.table(samples_summary, out, sep = "\t", quote = F, row.names = F, col.names = F)