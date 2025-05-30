#!/usr/bin/env Rscript
#This script reformats and prepares gCNV bed file for merging with bedtools

#load libraries
library(data.table)
library(tidyverse)
library(optparse)
options(scipen = 999)

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="comma sep list of callsets", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default=NULL,
              help="output file", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

#Define input parameters
input <- opt$input
output <- opt$output

#Read file
callset <- fread(input)

##Add extra info
callset$chrom_type_sample <- paste0(callset$chr, "_", callset$svtype, "_", callset$sample)

#Compute window size for merging
row_win <- 0.2 * (as.numeric(callset$end) - as.numeric(callset$start))

#Add the new columns directly
callset[, `:=` (
  start_win = as.integer(as.numeric(start) - row_win),
  end_win   = as.integer(as.numeric(end) + row_win)
)]

callset$start_win <- ifelse(callset$start_win < 0, 1, callset$start_win)

##Write subset for merging
write.table(callset, paste0(prefix, ".ref.bed"),
            sep = "\t", quote = F, row.names = F, col.names = F)

callset_ref <- subset(callset, select = c(chrom_type_sample, start_win, end_win, name))
callset_ref <- callset_ref[order(callset_ref$chrom_type_sample, callset_ref$start_win),]
write.table(callset_ref, paste0(prefix, "ref.for_merging.bed"),
            sep = "\t", quote = F, row.names = F, col.names = F)