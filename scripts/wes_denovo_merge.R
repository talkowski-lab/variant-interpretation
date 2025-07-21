#!/usr/bin/env Rscript
#Author: Alba Sanchis-Juan

##This script merges the CNV WES callsets

#load libraries
library(data.table)
library(optparse)

#parameters
option_list = list(
  make_option(c("-c", "--callsets"), type="character", default=NULL,
              help="comma sep list of callsets", metavar="character"),
  make_option(c("-r", "--release"), type="character", default=NULL,
              help="release date", metavar="character"),
  make_option(c("-o", "--outputdir"), type="character", default=NULL,
              help="output file", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

#Define input parameteres and read files
callsets <- opt$callsets
release <- opt$release
output_dir <- opt$outputdir

wes_denovos <- do.call(plyr::rbind.fill, lapply(strsplit(callsets, ",")[[1]], fread))

wes_denovos$sample_ori <- wes_denovos$sample
wes_denovos[!is.na(wes_denovos$sample_fix),]$sample <- wes_denovos[!is.na(wes_denovos$sample_fix),]$sample_fix

##Write output for merging
wes_denovos$chrom_type_sample <- paste0(wes_denovos$chr, "_", wes_denovos$svtype, "_", wes_denovos$sample)

wes_denovos_win <- cbind(wes_denovos, do.call(rbind, apply(wes_denovos, 1, function(row){
  row_win <- 0.2*(as.numeric(row['end']) - as.numeric(row['start']))
  data.table("start_win" = as.integer(as.numeric(row['start'])-row_win), "end_win" = as.integer(as.numeric(row['end'])+row_win))
})))

wes_denovos_win$start_win <- ifelse(wes_denovos_win$start_win < 0, 1, wes_denovos_win$start)

##write output
write.table(wes_denovos_win, paste0(output_dir, "/denovo_wes-", release, ".bed"), sep = "\t", quote = F, row.names = F)

##Write subset for merging
wes_denovos_win <- subset(wes_denovos_win, select = c(chrom_type_sample, start_win, end_win, name))
wes_denovos_win <- wes_denovos_win[order(wes_denovos_win$chrom_type_sample, wes_denovos_win$start_win),]
write.table(wes_denovos_win, paste0(output_dir, "/denovo_wes_for_merging-", release, ".bed"), sep = "\t", quote = F, row.names = F, col.names = F)