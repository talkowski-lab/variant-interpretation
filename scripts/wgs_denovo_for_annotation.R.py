#!/usr/bin/env Rscript
#Author: Alba Sanchis Juan

#This scripts reformats the SV denovo release for annotation

#Load libraries
library(data.table)
library(optparse)

# Define Input Arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Final denovo SVs", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Final denovo SVs prep for annotation", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

#Define input parameteres
input_file <- opt$input
output_file <- opt$output

#Read file
input <- fread(input_file)

#Reformat
input_no_CPX_manual <- subset(input, is_highly_complex != TRUE, select = c(chrom, start, end, name, sample, svtype, SVTYPE, CPX_TYPE, CPX_INTERVALS, CHR2, END2))
names(input_no_CPX_manual)[c(8,9)] <- c("CPX_TYPE_VCF", "CPX_INTERVALS_VCF")
input_no_CPX_manual[is.na(input_no_CPX_manual$SVTYPE),]$SVTYPE <- input_no_CPX_manual[is.na(input_no_CPX_manual$SVTYPE),]$svtype

input_CPX_manual <- subset(input, is_highly_complex == TRUE, select = c(chrom, start, end, name, sample, svtype, SVTYPE, CPX_TYPE_VCF, CPX_INTERVALS_VCF, CHR2, END2))

#Merge
input <- rbind(input_no_CPX_manual, input_CPX_manual)

#Write
write.table(input, output_file, sep = "\t", quote = F, row.names = F)