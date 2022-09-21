#!/usr/bin/env Rscript
#Author: Alba Sanchis Juan

#This script reformats the output bed files from module 01 with raw call evidence
#to match desired input for de novo SV filtering

#Define args and load libraries
args = commandArgs(trailingOnly=TRUE)
require(plyr)
library(data.table)
library(tidyverse)

input_bed <- args[1]
out_bed <- args[2]

#Read input bed
bed <- fread(input_bed)

#Split into one sample per row
bed %>%
  subset(select = c(`#chrom`, start, end, SVTYPE, samples)) %>%
  separate_rows(samples, sep = ",", convert = T) -> bed_split

#Reformat chromosome column
bed_split$CHROM <- paste0(bed_split$`#chrom`, "_", bed_split$SVTYPE, "_", bed_split$samples)

#Select specific samples
bed_subset <- subset(bed_split, select = c(CHROM, start, end, SVTYPE, samples))

#Write output file
write.table(bed_split, out_file, sep = "\t", quote = F, row.names = F, col.names = F)