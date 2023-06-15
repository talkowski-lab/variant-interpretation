#!/usr/bin/env Rscript
#Author: Alba Sanchis Juan and Nicole Calamari

#Define args and load libraries
args = commandArgs(trailingOnly=TRUE)
require(plyr)
library(data.table)
library(tidyverse)

input_bed <- args[1]

#Read input bed
bed <- fread(input_bed)

#adding back column names
colnames(bed) <- c("chrom", "start", "end", "name", "svtype", "samples")

#Split into one sample per row
bed %>%
  subset(select = c(chrom, start, end, name, svtype, samples)) %>%
  separate_rows(samples, sep = ",", convert = T) -> bed_split

#Write output file
write.table(bed_split, "single_sample.bed", sep = "\t", quote = F, row.names = F, col.names = T)
