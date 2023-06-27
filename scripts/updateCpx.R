#!/usr/bin/env Rscript

##Transpose CPX script##

##Load libraries
library(tidyr)
library(data.table)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

input_bed <- args[1]

##Read data
denovo <- fread(input_bed)
out_path <- "final.denovo.merged.cpx_split.bed"

##Split dataframe into complex and non-complex calls
denovo_non_cpx <- subset(denovo, svtype != "CPX", select = c("chrom", "start", "end", "name", "svtype", "sample"))
denovo_cpx <- subset(denovo, svtype == "CPX")

##Reformat cpx calls
denovo_cpx %>% 
  separate_rows(CPX_INTERVALS, sep = ",") %>%
  select(c("chrom", "start", "end", "name", "svtype", "sample", "CPX_INTERVALS")) %>%
  separate(CPX_INTERVALS, sep = "_", into = c("cpx_type", "cpx_coord"), ) %>%
  separate(cpx_coord, sep = ":", into = c("cpx_chr", "cpx_pos"), ) %>%
  separate(cpx_pos, sep = "-", into = c("cpx_start", "cpx_end"), ) -> denovo_cpx

##Update name of cpx calls and chane column names
denovo_cpx$name <- paste0(denovo_cpx$name, "_", denovo_cpx$cpx_type)
denovo_cpx <- subset(denovo_cpx, select = c("cpx_chr", "cpx_start", "cpx_end", "name", "cpx_type", "sample"))
names(denovo_cpx) <- names(denovo_non_cpx)

##Merge all de novo calls
denovo_split <- rbind(denovo_cpx, denovo_non_cpx)

##Write output
write.table(denovo_split, out_path, sep = "\t", quote = F, row.names = F, col.names = T)