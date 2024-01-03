#!/usr/bin/env Rscript

##This script reformats VCF file for PE evidence input
##Author: Alba Sanchis-Juan

#arguments
args = commandArgs(trailingOnly=TRUE)

input_file <- args[1]
prefix <- args[2]

#load libraries
library(data.table)
library(tidyr)

#read data
df <- fread(input_file)

df %>%
  subset(SVTYPE == "CTX") %>%
  separate_rows(samples, sep = ",") -> ctx

ctx_vcf_for_pe <- subset(ctx,
               select = c(`#chrom`, start, CHR2, END, samples, AC, name))

ctx_vcf_for_raw_ovl <- subset(ctx,
               select = c(`#chrom`, start, end, CHR2, END, samples, AC, name))

ctx_vcf_for_raw_ovl$`#chrom` <- paste0(ctx_vcf_for_raw_ovl$`#chrom`, "_", ctx_vcf_for_raw_ovl$samples)

#write output
write.table(ctx_vcf_for_pe, paste0(prefix, ".ctx.refForPE.split.bed"), sep = "\t", quote = F, row.names = F, col.names = F)
write.table(ctx_vcf_for_raw_ovl, paste0(prefix, ".ctx.refForOverlap.split.bed"), sep = "\t", quote = F, row.names = F, col.names = F)