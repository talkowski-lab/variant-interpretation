#!/usr/bin/env Rscript

##This script reformats raw files for PE evidence input
##Author: Alba Sanchis-Juan

#arguments
args = commandArgs(trailingOnly=TRUE)

input_file <- args[1]
output_file <- args[2]

#load libraries
library(data.table)
library(tidyverse)

#read data
raw_input <- fread(input_file, header = F)

#reformat columns
raw_input$V1 <- sapply(raw_input$V1, function(x) strsplit(x, "_")[[1]][1] )
raw_input$V2 <- as.numeric(raw_input$V2) + 20

raw_input$coords <- paste(raw_input$V1, raw_input$V2, raw_input$V8, raw_input$V9, sep = "_")
raw_input$num_samples <- sapply(raw_input$coords, function(x) nrow(subset(raw_input, coords == x)) )

#subset columns for output
raw_input %>%
  select(V1, V2, V8, V9, V6, num_samples, V4) -> raw_output

#write output
write.table(raw_output, output_file, sep = "\t", quote = F, row.names = F, col.names = F)