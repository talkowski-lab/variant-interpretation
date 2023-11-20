#!/usr/bin/env Rscript
#Author: Alba Sanchis-Juan
#This script prep files for running de novo pipeline

#Load libraries
library(optparse)
library(data.table)

##Arguments
option_list = list(
  make_option(c("-e", "--entity"), type="character", default=NULL,
              help="entity file name", metavar="character"),
  make_option(c("-m", "--membership"), type="character", default=NULL,
              help="membership file name", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL,
              help="sample file name", metavar="character"),
  make_option(c("-p", "--pedigree"), type="character", default=NULL,
              help="pedigree file name", metavar="character"),
  make_option(c("-o", "--out_dir"), type="character", default=".",
              help="output directory", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

entity <- fread(opt$entity)
membership <- fread(opt$membership)
sample <- fread(opt$sample)
ped <- fread(opt$pedigree)
out_dir <- opt$out_dir

#Make sample_batches file
sample_batches <- subset(membership, select = c(sample, `membership:sample_set_id`))

#Make batch_bincov_index file
batch_bincov_index <- subset(entity, select = c(`entity:sample_set_id`, merged_bincov, merged_bincov_index))
batch_bincov <- subset(entity, select = c(`entity:sample_set_id`, merged_bincov))

#Make batch_medianfile file
batch_medianfile <- subset(entity, select = c(`entity:sample_set_id`, median_cov))

#Make batch_depth_raw_file file
batch_depth_raw_file <- subset(entity, select = c(`entity:sample_set_id`, clustered_depth_vcf))

#Make batch_raw_file file
batch_raw_file <- subset(entity, select = c(`entity:sample_set_id`, clustered_manta_vcf, clustered_melt_vcf, clustered_wham_vcf))
batch_raw_file <- melt(batch_raw_file, id.vars = "entity:sample_set_id")
batch_raw_file$variable <- NULL

#Make sample_crai_cram file
sample_crai_cram <- subset(sample, select = c(`entity:sample_id`, crai, cram))

#Reformat ped file
ped <- ped[,c(1:6)]
names(ped) <- c("FamID", "IndividualID", "FatherID", "MotherID", "Gender", "Affected")

#Save files
write.table(sample_batches, file.path(out_dir, "sample_batches.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
write.table(batch_bincov_index, file.path(out_dir, "batch_bincov_index.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
write.table(batch_bincov, file.path(out_dir, "batch_bincov.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
write.table(batch_medianfile, file.path(out_dir, "batch_medianfile.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
write.table(batch_depth_raw_file, file.path(out_dir, "batch_depth_raw_files.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
write.table(batch_raw_file, file.path(out_dir, "batch_raw_files.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
write.table(sample_crai_cram, file.path(out_dir, "sample_crai_cram.txt"), sep = "\t", quote = F, row.names = F, col.names = F)
write.table(ped, file.path(out_dir, "pedigree.txt"), sep = "\t", quote = F, row.names = F, col.names = T)