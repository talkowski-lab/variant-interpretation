#!/usr/bin/env Rscript

##Overwrite callset with new annotations

#Load libraries
# library(tidyverse)
library(data.table)
library(optparse)

# Define Input Arguments
option_list = list(
  make_option(c("-s", "--svcallset"), type="character", default=NULL,
              help="De novo SV WGS file", metavar="character"),
  make_option(c("-b", "--bed"), type="character", default=NULL,
              help="Annotated bed file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Final annotated bed file", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

#Define input parameteres and read files
svcallset_file <- opt$svcallset
bed_file <- opt$bed
output_file <- opt$output

svcallset <- fread(svcallset_file)
bed <- fread(bed_file)

#Make annotations unique by name
bed2 <- unique(subset(bed, select = c("name", grep("PREDICTED", names(bed), value = T))))

##Separate svcallset with NA names before merging
svcallset_with_name <- subset(svcallset, !is.na(name))
svcallset_without_name <- subset(svcallset, is.na(name))

svcallset_with_name <- subset(svcallset_with_name, select = names(svcallset_with_name)[!names(svcallset_with_name) %in% names(bed2)[names(bed2) != "name"]])
svcallset_with_name_annot <- merge(svcallset_with_name, bed2, by = "name", all.x = T, all.y = F)

##Add annotations to sv calls without name (aneuploidies) - to fix this at some point by giving it a variant name
svcallset_without_name$chrom_coords_type <- paste(svcallset_without_name$chrom, svcallset_without_name$start, svcallset_without_name$end, svcallset_without_name$svtype,sep = "_")
bed$chrom_coords_type <- paste(bed$`#chrom`, bed$start+1, bed$end, bed$svtype,sep = "_")
bed3 <- unique(subset(bed, select = c("chrom_coords_type", grep("PREDICTED", names(bed), value = T))))

svcallset_without_name <- subset(svcallset_without_name, select = names(svcallset_without_name)[!names(svcallset_without_name) %in% names(bed2)[names(bed2) != "name"]])
svcallset_without_name_annot <- merge(svcallset_without_name, bed3, by = "chrom_coords_type", all.x = T, all.y = F)
svcallset_without_name_annot$chrom_coords_type <- NULL

##ANNOTATIONS ARE MISSING FOR ANEUPLOIDIES IN CHR X AND CHR Y
svcallset_annot <- rbind(svcallset_with_name_annot, svcallset_without_name_annot)

##Flag if is any manual
svcallset_annot$is_any_manual <- apply(svcallset_annot, 1, function(row){
  colnames <- c("is_mosaic_manual", "is_tloc_manual", "is_GD_manual", "is_aneuploidy")
  any(row[colnames] == TRUE)
})

##Reorder Column names
svcallset_annot2 <- subset(svcallset_annot, select = c(names(svcallset), names(svcallset_annot)[!names(svcallset_annot) %in% names(svcallset)]) )

#Write file
write.table(svcallset_annot2, output_file, sep = "\t", quote = F, row.names = F)