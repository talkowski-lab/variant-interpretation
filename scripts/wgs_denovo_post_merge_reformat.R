#!/usr/bin/env Rscript
#Author: Alba Sanchis Juan

#This scripts reformats the SV denovo release post merging
#Appends SV TLOCS and MOSAICS

#Load libraries
library(tidyverse)
library(data.table)
library(optparse)

# Define Input Arguments
option_list = list(
  make_option(c("-d", "--denovo"), type="character", default=NULL,
              help="De novo file", metavar="character"),
  make_option(c("-e", "--merged"), type="character", default=NULL,
              help="De novo pre-merge reformatted file", metavar="character"),
  make_option(c("-t", "--translocations"), type="character", default=NULL,
              help="Tlocs file", metavar="character"),
  make_option(c("-m", "--mosaics"), type="character", default=NULL,
              help="Mosaics file", metavar="character"),
  make_option(c("-w", "--mosaics2"), type="character", default=NULL,
              help="Mosaics file 2", metavar="character"),
  make_option(c("-g", "--genomicdisorders"), type="character", default=NULL,
              help="GDs file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="Output file name [default= %default]", metavar="character"),
  make_option(c("-a", "--aneuploidies"), type="character", default=NULL,
              help="Aneuploidies file name [default= %default]", metavar="character"),
  make_option(c("-r", "--remove"), type="character", default=NULL,
              help="Remove SVs file name [default= %default]", metavar="character"),
  make_option(c("-c", "--add"), type="character", default=NULL,
              help="Add SVs file name [default= %default]", metavar="character"),
  make_option(c("-q", "--qc"), type="character", default=NULL,
              help="QC file name [default= %default]", metavar="character"),
  make_option(c("-p", "--ped"), type="character", default=NULL,
              help="Pedigree file name [default= %default]", metavar="character"),
  make_option(c("-b", "--baf"), type="character", default=NULL,
              help="BAF file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

#Define input parameteres
denovo_file <- opt$denovo
merged_file <- opt$merged
tlocs_file <- opt$translocations
mosaics_file <- opt$mosaics
mosaics_file2 <- opt$mosaics2
gds_file <- opt$genomicdisorders
out_file <- opt$out
aneuploidies_file <- opt$aneuploidies
remove_file <- opt$remove
add_file <- opt$add
qc_file <- opt$qc
ped_file <- opt$ped
baf_file <- opt$baf

#Read files
denovo <- fread(denovo_file)
wgs_svs_merged <- fread(merged_file, col.names = c("chrom_sample", "start_fix", "end_fix", "sample_name_collapsed"))
tlocs <- fread(tlocs_file)
mosaics <- fread(mosaics_file)
mosaics2 <- fread(mosaics_file2)
gds <- fread(gds_file)
aneuploidies <- fread(aneuploidies_file)
remove_svs <- fread(remove_file)
add_svs <- fread(add_file, fill = TRUE)
qc <- fread(qc_file)
ped <- fread(ped_file)
baf <- fread(baf_file, fill = TRUE)

##Keep only real SVs
wgs_svs_real <- subset(denovo, result_final == "yes")

##Load after merging and reformat
wgs_svs_merged$sample_name <- wgs_svs_merged$sample_name_collapsed

wgs_svs_merged_ref <- separate_rows(wgs_svs_merged, sample_name, sep = ",", convert = TRUE)

wgs_svs_real_coll <- merge(wgs_svs_real, unique(subset(wgs_svs_merged_ref, select = c(start_fix, end_fix, sample_name_collapsed, sample_name))), by = "sample_name", all.x = T, all.y = F)

wgs_svs_real_coll$name_fix <- wgs_svs_real_coll$sample_name_collapsed

wgs_svs_real_coll$chrom_sample <- paste0(wgs_svs_real_coll$chrom, "_", wgs_svs_real_coll$sample)

wgs_svs_real_coll$name_fix <- apply(wgs_svs_real_coll, 1, function(row){
    paste(as.vector(unique(subset(wgs_svs_real_coll, sample_name_collapsed == row['sample_name_collapsed'],
                                  select = c(chrom_sample, start_fix, end_fix)))), collapse = "_")
})

wgs_svs_real_coll$is_fragmented <- wgs_svs_real_coll$name_fix %in% plyr::count(wgs_svs_real_coll$name_fix)[plyr::count(wgs_svs_real_coll$name_fix)$freq >1,]$x

subset(wgs_svs_real_coll) %>%
  group_by(name_fix) %>%
  summarise(chrom = unique(chrom), start = unique(start_fix)
            , end = unique(end_fix)
            , name_fix = unique(name_fix)
            , name = paste(unique(name), collapse = ",")
            , Path = paste(unique(Path), collapse = ",")
            , sample = unique(sample)
            , batch = unique(batch), AF = max(AF), AC = max(AC), AN = unique(AN), num_reviewers = max(num_reviewers), result_final = unique(result_final)
            , reviewers = paste(unique(unlist(strsplit(reviewers, "&"))), collapse = "&")
            , svtype = paste(unique(svtype), collapse = ",")
            , SVTYPE = paste(unique(SVTYPE), collapse = ",")
            , CHR2 = paste(unique(CHR2), collapse = ",")
            , END2 = paste(unique(END2), collapse = ",")
            , EVIDENCE_FIX = paste(unique(unlist(strsplit(EVIDENCE_FIX, ","))), collapse = ",")
            , chrom_sample = paste(unique(chrom_sample), collapse = ",")
            , PREDICTED_COPY_GAIN = paste(unique(unlist(strsplit(PREDICTED_COPY_GAIN, ","))), collapse = ",")
            , PREDICTED_DUP_PARTIAL = paste( unique(unlist(strsplit(PREDICTED_DUP_PARTIAL, ","))) , collapse = ",")
            , PREDICTED_INTERGENIC = all(PREDICTED_INTERGENIC)
            , PREDICTED_INTRAGENIC_EXON_DUP = paste( unique(unlist(strsplit(PREDICTED_INTRAGENIC_EXON_DUP, ","))) , collapse = ",")
            , PREDICTED_INTRONIC = paste( unique(unlist(strsplit(PREDICTED_INTRONIC, ","))) , collapse = ",")
            , PREDICTED_INV_SPAN = paste(unique(unlist(strsplit(PREDICTED_INV_SPAN, ","))), collapse = ",")
            , PREDICTED_LOF = paste(unique(unlist(strsplit(PREDICTED_LOF, ","))), collapse = ",")
            , PREDICTED_NEAREST_TSS = paste(unique(unlist(strsplit(PREDICTED_NEAREST_TSS, ","))), collapse = ",")
            , PREDICTED_PARTIAL_EXON_DUP = paste(unique(unlist(strsplit(PREDICTED_PARTIAL_EXON_DUP, ","))), collapse = ",")
            , PREDICTED_PROMOTER = paste(unique(unlist(strsplit(PREDICTED_PROMOTER, ","))), collapse = ",")
            , PREDICTED_TSS_DUP = paste(unique(unlist(strsplit(PREDICTED_TSS_DUP, ","))), collapse = ",")
            , PREDICTED_UTR = paste(unique(unlist(strsplit(PREDICTED_UTR, ","))), collapse = ",")
            , CPX_TYPE = paste(unique(unlist(strsplit(CPX_TYPE, ","))), collapse = ",")
            , CPX_INTERVALS = paste(unique(unlist(strsplit(CPX_INTERVALS, ","))), collapse = ",")
            ) -> wgs_svs_real_post_merge

wgs_svs_real_post_merge[grepl(",", wgs_svs_real_post_merge$svtype),]$SVTYPE <- "MULTIPLE"
wgs_svs_real_post_merge$SVLEN <- wgs_svs_real_post_merge$end - wgs_svs_real_post_merge$start


###############
##ADD MOSAICS##
###############
#Reformat MOSAICS input file
names(mosaics)[names(mosaics) %in% tail(names(mosaics), 5)] <- c("chr_ovl", "start_ovl", "end_ovl", "name_ovl", "size_ovl")
mosaics$chrom <- sapply(mosaics$chrom, function(x) strsplit(x, "_")[[1]][1])
names(mosaics)[names(mosaics) == "samples"] <- "sample"
mosaics$is_mosaic <- TRUE

#Get Mosaic variant names
mosaics_callset_names <- as.vector(unlist(sapply(subset(mosaics, size_ovl != "0")$name_ovl, function(x) unlist(strsplit(x, ",")[[1]]) )))
mosaics_manual_paths <- subset(mosaics2, `Is de novo` == "yes")$Path

#Flag variants in callset that overlap a GD, based on variant name
wgs_svs_real_post_merge$is_mosaic <- FALSE
wgs_svs_real_post_merge$is_mosaic <- apply(wgs_svs_real_post_merge, 1, function(row){
  any( paste0(row['sample'], "_", strsplit(row['name'], ",")[[1]]) %in% mosaics_callset_names |
         paste0(row['sample'], "_", strsplit(row['Path'], ",")[[1]]) %in% mosaics_manual_paths)
})
wgs_svs_real_post_merge$is_mosaic_manual <- FALSE

#Add new MOSAIC variants not previously in callset
mosaics_new <- subset(mosaics, size_ovl == "0", select = names(mosaics)[names(mosaics) %in% names(wgs_svs_real_post_merge)])
mosaics_new$is_mosaic_manual <- TRUE
mosaics_new$result_final <- "yes"
mosaics_new$reviewers <- "Arthur&Harrison"
mosaics_new$num_reviewers <- "2"
mosaics_new <- subset(mosaics_new, select = names(mosaics_new)[names(mosaics_new) %in% names(wgs_svs_real_post_merge)])

wgs_svs_real_post_merge_mosaics <- plyr::rbind.fill(wgs_svs_real_post_merge, mosaics_new)


#############
##ADD TLOCS##
#############
#Reformat TLOCS input file
names(tlocs)[names(tlocs) %in% tail(names(tlocs), 5)] <- c("chr_ovl", "start_ovl", "end_ovl", "name_ovl", "size_ovl")
tlocs$chrom <- sapply(tlocs$chrom, function(x) strsplit(x, "_")[[1]][1])
names(tlocs)[names(tlocs) == "samples"] <- "sample"
tlocs$is_tloc <- TRUE

#Get Tlocs variant names
tlocs_callset_names <- as.vector(unlist(sapply(subset(tlocs, size_ovl != "0")$name_ovl, function(x) unlist(strsplit(x, ",")[[1]]) )))

#Flag Tlocs in callset based on variant name
wgs_svs_real_post_merge_mosaics$is_tloc <- FALSE
wgs_svs_real_post_merge_mosaics$is_tloc <- apply(wgs_svs_real_post_merge_mosaics, 1, function(row){
  any( paste0(row['sample'], "_", strsplit(row['name'], ",")[[1]]) %in% tlocs_callset_names)
})
wgs_svs_real_post_merge_mosaics$is_tloc_manual <- FALSE

#Add new Tlocs not previously in callset
tlocs_new <- subset(tlocs, size_ovl == "0", select = names(tlocs)[names(tlocs) %in% names(wgs_svs_real_post_merge_mosaics)])
tlocs_new$is_tloc_manual <- TRUE
tlocs_new$result_final <- "yes"
tlocs_new$reviewers <- "Nehir"
tlocs_new$num_reviewers <- "1"
tlocs_new <- subset(tlocs_new, select = names(tlocs_new)[names(tlocs_new) %in% names(wgs_svs_real_post_merge_mosaics)])

wgs_svs_real_post_merge_mosaics_tlocs <- plyr::rbind.fill(wgs_svs_real_post_merge_mosaics, tlocs_new)

###########
##ADD GDS##
###########
#Reformat GD input file
names(gds)[names(gds) %in% tail(names(gds), 5)] <- c("chr_ovl", "start_ovl", "end_ovl", "name_ovl", "size_ovl")
gds$chrom <- sapply(gds$chrom, function(x) strsplit(x, "_")[[1]][1])
#names(gds)[names(gds) == "samples"] <- "sample"
gds$is_GD <- TRUE

#Get GD variant names
gds_callset_names <- as.vector(unlist(sapply(subset(gds, size_ovl != "0")$name_ovl, function(x) unlist(strsplit(x, ",")[[1]]) )))

#Flag variants in callset that overlap a GD, based on variant name
wgs_svs_real_post_merge_mosaics_tlocs$is_GD <- FALSE
wgs_svs_real_post_merge_mosaics_tlocs$is_GD <- apply(wgs_svs_real_post_merge_mosaics_tlocs, 1, function(row){
  any( paste0(row['sample'], "_", strsplit(row['name'], ",")[[1]]) %in% gds_callset_names)
})

#Add new GD variants not previously in callset
wgs_svs_real_post_merge_mosaics_tlocs$is_GD_manual <- FALSE

gds_new <- subset(gds, size_ovl == "0", select = names(gds)[names(gds) %in% names(wgs_svs_real_post_merge_mosaics_tlocs)])
gds_new$is_GD_manual <- TRUE
gds_new$result_final <- "yes"
gds_new$reviewers <- "Nehir"
gds_new$num_reviewers <- "1"
gds_new <- subset(gds_new, select = names(gds_new)[names(gds_new) %in% names(wgs_svs_real_post_merge_mosaics_tlocs)])

wgs_svs_real_post_merge_mosaics_tlocs_gds <- plyr::rbind.fill(wgs_svs_real_post_merge_mosaics_tlocs, gds_new)

# ##Add GD name
# # gd_names <-
#   apply(wgs_svs_real_post_merge_mosaics_tlocs_gds, 1, function(row){
# # apply(subset(wgs_svs_real_post_merge_mosaics_tlocs_gds, is_GD == TRUE & GD_name == "")[1,], 1, function(row){
#   if(row['is_GD'] == TRUE){
#     # row_name <- strsplit(as.character(row['name']), ",")[[1]]
#     row_sample <- row['sample']
#     row_chrom <- row['chrom']
#     paste0(subset(gds, sample == row_sample & chrom == row_chrom)$GD, collapse = ",")
#   }else{
#     ""
#   }
# })
#
# gd_names[lengths(gd_names) == 0] <- NA
# wgs_svs_real_post_merge_mosaics_tlocs_gds$GD_name <- unlist(gd_names, use.names = FALSE)

####################
##ADD ANEUPLOIDIES##
####################
##Remove aneuploidies in callset
wgs_svs_real_post_merge_mosaics_tlocs_gds$chrom_sample_type <- paste0(wgs_svs_real_post_merge_mosaics_tlocs_gds$chrom_sample, "_", wgs_svs_real_post_merge_mosaics_tlocs_gds$SVTYPE)
aneuploidies$chrom_sample_type <- paste0(aneuploidies$chrom, "_", aneuploidies$sample, "_", aneuploidies$svtype)

wgs_svs_real_post_merge_mosaics_tlocs_gds <- subset(wgs_svs_real_post_merge_mosaics_tlocs_gds, !chrom_sample_type %in% aneuploidies$chrom_sample_type)
wgs_svs_real_post_merge_mosaics_tlocs_gds$is_aneuploidy <- FALSE

##Add clean aneuploidies
aneuploidies$is_aneuploidy <- TRUE
aneuploidies$result_final <- "yes"
aneuploidies$reviewers <- "Alba"
aneuploidies$num_reviewers <- "1"

wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies <- plyr::rbind.fill(wgs_svs_real_post_merge_mosaics_tlocs_gds, aneuploidies)

######################
##Final clean up SVs##
######################
#Manual remove SVs
wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies_clean1 <- subset(wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies, !name %in% remove_svs$name)

#Manual add SVs
wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies_clean2 <- plyr::rbind.fill(wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies_clean1, add_svs)

##Fix is_highly_complex and is_unbalanced_tloc
wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies_clean2[is.na(wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies_clean2$is_highly_complex),]$is_highly_complex <- FALSE
wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies_clean2[is.na(wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies_clean2$is_unbalanced_tloc),]$is_unbalanced_tloc <- FALSE

##Remove samples not in batches considered in this callset
wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies_clean3 <- subset(wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies_clean2, sample %in% ped$sv_pipeline_id)

##Flag if sample has an issue
wgs_samples_parental_error <- subset(qc, maternal_error | paternal_error)$sample_id
wgs_samples_parental_error_ids <- subset(ped, snv_pipeline_id %in% wgs_samples_parental_error)$sv_pipeline_id

wgs_samples_sex_error <- subset(qc, sex_error == TRUE)$sample_id
wgs_samples_sex_error_ids <- subset(ped, snv_pipeline_id %in% wgs_samples_sex_error)$sv_pipeline_id

wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies_clean3$relatedness_error <- wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies_clean3$sample %in% wgs_samples_parental_error_ids
wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies_clean3$sex_error <- wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies_clean3$sample %in% wgs_samples_sex_error_ids

wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies_clean3$any_error <- wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies_clean3$relatedness_error | wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies_clean3$sex_error


##DEFINE FINAL FILE##
wgs_svs_final <- wgs_svs_real_post_merge_mosaics_tlocs_gds_aneuploidies_clean3

##ADD ALL GENES COLUMN
#wgs_svs_final[wgs_svs_final %in% c(",", "NA", "NA ,", NA, ", ", " ,", " ")] <- ""
genes_colnames <- grep("PREDICTED", names(wgs_svs_final), value = T)[!grep("PREDICTED", names(wgs_svs_final), value = T) %in% c("PREDICTED_INTERGENIC", "PREDICTED_NONCODING_BREAKPOINT", "PREDICTED_NONCODING_SPAN")]

wgs_svs_final$all_genes <- apply(wgs_svs_final, 1, function(row){
  paste(unique(sort(row[genes_colnames][!row[genes_colnames] %in% c("", ",", "NA", "NA ,", NA, ", ", " ,", " ") ])), collapse = ",")
})


##Add name fix for empty fields
wgs_svs_final[is.na(wgs_svs_final$name_fix),]$name_fix <- paste(
  wgs_svs_final[is.na(wgs_svs_final$name_fix),]$chrom,
  wgs_svs_final[is.na(wgs_svs_final$name_fix),]$sample,
  wgs_svs_final[is.na(wgs_svs_final$name_fix),]$start,
  wgs_svs_final[is.na(wgs_svs_final$name_fix),]$end,
  sep = "_")


##ADD BAF INFORMATION
wgs_svs_final$sample_coords <- paste0(wgs_svs_final$sample, "_", wgs_svs_final$chrom, "_", wgs_svs_final$start, "_", wgs_svs_final$end)

#Merge number of informative snps to master
wgs_svs_final <- merge(wgs_svs_final, baf, by = "sample_coords", all.x = T, all.y = F)

#Write file
write.table(wgs_svs_final, out_file, sep = "\t", quote = F, row.names = F)