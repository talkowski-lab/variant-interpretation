#!/usr/bin/env Rscript
#Author: Alba Sanchis Juan

#This scripts aggregates SV de novo calls across different cohorts and
#Appends SV adjudication results

#Load libraries
library(tidyverse)
library(data.table)
library(optparse)
# library(bedr)

# Define Input Arguments
option_list = list(
  make_option(c("-d", "--denovo"), type="character", default=NULL,
              help="De novo file of files", metavar="character"),
  # make_option(c("-l", "--outliers"), type="character", default=NULL,
  #             help="De novo outliers file of files", metavar="character"),
  make_option(c("-n", "--denovo_cohort_names"), type="character", default=NULL,
              help="De novo outliers file of files", metavar="character"),
  make_option(c("-f", "--flipbook"), type="character", default=NULL,
              help="Flipbook results file of files", metavar="character"),
  make_option(c("-l", "--flipbook_metadata"), type="character", default=NULL,
              help="Flipbook results file of files", metavar="character"),
  make_option(c("-o", "--out"), type="character", default="out.txt",
              help="Output file name [default= %default]", metavar="character"),
  make_option(c("-q", "--qc"), type="character", default=NULL,
              help="QC file from the pipeline", metavar="character"),
  make_option(c("-p", "--ped"), type="character", default="out.txt",
              help="Pedigree file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

denovo_fof <- opt$denovo
denovo_cohorts_file <- opt$denovo_cohort_names
# outliers_fof <- opt$outliers
flipbook_fof <- opt$flipbook
out_file <- opt$out
qc_file <- opt$qc
ped_file <- opt$ped
flipbook_metadata_file <- opt$flipbook_metadata

#Functions
detachAllPackages <- function() {
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  package.list <- setdiff(package.list,basic.packages)
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
}

#Read Input Files and Reformat
qc_file <- fread(qc_file)
ped <- fread(ped_file)
denovo_cohorts <- fread(denovo_cohorts_file, header = F)
flipbook_metadata <- fread(flipbook_metadata_file, header = F)

denovo_files <- strsplit(denovo_fof, ",")[[1]]

denovo <- do.call(plyr::rbind.fill, lapply(denovo_files, function(f){
  df <- fread(f)
  df$batch <- subset(denovo_cohorts, V2 == basename(f))$V1
  df
}))

denovo$Path <- paste0(denovo$sample, "_", denovo$name, "_denovo.png")
denovo$sample_name <- paste0(denovo$sample, "_", denovo$name)

# #Read denovo OUTLIERS files and reformat
# denovo_outliers_files <- read.table(outliers_fof)
#
# denovo_outliers <- do.call(plyr::rbind.fill, lapply(denovo_outliers_files$V2, function(f){
#   df <- fread(f)
#   df$batch <- subset(denovo_outliers_files, V2 == f)$V1
#   # df$batch <- paste(strsplit(basename(f), "_")[[1]][1:2], collapse = "_")
#   df
# }))
#
# # denovo_outliers$Path <- paste0(denovo_outliers$sample, "_", denovo_outliers$name, "_denovo.png")
# denovo_outliers$sample_name <- paste0(denovo_outliers$sample, "_", denovo_outliers$name)

#Read Flipbook files and reformat
flipbook_files <- strsplit(flipbook_fof, ",")[[1]]

responses <- do.call(plyr::rbind.fill, lapply(flipbook_files, function(f){
  tb <- fread(f)
  tb$analyst <- subset(flipbook_metadata, V2 == basename(f))$V1
  tb$tiebreaker_reviewed <- subset(flipbook_metadata, V2 == basename(f))$V3
  tb
}))

responses$Verdict <- NULL
responses$Confidence <- NULL

names(responses) <- c("Path", "Is_de_novo", "Reason_unsure_follow_up", "Notes", "analyst", "tiebreaker_reviewed")

if(nrow(responses[responses$Is_de_novo == "",]) >0){
  responses[responses$Is_de_novo == "",]$Is_de_novo <- "no"
}

responses$analyst <- factor(responses$analyst, levels = unique(responses$analyst))

responses <- unique(responses)

##Check that there is only 1 response per reviewer
responses %>%
  group_by(Path, analyst, tiebreaker_reviewed) %>%
  tally -> response_count

if(nrow(subset(response_count, n>1)) >0 ){
  print("Multiple entries for the same reviewer, need to fix")
}else{
  print("No repeats")
}

#Reformat data
responses %>%
  select(c(Path, analyst, Is_de_novo, Reason_unsure_follow_up, Notes, tiebreaker_reviewed)) %>%
  # select(c(Path, analyst, `Is de novo`)) %>%
  pivot_wider(names_from = analyst, values_from = Is_de_novo) -> responses2

reviewers_names <- unique(flipbook_metadata$V1)

##Prepare Path_fix for later
responses$Path_fix <- responses$Path
responses[grepl("CPX", responses$Path),]$Path_fix <- sapply(responses[grepl("CPX", responses$Path),]$Path, function(x) paste0(paste(strsplit(x, "_")[[1]][1:as.numeric(length(strsplit(x, "_")[[1]])-2)], collapse = "_"), "_denovo.png") )

#Add expression for upsetR plot
responses2$result <- apply(responses2, 1, function(row){
  result <- row[reviewers_names][!is.na(row[reviewers_names])]
  if (length(result) >1){
    paste(result, collapse = "&")
  }else{result}
})

#Add analysts
responses2$reviewers <- apply(responses2, 1, function(row){
  paste(names(row[reviewers_names][!is.na(row[reviewers_names])]), collapse = "&")
})

responses2$num_reviewers <- apply(responses2, 1, function(row){
  length(names(row[reviewers_names][!is.na(row[reviewers_names])]))
})

#Make entries unique

##If it fails, detach packages and load tidyverse only:
# detachAllPackages()
# library(tidyverse)
# library(dplyr)

responses2 %>%
  group_by(Path) %>%
  summarise(Path = unique(Path),
            result = paste(result, collapse = "&"),
            Reason_unsure_follow_up = paste(unique(Reason_unsure_follow_up), collapse = ";"),
            Notes = paste(unique(Notes), collapse = ";"),
            reviewers = paste(unique(reviewers), collapse = "&"),
            tiebreaker_reviewed = any(tiebreaker_reviewed),
            num_reviewers = paste(unique(num_reviewers), collapse = ";")
            ) -> responses2_ref

#Fix num reviewers
responses2_ref$num_reviewers <- apply(responses2_ref, 1, function(row){
  length(strsplit(as.character(row["reviewers"]), "&")[[1]])
})

#Keep only those with at least one reviewer
responses3 <- subset(responses2_ref, !is.na(num_reviewers))

#Fix CPX responses
responses_cpx <- responses3[grepl("CPX", responses3$Path),]

if(nrow(responses_cpx) >2){
  # responses_cpx$name <- sapply(responses_cpx$Path, function(x) paste(strsplit(x, "_")[[1]][1:as.numeric(length(strsplit(x, "_")[[1]])-2)], collapse = "_"))
  # responses_cpx$name <- sapply(responses_cpx$Path, function(x) paste(strsplit(x, "_")[[1]][1:as.numeric(length(strsplit(x, "_")[[1]])-2)], collapse = "_"))

  responses_cpx %>% #I had lots of issues with the library loading
  # subset(responses_cpx, name == "__202030105__591583_second_run_cp_cohort.chr10.final_cleanup_CPX_chr10") %>%
# rbind(subset(responses_cpx, name == "__202030105__591583_second_run_cp_cohort.chr10.final_cleanup_CPX_chr10"), head(subset(responses_cpx, name != "__202030105__591583_second_run_cp_cohort.chr10.final_cleanup_CPX_chr10"), 100)) %>%
    group_by(Path) %>%
    # group_by(name) %>%
    summarise(Path = unique(Path),
    # summarise(Path = paste0(unique(name), "_denovo.png"),
              result = paste(unique(result), collapse = ";"),
              Reason_unsure_follow_up = paste(unique(Reason_unsure_follow_up), collapse = ";"),
              Notes = paste(unique(Notes), collapse = ";"),
              reviewers = paste(unique(reviewers), collapse = ";"),
              tiebreaker_reviewed = any(tiebreaker_reviewed),
              num_reviewers = paste(unique(num_reviewers), collapse = ";")
              ) %>%
    select(Path, result, Reason_unsure_follow_up, Notes, reviewers, tiebreaker_reviewed, num_reviewers) -> responses_cpx2

  responses_cpx2$Path <- sapply(responses_cpx2$Path, function(x) paste0(paste(strsplit(x, "_")[[1]][1:as.numeric(length(strsplit(x, "_")[[1]])-2)], collapse = "_"), "_denovo.png") )

  responses_non_cpx <- subset(responses3, !grepl("CPX", Path), select = c(Path, result, Reason_unsure_follow_up, Notes, reviewers, tiebreaker_reviewed, num_reviewers))

  responses_fix <- rbind(responses_cpx2, responses_non_cpx)

}else{
  responses_fix <- responses3
}

##Add discrepant result conclusion
tiebreaker_results <- unique(subset(responses, tiebreaker_reviewed == TRUE, select = c(Path_fix, analyst, Is_de_novo)))
names(tiebreaker_results) <- c("Path", "tiebreaker_reviewer", "tiebreaker_result")

responses_fix_disc <- merge(responses_fix, tiebreaker_results, by = "Path", all.x = T, all.y = F)

#Append responses to denovo table
denovo_merge <- merge(denovo, responses_fix_disc, by = "Path", all.x = T, all.y = F)

#Create discrepant results columns
denovo_merge$SV_plots_reviewed <- !is.na(denovo_merge$num_reviewers)

#Add final result column for those that have been reveiwed by 2 and are not discrepant
# + those with tiebreaker
denovo_merge$result_final <- ""

denovo_merge$result_final <- apply(denovo_merge, 1, function(row){
  ifelse(all(strsplit(row['result'], "&|;")[[1]] == "yes"), "yes", row['result_final'])
})

denovo_merge$result_final <- apply(denovo_merge, 1, function(row){
  ifelse(all(strsplit(row['result'], "&|;")[[1]] %in% c("no", "unsure")), "no_unsure", row['result_final'])
})

denovo_merge[!is.na(denovo_merge$tiebreaker_result),]$result_final <- denovo_merge[!is.na(denovo_merge$tiebreaker_result),]$tiebreaker_result
denovo_merge[denovo_merge$result_final %in% c("no", "unsure"),]$result_final <- "no_unsure"

# denovo_merge[!denovo_merge$all_yes & !denovo_merge$all_unsure_no,]$result_discrepant

#Flag if batch is from INS_recover and fix batch/cohort columns
#Fix batch for plot
denovo_merge$is_ins_recover <- denovo_merge$batch == "INS_recover"
denovo_merge[is.na(denovo_merge$cohort),]$cohort <- denovo_merge[is.na(denovo_merge$cohort),]$batch
denovo_merge[denovo_merge$batch == "INS_recover",]$batch <- denovo_merge[denovo_merge$batch == "INS_recover",]$cohort



##Flag if sample has an issue
samples_parental_error <- subset(qc_file, maternal_error | paternal_error)$sample_id
samples_parental_error_ids <- subset(ped, snv_pipeline_id %in% samples_parental_error)$sv_pipeline_id

samples_sex_error <- subset(qc_file, sex_error)$sample_id
samples_sex_error_ids <- subset(ped, snv_pipeline_id %in% samples_sex_error)$sv_pipeline_id


denovo_merge$relatedness_error <- denovo_merge$sample %in% samples_parental_error_ids
denovo_merge$sex_error <- denovo_merge$sample %in% samples_sex_error_ids

denovo_merge$any_error <- denovo_merge$relatedness_error | denovo_merge$sex_error


##Write final output
write.table(denovo_merge, out_file, sep = "\t", quote = F, row.names = F)

denovo_merge$chrom_sample_type <- paste0(denovo_merge$chrom, "_", denovo_merge$sample, "_", denovo_merge$SVTYPE)
denovo_for_merging <- subset(denovo_merge, select = c(chrom_sample_type, start, end, sample_name))
write.table(denovo_for_merging, paste0(gsub(".txt", "", out_file), "_ref_for_merge.txt"), sep = "\t", quote = F, row.names = F, col.names = F)

#Without taking type into account
# denovo_merge$chrom_sample <- paste0(denovo_merge$chrom, "_", denovo_merge$sample)
# denovo_for_merging <- subset(denovo_merge, select = c(chrom_sample, start, end, sample_name))
# write.table(denovo_for_merging, paste0(gsub(".txt", "", out_file), "_ref_for_merge.txt"), sep = "\t", quote = F, row.names = F, col.names = F)


##WORK IN PROGRESS
# ##Merge fragmented calls
# denovo_merge$chrom_sample <- paste0(denovo_merge$chrom, "_", denovo_merge$sample)
# denovo_merge_coords$chrom_sample <- gsub("-", "_", denovo_merge$chrom_sample)
# denovo_merge$region <- paste0(denovo_merge$chrom, ":", denovo_merge$start, "-", denovo_merge$end)
#
# denovo_merge_coords <- subset(denovo_merge, select = c(chrom, start, end, sample_name))
# names(denovo_merge_coords)[1] <- "chr"
#
# denovo_merge_coords_sort <- bedr.sort.region(denovo_merge_coords, method = "natural")
#
# denovo_merge_coords_sort <- bedr(
#   engine = "bedtools",
#   input = list(i = denovo_merge_coords),
#   method = "sort",
#   params = ""
# )

##Write real cleaned only
# write.table(denovo_merge, out_file, sep = "\t", quote = F, row.names = F)
