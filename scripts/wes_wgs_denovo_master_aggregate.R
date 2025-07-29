#!/usr/bin/env Rscript

##Merge WES and WGS callsets

#Load libraries
library(data.table)
library(optparse)
library(parallel)

# Define Input Arguments
option_list = list(
  make_option(c("-g", "--snvgenomes"), type="character", default=NULL,
              help="SNV/indels WGS file", metavar="character"),
  make_option(c("-e", "--snvexomes"), type="character", default=NULL,
              help="SNV/indels WES file", metavar="character"),
  make_option(c("-c", "--cnvexomes"), type="character", default=NULL,
              help="CNV WES file", metavar="character"),
  make_option(c("-s", "--svgenomes"), type="character", default=NULL,
              help="SV WGS file", metavar="character"),
  make_option(c("-p", "--ped"), type="character", default=NULL,
              help="Pedigree file", metavar="character"),
  make_option(c("-d", "--date"), type="character", default=NULL,
              help="Release date", metavar="character"),
  make_option(c("-a", "--additional"), type="character", default=NULL,
              help="Additional/DDD SNV/indels denovo callset", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

#Define input parameteres and read files
snvgenomes_file <- opt$snvgenomes
snvexomes_file <- opt$snvexomes
cnvexomes_file <- opt$cnvexomes
svgenomes_file <- opt$svgenomes
ped_file <- opt$ped
release_date <- opt$date
additional_file <- opt$additional

snv_wgs <- fread(snvgenomes_file)
snv_wes <- fread(snvexomes_file)
cnv_wes <- fread(cnvexomes_file)
sv_wgs <- fread(svgenomes_file)
ped <- fread(ped_file)
additional_denovo <- fread(additional_file)

#############################
##Merge WES+WGS SNVs/indels##
#############################

##Reformat SNV WES and exclude outliers
names(snv_wes)[names(snv_wes) == "id"] <- "SAMPLE"
names(snv_wes)[!names(snv_wes) %in% names(snv_wgs)] <- paste0("wes_", names(snv_wes)[!names(snv_wes) %in% names(snv_wgs)])
snv_wes$data_type <- "WES"

##Exclude SNV Outlier samples - no lower bound for WES
snv_wes_q0_5 <- quantile(plyr::count(snv_wes$SAMPLE)$freq, 0.95) #Keeping quartile below which 95% of the data falls
snv_wes_iqr_value <- IQR(plyr::count(snv_wes$SAMPLE)$freq)
snv_wes_upper_bound <- snv_wes_q0_5 + 1.5 * snv_wes_iqr_value

snv_wes_outliers <- plyr::count(snv_wes$SAMPLE)[plyr::count(snv_wes$SAMPLE)$freq > snv_wes_upper_bound, ]$x
snv_wes_samples_keep <- unique(snv_wes$SAMPLE)[!unique(snv_wes$SAMPLE) %in% snv_wes_outliers]

##Reformat SNV WGS and exclude outliers
names(snv_wgs)[!names(snv_wgs) %in% names(snv_wes)] <- paste0("wgs_", names(snv_wgs)[!names(snv_wgs) %in% names(snv_wes)])
snv_wgs$data_type <- "WGS"

##Exclude SNV Outlier samples for upper and lower bond
snv_wgs_q0_5 <- quantile(plyr::count(snv_wgs$SAMPLE)$freq, 0.95) #Keeping quartile below which 95% of the data falls
snv_wgs_iqr_value <- IQR(plyr::count(snv_wgs$SAMPLE)$freq)
snv_wgs_upper_bound <- snv_wgs_q0_5 + (1.5 * snv_wgs_iqr_value)
snv_wgs_lower_bound <- snv_wgs_q0_5 - (3 * snv_wgs_iqr_value)

snv_wgs_outliers <- plyr::count(snv_wgs$SAMPLE)[plyr::count(snv_wgs$SAMPLE)$freq > snv_wgs_upper_bound | plyr::count(snv_wgs$SAMPLE)$freq < snv_wgs_lower_bound, ]$x
snv_wgs_samples_keep <- unique(snv_wgs$SAMPLE)[!unique(snv_wgs$SAMPLE) %in% snv_wgs_outliers]

snv_all <- rbind(snv_wes, snv_wgs, fill = TRUE)

##Add additional denovo calls from papers
names(additional_denovo)[2] <- "SAMPLE"
# additional_denovo$cohort <- "DDD"
# additional_denovo2 <- subset(additional_denovo, IN_KAPLANIS == TRUE)

names(additional_denovo)[!names(additional_denovo) %in% names(snv_all)] <- paste0("wes_additional_", names(additional_denovo)[!names(additional_denovo) %in% names(snv_all)])
additional_denovo$data_type <- "WES"

snv_all_additional <- rbind(snv_all, additional_denovo, fill = TRUE)


#####################
##Merge WES+WGS SVs##
#####################

##Reformat CNV WES and exclude outliers
names(cnv_wes)[names(cnv_wes) == "sample"] <- "SAMPLE"
cnv_wes$name_fix <- cnv_wes$name
names(sv_wgs)[names(sv_wgs) == "sample"] <- "SAMPLE"

names(cnv_wes)[!names(cnv_wes) %in% names(sv_wgs)] <- paste0("wes_", names(cnv_wes)[!names(cnv_wes) %in% names(sv_wgs)])
cnv_wes$data_type <- "WES"

##Exclude CNV Outlier samples - no lower bound
cnv_wes_q0_5 <- quantile(plyr::count(cnv_wes$SAMPLE)$freq, 0.995) #Keeping quartile below which 99.5% of the data falls
cnv_wes_iqr_value <- IQR(plyr::count(cnv_wes$SAMPLE)$freq)
cnv_wes_upper_bound <- cnv_wes_q0_5 + 1.5 * cnv_wes_iqr_value

cnv_wes_outliers <- plyr::count(cnv_wes$SAMPLE)[plyr::count(cnv_wes$SAMPLE)$freq > cnv_wes_upper_bound, ]$x
cnv_wes_samples_keep <- unique(cnv_wes$SAMPLE)[!unique(cnv_wes$SAMPLE) %in% cnv_wes_outliers]


##Reformat SV WGS and exclude outliers
##Exclude AF>0.01 for WGS de novo calls
# sv_wgs <- subset(sv_wgs, !(AF > 0.01 & is_GD == FALSE))
names(sv_wgs)[!names(sv_wgs) %in% names(cnv_wes)] <- paste0("wgs_", names(sv_wgs)[!names(sv_wgs) %in% names(cnv_wes)])
sv_wgs$data_type <- "WGS"

##Exclude sv Outlier samples - no lower bound
sv_wgs_q0_5 <- quantile(plyr::count(sv_wgs$SAMPLE)$freq, 0.995) #Keeping quartile below which 99.5% of the data falls
sv_wgs_iqr_value <- IQR(plyr::count(sv_wgs$SAMPLE)$freq)
sv_wgs_upper_bound <- sv_wgs_q0_5 + 1.5 * sv_wgs_iqr_value
sv_wgs_lower_bound <- sv_wgs_q0_5 - (3 * sv_wgs_iqr_value) #Pending to check the value and continue from here to update callset

sv_wgs_outliers <- plyr::count(sv_wgs$SAMPLE)[plyr::count(sv_wgs$SAMPLE)$freq > sv_wgs_upper_bound | plyr::count(sv_wgs$SAMPLE)$freq > sv_wgs_lower_bound, ]$x
sv_wgs_samples_keep <- unique(sv_wgs$SAMPLE)[!unique(sv_wgs$SAMPLE) %in% sv_wgs_outliers]

##MERGE FILES
sv_all <- rbind(cnv_wes, sv_wgs, fill = TRUE)

##Remove outliers SNVs/indels and SVs in samples with any error = TRUE
# sv_all$KEEP <- ((sv_all$SAMPLE %in% c(subset(ped, v01_freeze == TRUE)$sv_pipeline_id)) &
#                   (!sv_all$SAMPLE %in% c(cnv_wes_outliers, sv_wgs_outliers)))
# snv_all_additional$KEEP <- ((snv_all_additional$SAMPLE %in% c(subset(ped, v01_freeze == TRUE)$snv_pipeline_id)) &
#                        (!snv_all_additional$SAMPLE %in%  c(snv_wes_outliers, snv_wgs_outliers)))
sv_all$KEEP <- !sv_all$SAMPLE %in% c(cnv_wes_outliers, sv_wgs_outliers)
snv_all_additional$KEEP <- !snv_all_additional$SAMPLE %in%  c(snv_wes_outliers, snv_wgs_outliers)

#Fix missing IDS in SNV data for additional cohorts
snv_all_additional[snv_all_additional$SAMPLE == "",]$SAMPLE <- snv_all_additional[snv_all_additional$SAMPLE == "",]$`wes_additional_Blinded ID`
snv_all_additional[snv_all_additional$SAMPLE == "",]$SAMPLE <- snv_all_additional[snv_all_additional$SAMPLE == "",]$wes_additional_id

##snv_all_additional[snv_all_additional$cohort == "DDD_Kaplanis_et_al_2020_GRCh38_liftover_noGeneDxASC",]$SAMPLE <- snv_all_additional[snv_all_additional$cohort == "DDD_Kaplanis_et_al_2020_GRCh38_liftover_noGeneDxASC",]$wes_additional_id

##Get final callset
snv_all_final <- subset(snv_all_additional, KEEP == TRUE)
sv_all_final <- subset(sv_all, KEEP == TRUE)

##Write outputs
write.table(snv_all_additional, gzfile(paste0("DENOVO_SNVS_INDELS-", release_date, ".txt.gz"), "w"), sep = "\t", quote = F, row.names = F)
close(gzfile(paste0("DENOVO_SNVS_INDELS-", release_date, ".txt.gz"), "w"))

write.table(sv_all, gzfile(paste0("DENOVO_SVS-", release_date, ".txt.gz"), "w"), sep = "\t", quote = F, row.names = F)
close(gzfile(paste0("DENOVO_SVS-", release_date, ".txt.gz"), "w"))

write.table(snv_all_final, gzfile(paste0("DENOVO_SNVS_INDELS-", release_date, ".final.txt.gz"), "w"), sep = "\t", quote = F, row.names = F)
close(gzfile(paste0("DENOVO_SNVS_INDELS-", release_date, ".final.txt.gz"), "w"))

write.table(sv_all_final, gzfile(paste0("DENOVO_SVS-", release_date, ".final.txt.gz"), "w"), sep = "\t", quote = F, row.names = F)
close(gzfile(paste0("DENOVO_SVS-", release_date, ".final.txt.gz"), "w"))

write.table(c(snv_wes_outliers, snv_wgs_outliers), gzfile(paste0("DENOVO_SNVS_INDELS-", release_date, ".outliers.txt.gz"), "w"), sep = "\n", quote = F, row.names = F, col.names = F)
close(gzfile(paste0("DENOVO_SNVS_INDELS-", release_date, ".outliers.txt.gz"), "w"))

write.table(c(cnv_wes_outliers, sv_wgs_outliers), gzfile(paste0("DENOVO_SVS-", release_date, ".outliers.txt.gz"), "w"), sep = "\n", quote = F, row.names = F, col.names = F)
close(gzfile(paste0("DENOVO_SVS-", release_date, ".outliers.txt.gz"), "w"))