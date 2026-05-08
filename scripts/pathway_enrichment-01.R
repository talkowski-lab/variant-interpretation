#!/usr/bin/env Rscript

# Generate annotated snv_indels_info file with phenotype counts
# Author: Alba Sanchis-Juan

#############
##LIBRARIES##
#############
library(data.table)
library(purrr)

#############
##ARGUMENTS##
#############
option_list = list(
  make_option(c("-i", "--info"), type="character", default=".",
              help="snvs indels info file", metavar="character"),
  make_option(c("-p", "--pedigree"), type="character", default=".",
              help="final pedigree with controls", metavar="character"),
  make_option(c("-d", "--denovo"), type="character", default=".",
              help="de novo file", metavar="character"),
  make_option(c("-a", "--phenotype1"), type="character", default=".",
              help="phenotype1 [default= %default]", metavar="character"),
  make_option(c("-b", "--phenotype2"), type="character", default=".",
              help="phenotype2 [default= %default]", metavar="character"),
)

opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt <- parse_args(opt_parser)

snvs_indels_info <- fread(opt$info)
pedigree_final <- fread(opt$pedigree)
dn_snvs_indels <- fread(opt$denovo)
phenotype1 <- opt$phenotype1
phenotype2 <- opt$phenotype2

###Hpo maker phenotypes
snvs_indels_info_short <- subset(snvs_indels_info, select = c("gene_id", "gene", grep("prior", names(snvs_indels_info), value = T), "LOEUF", grep("mu.", names(snvs_indels_info), value = T), "chrom", "start", "end", "symbol_mart", "hgnc_id", "gene_type", "gene_id_alt"))

##Marker phenos
# dn_snvs_indels[dn_snvs_indels$sample %in% subset(pedigree_final, affected == 1)$snv_pipeline_id,]$hpo_marker_name <- "Control"
make_suffix <- function(x) {
  x %>%
    tolower() %>%
    gsub("[^a-z0-9]+", "_", .) %>%  # replace non-alphanumeric with _
    gsub("^_|_$", "", .)             # strip leading/trailing underscores
}

#keep only detailed phenotypes or controls
pedigree_final <- subset(pedigree_final, affected == 1 | phenotype_level == "detailed")

pedigree_final$phenotype1 <- grepl(phenotype1, pedigree_final$hpo_marker_name)
pedigree_final$phenotype2 <- grepl(phenotype2, pedigree_final$hpo_marker_name)

pedigree_final_unique <- subset(pedigree_final, (phenotype1 & !phenotype2) | (!phenotype1 & phenotype2))

hpo_list <- c(phenotype1, phenotype2)

# For each phenotype, compute counts and rename columns with suffix
hpo_tables <- map(hpo_list, function(hpo) {
  suffix <- make_suffix(hpo)
  
  dn_snvs_indels %>%
    # subset(grepl(hpo, hpo_marker_name)) %>%          # <-- filter to this phenotype
    subset(sample %in% subset(pedigree_final_unique, grepl(hpo, hpo_marker_name))$snv_pipeline_id) %>%
    separate_rows(symbol_fix, sep = "\\|") %>%
    subset(symbol_fix != "") %>%
    distinct(symbol_fix, GENE_SAMP, CSQ_highest, mis1, mis2) %>%
    group_by(symbol_fix) %>%
    summarise(
      dn.ptv  = sum(CSQ_highest == "PTV"),
      dn.mis  = sum(CSQ_highest == "MIS"),
      dn.syn  = sum(CSQ_highest == "SYN"),
      dn.mis1 = sum(mis1, na.rm = TRUE),
      dn.mis2 = sum(mis2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename_with(~ paste0(., ".", suffix), -symbol_fix)  # add suffix to all cols except symbol_fix
})

# Join all phenotype tables into one wide matrix
result <- purrr::reduce(hpo_tables, full_join, by = "symbol_fix") %>%
  mutate(across(-symbol_fix, ~ replace_na(., 0)))

#Merge
snvs_indels_info_counts <- merge(snvs_indels_info_short, 
                                 result, 
                                 by.x = "gene", 
                                 by.y="symbol_fix", 
                                 all.x=TRUE, 
                                 all.y=FALSE)

##create table with correspondance and counts
marker_corresp <- data.table("hpo_marker_name" = hpo_list, "hpo_suffix" = make_suffix(hpo_list))

marker_corresp$num_samples <- apply(marker_corresp, 1, function(row){
  nrow(subset(pedigree_final_unique, grepl(row['hpo_marker_name'], hpo_marker_name)))
})

marker_corresp_with_controls <- rbind(marker_corresp, 
                                      data.table("hpo_marker_name" = "Control", 
                                                 "hpo_suffix" = "control", 
                                                 "num_samples" = nrow(subset(pedigree_final, affected ==1))))

#save table
write.table(marker_corresp, "hpo_marker_pair_counts.txt",
            sep = "\t", quote = F, row.names = F)

write.table(snvs_indels_info_counts, "crossdev_pair_dn_counts.txt",
            sep = "\t", quote = F, row.names = F)

# write.table(subset(snvs_indels_info_counts, gene_type == "protein_coding"), 
#             "crossdev_master_dn_counts.coding_only.txt",
#             sep = "\t", quote = F, row.names = F)