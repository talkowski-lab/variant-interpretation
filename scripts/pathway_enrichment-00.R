#!/usr/bin/env Rscript

# CrossDEV Pathway Analysis
# Author: Nehir E. Kurtas

################
##LOAD LIBRARY##
################
library(optparse)
library(dplyr)
library(stringr)
library(purrr)
library(AnnotationDbi)
library(org.Hs.eg.db)


#############
##ARGUMENTS##
#############
option_list = list(
  make_option(c("-g", "--genes"), type="character", default=".",
              help="list of genes, inluding mutation and phenotype table, to run pathway on [default= %default]", metavar="character"),
  make_option(c("-i", "--id_hgnc"), type="character", default=".",
              help="HGNC id file", metavar="character")

  )


opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt <- parse_args(opt_parser)


info_format <- read.delim(opt$genes, header = TRUE, sep = "\t")

#UPLOAD hgnc database
hgnc_ids <- read.delim(opt$id_hgnc, header = TRUE, sep = "\t")

info_select <- info_format %>% dplyr::select(gene,gene_id,LOEUF)

#name <- sub("\\.txt$", "", basename(opt$genes))


################
#GET ENTREZ ID #
################

#1. FIRST Method: USE ENSEMBL ID in the info_select_updated
background_gene_ids <- mapIds(org.Hs.eg.db, keys = info_select$gene_id, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
background_gene_ids_df_first <- data.frame( ensembl_id = names(background_gene_ids),
                                            entrez_id = as.character(background_gene_ids) # gene IDs will be the values
)
#to get gene symbols
info_select_updated_short <- info_select %>% dplyr::select(gene,gene_id)
background_gene_ids_df_first_merged <- merge(background_gene_ids_df_first,info_select_updated_short,by.x="ensembl_id",by.y="gene_id")

#however still a lot of genes are missing from mapping to extract entrez id
background_gene_names_NAs <- subset(info_select, gene_id %in% subset(background_gene_ids_df_first,is.na(entrez_id)==TRUE)$ensembl_id)$gene


#2. SECOND Method: USE GENE SYMBOL in info_select_updated and hgnc_ids
background_gene_ids_df_second <- data.frame(
  gene_name = subset(hgnc_ids, symbol %in% background_gene_names_NAs)$symbol,  # gene names will be the row names
  entrez_id_second = subset(hgnc_ids, symbol %in% background_gene_names_NAs)$entrez_id  # gene IDs will be the values
)

#merge with first entrez id map
background_gene_ids_df_first_second_merged <- merge(background_gene_ids_df_first_merged,background_gene_ids_df_second,by.y="gene_name",by.x="gene",all.x=TRUE)
background_gene_ids_df_first_second_merged$combined_entrez_ids <- coalesce(as.character(background_gene_ids_df_first_second_merged$entrez_id),
                                                                           as.character(background_gene_ids_df_first_second_merged$entrez_id_second))


#3. THIRD Methods: use previous symbols
background_gene_ids_df_first_second_merged_NA <- subset(background_gene_ids_df_first_second_merged,is.na(combined_entrez_ids))$gene

hgnc_ids <- hgnc_ids %>%
  mutate(matched_gene = map_chr(prev_symbol, function(symbols) {
    matched <- background_gene_ids_df_first_second_merged_NA[str_detect(symbols, paste0("(^|\\|)", background_gene_ids_df_first_second_merged_NA, "(\\||$)"))]
    if (length(matched) > 0) matched[1] else FALSE
  }))


background_gene_ids_df_third <- data.frame(
  gene_name = subset(hgnc_ids, matched_gene %in% background_gene_ids_df_first_second_merged_NA)$matched_gene,  # gene names will be the row names
  entrez_id_third = subset(hgnc_ids, matched_gene %in% background_gene_ids_df_first_second_merged_NA)$entrez_id  # gene IDs will be the values
)

background_gene_ids_df_first_second_third_merged <- merge(background_gene_ids_df_first_second_merged,background_gene_ids_df_third,by.y="gene_name",by.x="gene",all.x=TRUE)

background_gene_ids_df_first_second_third_merged$combined_entrez_ids_second <- coalesce(as.character(background_gene_ids_df_first_second_third_merged$combined_entrez_ids),
                                                                                        as.character(background_gene_ids_df_first_second_third_merged$entrez_id_third))


background_gene_ids_df <- background_gene_ids_df_first_second_third_merged %>% dplyr::select(gene,combined_entrez_ids_second)


#REMOVE DUPLICATES
info_select_updated_entrez_ids <- merge(info_select,background_gene_ids_df,by.x="gene",by.y="gene",all.x=TRUE)

info_select_2 <- info_select_updated_entrez_ids %>% dplyr::select(gene,combined_entrez_ids_second)
colnames(info_select_2) <-c("gene","entrez_id")

info_select_updated <- merge(info_format,info_select_2,by.x="gene",by.y="gene")

write.table(info_select_updated,
            file =  "genelist_fixed_entrez_id.txt", 
            sep = "\t",
            row.names = FALSE,col.names = TRUE,quote=FALSE)



