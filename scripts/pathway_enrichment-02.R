#!/usr/bin/env Rscript

# CrossDEV Pathway Analysis
# Author: Nehir E. Kurtas

################
##LOAD LIBRARY##
################
library(optparse)
library(dplyr)
library(tidyr)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GO.db)
library(ggplot2)
library(ggrepel)


#############
##ARGUMENTS##
#############
option_list = list(
  make_option(c("-m", "--marker"), type="character", default=".",
              help="marker counts file", metavar="character"),
  make_option(c("-g", "--genes"), type="character", default=".",
              help="list of genes, inluding mutation and phenotype table, to run pathway on [default= %default]", metavar="character"),
  make_option(c("-a", "--phenotype1"), type="character", default=".",
              help="phenotype1 [default= %default]", metavar="character"),
  make_option(c("-b", "--phenotype2"), type="character", default=".",
              help="phenotype2 [default= %default]", metavar="character"),
  make_option(c("-u", "--mutation"), type="character", default=".",
              help="choose dn.syn or dn.mis2.ptv", metavar="character"),
  make_option(c("-e", "--eigenvalue"), type="character", default=".",
              help="eigen value used to define significance threshold in volcano plots", metavar="character")
)


opt_parser <- OptionParser(option_list=option_list, add_help_option=FALSE)
opt <- parse_args(opt_parser)

genes <- read.delim(opt$genes, sep = "\t")
counts_df <- read.delim(opt$marker, sep = "\t")
phenotype1 <- opt$phenotype1
phenotype2 <- opt$phenotype2
mutation <- opt$mutation
eigenvalue <- opt$eigenvalue



#Upload input file
#This file "crossdev_master_dn_counts.v4_entrez_id.txt" is created by "09_Pathway_step1.R"
# genes <- read.delim("crossdev_master_dn_counts.v4_entrez_id.txt", sep = "\t")
# counts_df <- read.delim("hpo_marker_counts.txt", header = TRUE, sep = "\t")
# phenotype1 <- "autism"
# phenotype2 <- "abnormality_of_the_face"
# mutation <- "dn.mis2.ptv"


######################
#GET ALL GO DATABASE#
######################
all_go <- keys(GO.db, keytype = "GOID")
# Get the ontology category (BP, MF, CC)
go_ont <- AnnotationDbi::select(GO.db, keys = all_go, columns = c("TERM", "ONTOLOGY"), keytype = "GOID")

#upload functions
convert_to_gene_symbols <- function(entrez_ids) {
  # Split the Entrez IDs by "/"
  entrez_ids_split <- strsplit(entrez_ids, "/")[[1]]

  # Map the Entrez IDs to gene symbols
  gene_symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids_split, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

  # Combine the gene symbols back with "/"
  gene_names <- paste(gene_symbols, collapse = "/")

  return(gene_names)
}


collect_stats <-function(data,
                         phenotype,
                         mutation){


  # data <- file_subset
  # phenotype <-"abnormality_of_the_nervous_system"
  # mutation <-"dn.mis2.ptv"

  #go_ont
  perform_fisher_test <- function(row) {
    matrix_data_i <- matrix(c(as.numeric(row$success),
                              as.numeric(row$annotated)-as.numeric(row$success),
                              as.numeric(row$mutated)-as.numeric(row$success),
                              as.numeric(row$background)-as.numeric(row$mutated)-as.numeric(row$annotated)+as.numeric(row$success)),  # Replace with your own data
                            nrow = 2,
                            byrow = TRUE,
                            dimnames = list(
                              c("Annotated", "Not_annotated"),
                              c("mutated", "not_mutated")
                            ))

    fisher_test_result <- fisher.test(matrix_data_i,alternative="greater")
    test_result <- data.frame(fisher_test_result$estimate[[1]],fisher_test_result$conf.int[[1]],fisher_test_result$conf.int[[2]],fisher_test_result$p.value)
    test_result$type_new <- "fisher.exact.greater"
    colnames(test_result) <- c("OR_Fisher","CI_low_Fisher","CI_up_Fisher","p_value_Fisher","type")
    return(test_result)
  }
  gene_ids_df <- data.frame(
    ENSEML = data$gene_id,  # gene names will be the row names
    gene_id = as.character(data$entrez_id)  # gene IDs will be the values
  )
  ######################
  #Biological Process###
  ######################
  print("calculating biological process")

  #all categories
  go_ont_BP <- go_ont %>% filter (ONTOLOGY=="BP")

  GO_BP <- enrichGO(gene = unique(gene_ids_df$gene_id),
                    OrgDb = org.Hs.eg.db,    # For human genes (replace for other organisms)
                    universe=unique(as.character(subset(genes,is.na(entrez_id)==FALSE)$entrez_id)),
                    keyType = "ENTREZID",    # Specify the gene ID type used
                    ont = "BP",              # You can also use "MF" or "CC" for Molecular Function or Cellular Component
                    pAdjustMethod = "bonferroni",    # Adjust p-values using the Benjamini-Hochberg method
                    pvalueCutoff = 1,     # Adjust the p-value cutoff threshold
                    qvalueCutoff = 1)
  GO_BP_df <- as.data.frame(GO_BP)
  GO_BP_df$GO <- "Biological Process"
  GO_BP_df[c('success', 'mutated')] <- str_split_fixed(GO_BP_df$GeneRatio, '/', 2)
  GO_BP_df[c('annotated', 'background')] <- str_split_fixed(GO_BP_df$BgRatio, '/', 2)

  #re-calculate p-values
  results_list <- lapply(seq_len(nrow(GO_BP_df)), function(i) perform_fisher_test(GO_BP_df[i, ]))
  final_results <- do.call(rbind, results_list)
  GO_BP_df <- cbind(GO_BP_df,final_results)
  GO_BP_df$p_fisher_corrected <- p.adjust(GO_BP_df$p_value_Fisher, method = "bonferroni", n = nrow(go_ont_BP))

  GO_BP_df$Number_of_GO_categories <- nrow(GO_BP_df)

  #extract gene name for significant results
  # GO_BP_df_select <- GO_BP_df %>% filter(p_value_Fisher <0.05)
  # GO_BP_df_select$gene_name<- sapply(GO_BP_df_select$geneID, convert_to_gene_symbols)
  # GO_BP_df_select2 <- GO_BP_df_select %>% dplyr::select(ID,gene_name)
  # GO_BP_df <- merge(GO_BP_df,GO_BP_df_select2,by.x="ID",by.y="ID",all.x=TRUE)

  #add the remaining BP categories that were missing due to zero intersect
  number_of_mutations <- as.numeric(strsplit(GO_BP_df$GeneRatio[1], "/")[[1]][2])
  go_ont_BP_missing <- go_ont_BP %>% filter (!GOID %in% unique(GO_BP_df$ID))
  missing_BP <- data.frame(ID = go_ont_BP_missing$GOID,
                           Description = go_ont_BP_missing$TERM,
                           GeneRatio= paste("0",number_of_mutations,sep="/"),
                           GO="Biological Process",
                           success=as.numeric("0"),
                           mutated=number_of_mutations)

  missing_cols <- setdiff(names(GO_BP_df), names(missing_BP))
  missing_BP[missing_cols] <- NA
  missing_BP <- missing_BP[, names(GO_BP_df)]

  #merge tested and missing GO terms
  GO_BP_df_final <- rbind(GO_BP_df, missing_BP)


  GO_BP_df_final$mutation <- mutation
  GO_BP_df_final$phenotype <- phenotype
  return(GO_BP_df_final)
}

#Repeat this in all OFC categories
collect_stats_catgories <- function(data,
                                    mutation,
                                    phenotype){

  # data <- info_select_updated
  # phenotype <-"abnormality_of_the_nervous_system"
  # mutation <-"dn.mis2.ptv"
  # file_subset <- info_select_updated %>% filter (get(mutation_name)>1 | get(mutation_name)==1)

  mutation_name <- paste0(mutation, ".",phenotype, sep="")

  file_subset <- data %>% filter (get(mutation_name)>1 | get(mutation_name)==1)
  file_stat <- collect_stats(file_subset,
                             mutation=mutation,
                             phenotype=phenotype)
  return(file_stat)
}



phenotypes <- counts_df$hpo_suffix
for (pheno in phenotypes) {
  cols_to_sum <- c(paste0("dn.ptv.", pheno),
                   paste0("dn.mis2.", pheno))
  # Create new column with the sum
  genes[[paste0("dn.mis2.ptv.", pheno)]] <- rowSums(genes[, cols_to_sum], na.rm = TRUE)
  genes[[paste0("dn.syn.", pheno)]] <-
    ifelse(is.na(genes[[paste0("dn.syn.", pheno)]]), 0,
           genes[[paste0("dn.syn.", pheno)]])

}


######################
# PATHWAY ENRICHMENT##
######################
print("Step 1/5 is processing")
#phenotype1
GO_phenotype1 <- collect_stats_catgories(genes,
                                                  mutation=mutation,
                                                  phenotype=phenotype1)


#phenotype2
GO_phenotype2 <- collect_stats_catgories(genes,
                                                  mutation=mutation,
                                                  phenotype=phenotype2)

########################################
#COMBINE  phenotype1 and phenotype2 ####
########################################

merge_pheno1_pheno2 <-function (pheno1_input,
                                pheno2_input
) {

  # pheno1_input <-GO_phenotype1_mis2_ptv
  # pheno2_input <- GO_phenotype2_mis2_ptv
  #
  pheno2_input_simple <- pheno2_input %>% dplyr::select(-Description,-GO,-type)

  names(pheno1_input)[!names(pheno1_input) %in% c("ID", "Description", "GO","type")] <- paste0("pheno1_", names(pheno1_input)[!names(pheno1_input) %in% c("ID", "Description", "GO","type")])
  names(pheno2_input_simple)[!names(pheno2_input_simple) %in% c("ID", "Description", "GO","type")] <- paste0("pheno2_", names(pheno2_input_simple)[!names(pheno2_input_simple) %in% c("ID", "Description", "GO","type")])

  pheno1_pheno2_merged <- merge(pheno1_input,pheno2_input_simple,by.x="ID",by.y="ID")

  #pheno1_pheno2_final
  return(pheno1_pheno2_merged)
}

pheno1_pheno2_merged <- merge_pheno1_pheno2(GO_phenotype1,GO_phenotype2)


###########################
# EXTRACT CASE  NUMBER ####
###########################

print("Step 2/5 is processing")

extract_carriers <- function(data,
                             mutation
){
  #data<- pheno1_pheno2_merged
  #mutation<-"dn.mis2.ptv"

  mutation_pheno1 <- paste0(mutation, ".",data$pheno1_phenotype[1], sep="")
  mutation_pheno2 <- paste0(mutation, ".",data$pheno2_phenotype[1], sep="")

  data_NA<- data %>% filter (is.na(pheno2_geneID)=="TRUE" & is.na(pheno1_geneID)=="TRUE" )
  data_NA$pheno1_carrier_counts<- as.numeric("0")
  data_NA$pheno2_carrier_counts<- as.numeric("0")

  data_red <- data %>% filter (!ID %in% data_NA$ID)

  for (i in 1:nrow(data_red)){
    # df <- data.frame(geneID = data$geneID[1]) %>%
    #   separate_rows(geneID, sep = "/")
    # i<-as.numeric("2")

    ID_i <- data_red$ID[i]
    #ID_i <- "GO:0040029"
    gene_array_pheno1 <- unlist(strsplit(subset(data_red,ID==ID_i)$pheno1_geneID, "/"))
    gene_array_pheno2 <- unlist(strsplit(subset(data_red,ID==ID_i)$pheno2_geneID, "/"))

    pathway_genes <- unique(c(gene_array_pheno1,gene_array_pheno2))
    pathway_genes <- pathway_genes[!is.na(pathway_genes) & pathway_genes != ""]

    sum_all_pheno1 <- as.numeric(sum(subset(genes,entrez_id %in% pathway_genes)[[mutation_pheno1]]))
    sum_all_pheno2 <- as.numeric(sum(subset(genes,entrez_id %in% pathway_genes)[[mutation_pheno2]]))
    data_red$pheno1_carrier_counts[i] <- as.numeric(sum_all_pheno1)
    data_red$pheno2_carrier_counts[i] <- as.numeric(sum_all_pheno2)
  }

  # Combine all rows into one data.frame
  results_df <- rbind(data_red,data_NA)
  return(results_df)
}


pheno1_pheno2_merged_carriers <-extract_carriers(pheno1_pheno2_merged,
                                                             mutation=mutation)



######################################
#COMPARE phenotype1 vs phenotype2 ####
######################################
print("Step 3/5 is processing")

stat_compare_carriers <- function(data,
                                  phenotype1,
                                  phenotype2,
                                  mutation
) {
  # data <- pheno1_pheno2_merged_syn_carriers
  # phenotype1 <- phenotype1
  # phenotype2 <- phenotype2
  # mutation <- "dn.syn"

  #Description<-"SWI/SNF complex"

  # n_cores <- 4
  data_red <- data %>% filter(mutation == {{mutation}})
  categories <- unique(data_red$Description)

  run_test <- function(category) {
    data_i <- data_red %>% filter(Description == category)
    #data_i <- data_red %>% filter(Description == "behavioral fear response")

    pheno1_carrier_counts <- as.numeric(data_i$pheno1_carrier_counts)
    total_pheno1 <- as.numeric(subset(counts_df,hpo_suffix==phenotype1)$num_samples)
    pheno2_carrier_counts <- as.numeric(data_i$pheno2_carrier_counts)
    total_pheno2 <- as.numeric(subset(counts_df,hpo_suffix==phenotype2)$num_samples)

    matrix_data_i <- matrix(c(pheno1_carrier_counts,
                              total_pheno1 - pheno1_carrier_counts,
                              pheno2_carrier_counts,
                              total_pheno2 - pheno2_carrier_counts),
                            nrow = 2,
                            byrow = TRUE,
                            dimnames = list(
                              c("pheno1", "pheno2"),
                              c("mutated", "not_mutated")
                            ))
    fisher_test_result <- fisher.test(matrix_data_i)
    data.frame(
      pheno1_non_carriers = total_pheno1 - pheno1_carrier_counts,
      pheno2_non_carriers=total_pheno2 -pheno2_carrier_counts,
      OR_comparison_pheno1_pheno2 = fisher_test_result$estimate[[1]],
      CI_low_comparison_pheno1_pheno2 = fisher_test_result$conf.int[[1]],
      CI_up_comparison_pheno1_pheno2 = fisher_test_result$conf.int[[2]],
      p_value_comparison_pheno1_pheno2= fisher_test_result$p.value,
      Description = data_i$Description[1]

    )
  }

  result_list <- lapply(categories, run_test)

  result_df <- do.call(rbind, result_list[!sapply(result_list, is.null)])

  result_df$p_value_comparison_pheno1_pheno2_corrected <- p.adjust(
    result_df$p_value_comparison_pheno1_pheno2,
    method = "bonferroni",
    n = nrow(result_df)
  )
  result_df_merged <- merge(data,result_df,by.x="Description",by.y="Description")
  return(result_df_merged)
}

pheno1_pheno2_merged_carriers_stat <- stat_compare_carriers(pheno1_pheno2_merged_carriers,
                                                                        phenotype1=phenotype1,
                                                                        phenotype2= phenotype2,
                                                                        mutation=mutation )


#################
# WRITE ALL #####
#################
print("Step 4/5 is processing")

write.table(pheno1_pheno2_merged_carriers_stat,
             file = paste0("GO_pathway_", phenotype1, "_", phenotype2,"_", mutation, ".tsv"), 
             sep = "\t",
             row.names = FALSE,
             col.names = TRUE,
             quote=FALSE)

#######################
# MAKE A PLOT ALL #####
#######################
print("Step 5/5 is processing")

make_volano_plot <- function(data_input,
                                 mutation,
                                eigenvalue,
                                 title)
{
  # data_input <- GO_ALL_combined
  # mutation<-"dn.syn"
  # mutation<-"dn.mis2.ptv"
  #title<-"den"


  plot_input<- data_input %>% filter(pheno1_mutation=={{mutation}})
  phenotype1 <- data_input$pheno1_phenotype[1]
  phenotype2 <- data_input$pheno2_phenotype[1]

  eigenvalue <- as.numeric(eigenvalue)

  threshold <- as.numeric(0.05/eigenvalue)

  #use direct p-values
  plot_input$negLog10_p_fisher <- -log10(plot_input$p_value_comparison_pheno1_pheno2)
  plot_input$Log2_OR_Fisher <- log2(plot_input$OR_comparison_pheno1_pheno2)

  plot_input <- plot_input %>% mutate(
    significance_comparison_pheno1_pheno2= ifelse(p_value_comparison_pheno1_pheno2 <0.05,"significant","not_significant")
  )
  plot_input <- plot_input %>%
    mutate(
      significance_direction_pheno1_pheno2 = ifelse(
        OR_comparison_pheno1_pheno2 > 1, "phenotype1",
        ifelse(
          OR_comparison_pheno1_pheno2 == 0, "none", "phenotype2"
        )
      )
    )

  #colour Bonferroni corrected p-values
  plot_input$p_fisher_comparison_corrected_significance <- plot_input$p_value_comparison_pheno1_pheno2 < threshold


  plot_input <- plot_input %>%
    mutate(significance_label = case_when(
      p_fisher_comparison_corrected_significance & significance_comparison_pheno1_pheno2 == "significant" ~ "Significant (pheno1 vs pheno2 after correction)",
      p_fisher_comparison_corrected_significance == FALSE & significance_comparison_pheno1_pheno2 == "significant" ~ "Nominally significant (pheno1 vs pheno2)",
      TRUE ~ "Not significant"
    ))

  top5_input_1 <- plot_input %>%
    filter ( !is.na(Log2_OR_Fisher),
             !is.na(negLog10_p_fisher),
             is.finite(Log2_OR_Fisher),
             is.finite(negLog10_p_fisher), OR_comparison_pheno1_pheno2>1,p_value_comparison_pheno1_pheno2< threshold
    ) %>%  arrange(p_value_comparison_pheno1_pheno2) %>% slice_head(n = 15)

  top5_input_2 <- plot_input %>%
    filter ( !is.na(Log2_OR_Fisher),
             !is.na(negLog10_p_fisher),
             is.finite(Log2_OR_Fisher),
             is.finite(negLog10_p_fisher), OR_comparison_pheno1_pheno2<1,p_value_comparison_pheno1_pheno2<threshold
    ) %>%  arrange(p_value_comparison_pheno1_pheno2) %>% slice_head(n = 15)
  top5_input <- rbind(top5_input_1,top5_input_2)


  volcano_plot_pheno1_vs_pheno2 <- ggplot(subset(plot_input),
                                          aes(x = Log2_OR_Fisher, y = negLog10_p_fisher,
                                          )) +
    scale_fill_manual(
      values = c(
        "Significant (pheno1 vs pheno2 after correction)" = "darkred",
        "Nominally significant (pheno1 vs pheno2)" = "grey30",
        "Not significant" = "gray70"
      )
    ) +
    geom_jitter(
      aes(shape = mutation,
          #color = stroke_group,
          fill = significance_label
      ),
      alpha = 0.6,
      size = 3,
      width = 0.15,
      height = 0
    ) +
    scale_shape_manual(values = c(
      "dn.syn" = 24,       # circle
      "dn.mis2.ptv" = 21        # triangle
    )) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(threshold), linetype = "dashed") +
    theme_minimal() +
    labs(
      x = paste(phenotype2,"<---    ","log2(Odds Ratio)","    --->",phenotype1,sep=" "),
      y = "-log10(p-value)",
      title=title,
      shape = "Mutation",
      fill= "Significance Label"
    )  +
    guides(
      fill = guide_legend(
        override.aes = list(
          shape = 21,
          colour = "black"
        )
      )
    ) +
    #coord_cartesian(ylim = c(0, 8)) +
    theme(strip.text = element_text(face = "bold", color ="black", size = 12),
          legend.position = "bottom") +
    geom_text_repel(data = top5_input, aes(label = Description),
                    box.padding = 0.3, point.padding = 0.2, max.overlaps = Inf,
                    color = "black", size = 3)

  return(volcano_plot_pheno1_vs_pheno2)
}

plot1<- make_volano_plot(data_input=pheno1_pheno2_merged_carriers_stat,
                             mutation=mutation,
                             eigenvalue=eigenvalue,
                             title=paste(subset(counts_df,hpo_suffix == phenotype2)$hpo_marker_name,"vs",
                                         subset(counts_df,hpo_suffix == phenotype1)$hpo_marker_name,"\n de novo mis2 or PTV"))

pdf(paste("Volcano_",phenotype1,"_",phenotype2,"_",mutation,".pdf", sep = ""), width = 10, height = 7)
print(plot1)
dev.off()