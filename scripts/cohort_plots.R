#!/usr/bin/env Rscript
#Author: Alba Sanchis-Juan
##This script makes general QC plots for a specific cohort given a PED file

#############
##Libraries##
#############

library(ggplot2)
library(ontologyIndex)
library(RSQLite)
library(reshape)
library(data.table)
library(dplyr)
library(stringr)
data(hpo)

hpo_db <- fread("/Users/an436/Documents/work/resources/hpo/hpo_names.txt")
top_nodes <- setdiff(hpo$children[["HP:0000118"]], c("HP:0045027","HP:0001608","HP:0002664","HP:0000152","HP:0000769","HP:0001197","HP:0040064"))

#############
##Load data##
#############

manifest <- fread("~/Documents/work/projects/CMG/metadata/Talkowski_Metadata_Export.ped.txt", stringsAsFactors = T)
manifest$hpo_present_fix <- sapply(manifest$hpo_present, function(x) 
  paste(str_extract_all(string = x, pattern = "HP:[0-9]+")[[1]], collapse = ",")
)
all_samples <- subset(manifest, affected == "Affected" & data_type == "WGS")$subject_id
hpos <- subset(manifest, affected == "Affected" & data_type == "WGS")$hpo_present_fix

ped <- fread("~/Documents/work/projects/CMG/phase2/pedigrees/CMG_phase1_2-20210622_excl_repeat.ped",
             col.names = c("family_id", "subject_id", "father", "mother", "sex", "affected"))

#########
##PLOTS##
#########

#Number of HPO/sample

hpo_counts <- sapply(hpos, function(s) length(strsplit(as.character(s), ",")[[1]]))
# boxplot(hpo_counts)
hpo_counts <- data.frame("n_hpos" = hpo_counts)
p_hpo <- ggplot(hpo_counts, aes(x = "", y = n_hpos)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    # axis.text.x = element_text(angle = 90, hjust = 1),
    # axis.text.x = element_blank(),
    # legend.position = "none",
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 16),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black")) +
  labs(x = "Sample", y="#HPO")
p_hpo

#HPOS main distribution

hpo_names <- data.frame("hpos" = unlist(sapply(hpos, function(s) strsplit(as.character(s), ",")[[1]])))

hpo_db_sub <- unique(subset(hpo_db, select = c(HPO_ID, HPO_DESCRIPTION)))
hpo_names_ids <- merge(hpo_names, hpo_db_sub, all.x=T, all.y=F, by.x="hpos", by.y="HPO_ID")
hpo_names_count <- plyr::count(hpo_names_ids$HPO_DESCRIPTION)
hpo_names_count <- hpo_names_count[order(hpo_names_count$freq),]
hpo_names_count$x <- factor(hpo_names_count$x, levels = hpo_names_count$x)

# p_phen <- ggplot(subset(hpo_names_count, freq >200 & is.na(x) == F), aes(x = x, y = freq)) +
p_phen <- ggplot(subset(hpo_names_count, freq >20 & is.na(x) == F), aes(x = x, y = freq)) +
  geom_col() +
  labs(x = "Phenotype", y = "Count") +
  coord_flip() +
  theme_bw() +
  theme(
    # axis.text.x = element_text(angle = 90, hjust = 1),
    # axis.text.x = element_blank(),
    # legend.position = "none",
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 16),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"))

p_phen

##Heatmap

wgs <- subset(manifest, affected == "Affected" & data_type == "WGS")$subject_id

wgs_count <- do.call(rbind, lapply(wgs, function(s){
  #s <- samples[1]
  indv_hpos <- strsplit(as.character(manifest[manifest$subject_id == s,]$hpo_present_fix), ",")[[1]]
  indv_terms_with_ancs <- lapply(indv_hpos, get_ancestors, ontology=hpo)
  indv_terms_with_ancs_unlist <- unique(unlist(lapply(indv_terms_with_ancs, function(x) x[x %in% top_nodes])))
  
  if(length(indv_terms_with_ancs_unlist) >0){
    return(data.frame(sample = s, 
                      hpo = indv_terms_with_ancs_unlist))
  }
}))

wgs_count$hpo_name <- hpo$name[as.character(wgs_count$hpo)]

hpo_categories <- unique(wgs_count$hpo_name)

wgs_categories <- do.call(rbind, lapply(hpo_categories, function(c){
  # c<- hpo_categories[1]
  category_samples <- unique(subset(wgs_count, hpo_name == c)$sample)
  df <- data.frame(table(subset(wgs_count, sample %in% category_samples)$hpo_name))
  if(length(hpo_categories[!hpo_categories %in% df$Var1]) > 0){
    df <- rbind(df, data.frame("Var1" = hpo_categories[!hpo_categories %in% df$Var1], "Freq" = 0))
  }
  df$y_axis <- c
  return(df)
})
)

wgs_categories$both_categories <- apply(wgs_categories, 1, function(r) {
  paste(sort(c(as.character(r['Var1']), as.character(r['y_axis']))), collapse = ";")
})

wgs_categories2 <- unique(wgs_categories[,c(2,4)])
wgs_categories2$x_sort <- sapply(wgs_categories2$both_categories, function(x) strsplit(x, ";")[[1]][1])
wgs_categories2$y_sort <- sapply(wgs_categories2$both_categories, function(x) strsplit(x, ";")[[1]][2])

wgs_categories2$Label <- wgs_categories2$Freq
font_threshold <- 200
wgs_categories2$font_col <- ifelse(wgs_categories2$Freq < font_threshold, 0, 1)
wgs_categories2$font_col <- factor(wgs_categories2$font_col, levels = unique(wgs_categories2$font_col))

p_wgs_heat <- ggplot(wgs_categories2, aes(x_sort, y_sort)) + 
  geom_tile(aes(fill = Freq)) + 
  # scale_fill_gradient(low = "blue", high = "yellow", name = "Number of HPOterms per category"
  #                     # , limits=c(0,50)
  # ) +
  scale_fill_gradientn(colours = c("#FFFAC2", "blue", "black"), name = "Number of cases",
                       # values = scales::rescale(c(-0.5, -0.05, 0, 0.05, 0.5))
  ) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size=10)
        , axis.text.y  = element_text(size=10)
        , axis.title = element_text(size = 16)
        , legend.position = c(0.8, 0.2)
        , legend.box.just = "right"
        # , legend.justification = c(0.8,0.2)
        # , legend.text = element_text(margin = margin(0, 0, 0, 0, "pt"), size = 10)
        , panel.border = element_blank() 
        , panel.grid.major = element_blank()
        , panel.grid.minor = element_blank()
        , axis.line = element_line(colour = "black")) +
  # geom_text(aes(label = Freq)) +
  geom_text(aes(label = Label, colour=font_col), show.legend = FALSE) +
  scale_colour_manual(values=c("black", "white")) +
  labs(x="HPO category", y="HPO category", title = "MGRC WGS - Phase 2")

p_wgs_heat

# Sex distribution

ped[ped$affected %in% c("Unknown", 0),]$affected <- -9
ped[ped$sex == "pending",]$sex <- "-9"
ped[ped$affected == 0,]$affected <- "-9"

ped$affected <- factor(ped$affected, levels = c(-9, 1, 2))
ped$sex <- factor(ped$sex, levels = c(-9, 1, 2))


df <- data.frame(reshape2::melt(table(ped$sex)), "variable" = "sex")
df2 <- data.frame(reshape2::melt(table(ped$affected)), "variable" = "affected")

df$Var1 <- factor(df$Var1, levels = c("1", "2", "-9"))
df2$Var1 <- factor(df2$Var1, levels = c("1", "2", "-9"))

df_p <- rbind(df, df2)
df_p$Var1 <- as.factor(df_p$Var1)
p1 <- ggplot(df_p, aes(x = variable, y = value, group = Var1)) +
  geom_col(aes(fill = Var1)) +
  # coord_flip() +
  scale_fill_manual(values = c("#F0C808", "#086788", "#B56576")) +
  theme_bw() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),
    # axis.text.x = element_blank(),
    # legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"))
p1

# Family size counts

ped_counts <- plyr::count(ped$family_id)
ped_counts$freq <- factor(ped_counts$freq, levels = 1:9)

p2 <- ggplot(ped_counts, aes(x = freq)) +
  geom_bar() +
  labs(x = "Family size", y = "Count") +
  theme_bw() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1),
    # axis.text.x = element_blank(),
    # legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    panel.border = element_blank(), 
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"))

p2


##PLOTS
ggsave(plot = p_hpo, filename = "~/Documents/work/projects/CMG/phase2/plots/WGS_num_hpos.png"
       , width = 4, height = 6)

ggsave(plot = p_phen, filename = "~/Documents/work/projects/CMG/phase2/plots/WGS_phenotypes_count.png"
       , width = 8, height = 6)

ggsave(plot = p_wgs_heat, filename = "~/Documents/work/projects/CMG/phase2/plots/WGS_phenotypes.png"
       , width = 12, height = 8)

ggsave(plot = p1, filename = "~/Documents/work/projects/CMG/phase2/plots/WGS_sex_affected.png"
       , width = 4, height = 6)

ggsave(plot = p2, filename = "~/Documents/work/projects/CMG/phase2/plots/WGS_family_size.png"
       , width = 8, height = 6)

