#!/usr/bin/env Rscript

##Author: Alba Sanchis-Juan
##Description: Post-hoc QC analysis of filtered de novo SVs


##Load libraries
library(tidyr)
library(ggplot2)
library(UpSetR)
library(data.table)
library(dplyr)
library(gridExtra)
library(grid)
library(tidyverse)
library("grid")
library("ggplotify")

args = commandArgs(trailingOnly=TRUE)

input_bed <- args[1]
input_outliers <-args [2]
input_ped <- args[3]
out_file <- args[4]
#known_out_file <- args[5]


ped <- fread(input_ped)


##De novo CMG calls
date <- "20220427"
denovo <- as.data.frame(read.table(input_bed,header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))
outliers <- as.data.frame(read.table(input_outliers,header = TRUE, sep="\t",stringsAsFactors=FALSE, quote=""))


denovo %>%  na.omit()
denovo<- subset(denovo, chrom!= 'chrom')

denovo$SVTYPE <- factor(denovo$SVTYPE, levels = rev(c("DEL", "DUP", "INS", "INV", "CPX", "CTX")))

denovo_ins <- subset(denovo, SVTYPE == "INS")
denovo_del <- subset(denovo, SVTYPE == "DEL")
denovo_dup <- subset(denovo, SVTYPE == "DUP")
denovo_inv <- subset(denovo, SVTYPE == "INV")
denovo_cpx <- subset(denovo, SVTYPE == "CPX")
denovo_ctx <- subset(denovo, SVTYPE == "CTX")


##Get some numbers
length(unique(denovo$name)) #of unique denovos
length(unique(denovo$sample)) #number of samples with denovos

length(unique(outliers$name))
length(unique(outliers$sample))
outliers_in_gd <- subset(outliers, select= c(name, in_gd))



# denovo2 <-subset(denovo, !name %in% subset(denovo, GQ == 1 & ALGORITHMS == "wham")$name)
#known_denovo <- c("phase2_DEL_chr15_52", "phase2_INS_chr16_366", "phase2_DEL_chr16_3125", "phase2_DEL_chr16_4510",
                  #"phase2_DEL_chr19_484", "phase2_DEL_chr2_2795", "phase2_DEL_chr5_2228", "phase2_DUP_chr1_1428", "phase2_DEL_chr1_3303",
                  #"phase2_DEL_chr11_64", "phase2_INV_chr11_4", "phase2_DUP_chr12_4348", "phase2_DEL_chr16_590", 
                  #"phase2_DEL_chr16_3432", "phase2_DEL_chr17_1588", "phase2_DEL_chr19_2844", "phase2_DEL_chr2_8361", "phase2_DEL_chr3_593",
                  #"phase2_CTX_chr5_1", "phase2_DEL_chr9_1208", "phase2_DUP_chr9_4575", "phase2_DUP_chrX_1004")
#phase2_DEL_chr22_384 not denovo, not trio

#known_denovos <- subset(denovo, name == "phase2_INS_chr16_366" | name == "phase2_DEL_chr16_3125" | name == "phase2_DEL_chr16_4510" | name == "phase2_DEL_chr16_590" | name == "phase2_DEL_chr16_3432", select= c(name, in_gd, AF))

#known_denovos_in_outliers <- subset(outliers, name == "phase2_INS_chr16_366" | name == "phase2_DEL_chr16_3125" | name == "phase2_DEL_chr16_4510" | name == "phase2_DEL_chr16_590" | name == "phase2_DEL_chr16_3432", select= c(name, in_gd, AF))


#table(known_denovo %in% denovo$name)
#table(known_denovo %in% outliers$name)
#known_denovo[!known_denovo %in% denovo$name] #gives the known denovos that we do not find in this data

library(dplyr)

denovo %>%
  # subset(chrom != "chrX") %>%
  select(sample, name) %>%
  group_by(sample) %>%
  tally() -> sample_count

sample_count %>%
  ggplot(aes(y = n)) +
  geom_boxplot() -> boxplot

#exclude tmporarily chrX
# denovo <- subset(denovo, chrom != "chrX")

##Define colors
del_col <- "#D43925"
dup_col <- "#2376B2"
ins_col <- "#D474E0"
inv_col <- "#FA9627"
ctx_col <- "#638E6C"
cpx_col <- "#4DA1A9"

##Make plots
denovo$chrom <- factor(denovo$chrom, levels = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"))


##Count by chromosome
denovo %>%
  select(chrom, name, SVTYPE) %>%
  unique() %>%
  ggplot(aes(x = chrom, fill = SVTYPE)) +
  geom_bar() +
  labs(x="Chromosome", y="Number de novo SVs", title= "Number of De Novo SVs per Chromosome") +
  scale_fill_manual(values = rev(c(del_col, dup_col, ins_col, inv_col, cpx_col, ctx_col))) + 
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size=25),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")) -> p_chr_count
p_chr_count

# #Count by sample
# denovo %>%
#   select(sample, name, SVTYPE) %>%
#   unique() %>%
#   ggplot(aes(x = sample, fill = SVTYPE)) +
#   geom_bar() +
#   labs(x="Sample", y="Number de novo SVs", title= "Number of De Novo SVs per Sample") +
#   scale_fill_manual(values = rev(c(del_col, dup_col, ins_col, inv_col, cpx_col, ctx_col))) + 
#   theme_classic() +
#   theme(
#     legend.position = "right",
#     axis.text = element_text(size = 10),
#     axis.title = element_text(size = 11),
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.border = element_blank(),
#     axis.line = element_line(colour = "black")) -> p_sample_count
# p_sample_count

#count by sample

df <- as.data.frame(table(denovo$sample))
df$samples <- as.character(df$Var1)
df$num_of_denovos <- as.numeric(as.character(df$Freq))
df$samples <- as.factor(df$samples)


p <- ggplot(df, aes(y=num_of_denovos)) + geom_boxplot() + labs(title = "Number of de novo SVs per sample", y = "Number of de novos", x = "Samples") + theme_classic() + theme(
  axis.text = element_text(size = 20),
  axis.title = element_text(size = 20),
  plot.title = element_text(size=25))

##De novo count by allele frequency
denovo$AF <- as.numeric(denovo$AF)  
denovo$AC <- as.numeric(denovo$AC) 
ac_1_freq <- min(subset(denovo, AC == 1)$AF)
max_AF_str <- toString(round(max(denovo$AF), digits=3))
# denovo$AF_interv <- cut(denovo$AF, breaks = c(0, ac_1_freq, 0.001, 0.01, 0.03, max(denovo$AF)),
#                         labels = c("AC=1", "AC 1 - AF<=0.001", "AF 0.001-0.01", "AF 0.01-0.03", "AF>0.03"))
denovo$AF_interv <- cut(denovo$AF, breaks = c(0, ac_1_freq, 0.001, max(denovo$AF)),
                        labels = c("AC=1", paste("AC 1 - AF<=", max_AF_str), paste("AF >", max_AF_str)))

denovo %>%
  select(AF_interv, name, SVTYPE) %>%
  unique() %>%
  ggplot(aes(x = AF_interv, fill = SVTYPE)) +
  geom_bar() +
  # scale_y_log10()+
  labs(x="Allele Frequency", y="Number de novo SVs", title="Number of De Novo\n SVs by Frequency") +
  scale_fill_manual(values = rev(c(del_col, dup_col, ins_col, inv_col, cpx_col, ctx_col))) + 
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size=25),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")) -> p_af_count
p_af_count


denovo_in_gd <- subset(denovo, in_gd == "True")
denovo_not_in_gd <- subset(denovo, in_gd == "False")


subset(denovo_in_gd, AF_interv == "AF 0.001-0.01")$chrom_type_sample


denovo_in_gd_sub <- subset(denovo_in_gd, select= c(name,AF_interv))

#denovo_not_in_gd_sub <- subset(denovo_not_in_gd, select= c(name,AF,gnomAD_V2_AF))


denovo_in_gd %>%
  select(AF_interv, name, SVTYPE) %>%
  unique() %>%
  ggplot(aes(x = AF_interv, fill = SVTYPE)) +
  geom_bar() +
  # scale_y_log10()+
  labs(x="Allele Frequency", y="Number de novo SVs", title= "Number of De Novo SVs \nin Genomic Disorder Regions by Frequency") +
  scale_fill_manual(values = rev(c(del_col, dup_col, ins_col, inv_col, cpx_col, ctx_col))) + 
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size=25),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")) -> p_af_count_in_gd
p_af_count_in_gd

denovo_not_in_gd %>%
  select(AF_interv, name, SVTYPE) %>%
  unique() %>%
  ggplot(aes(x = AF_interv, fill = SVTYPE)) +
  geom_bar() +
  # scale_y_log10()+
  labs(x="Allele Frequency", y="Number de novo SVs", title= "Number of De Novo SVs not in \n Genomic Disorder Regions by Frequency") +
  scale_fill_manual(values = rev(c(del_col, dup_col, ins_col, inv_col, cpx_col, ctx_col))) + 
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size=25),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")) -> p_af_count_not_in_gd
p_af_count_not_in_gd

# 
# lay <- rbind(c(1,1,2,2,3,3))
# ml <- grid.arrange(p_af_count, p_af_count_in_gd, p_af_count_not_in_gd, layout_matrix = lay, top=textGrob("Chr16 AF Count CMG DeNovo SV Data"))
# ggsave("chr16.p_af_count.pdf", ml, width = 9, height = 5)
# 



##De novo count by size
denovo$SVLEN <- as.numeric(denovo$SVLEN)  
#denovo[denovo$SVLEN < 1,]$SVLEN <- 1
denovo$SVLEN_interv <- cut(denovo$SVLEN, breaks = c(0, 100, 500, 2500, 10000, 50000, max(denovo$SVLEN)),
                           labels = c("<100bp", "101-500bp", "501bp-2.5kb", "2.5kb-10kb", "10kb-50kb", ">50kb"))

denovo %>%
  select(SVLEN_interv, name, SVTYPE) %>%
  unique() %>%
  ggplot(aes(x = SVLEN_interv, fill = SVTYPE)) +
  geom_bar() +
  # scale_y_log10()+
  labs(x="Size of SV", y="Number de novo SVs", title= "De Novo SVs by Size") +
  scale_fill_manual(values = rev(c(del_col, dup_col, ins_col, inv_col, cpx_col, ctx_col))) + 
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size=25),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")) -> p_size_count
p_size_count

##Evidence plot

SR_denovo <- subset(denovo, EVIDENCE_FIX == "SR", select = c(ALGORITHMS, GQ))
caller_quality_denovos <- subset(denovo, select = c(name,ALGORITHMS, GQ, EVIDENCE))


denovo %>%
  select(name, SVTYPE, EVIDENCE_FIX) %>%
  unique() %>%
  ggplot(aes(x = forcats::fct_infreq(EVIDENCE_FIX), fill = SVTYPE)) +
  geom_bar() +
  # scale_y_log10()+
  labs(x="Evidence", y="Number de novo SVs", title= "De Novo SVs by Evidence") +
  scale_fill_manual(values = rev(c(del_col, dup_col, ins_col, inv_col, cpx_col, ctx_col))) + 
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 20),
    plot.title = element_text(size=25),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black")) -> p_evidence
p_evidence

##Annotation plot
annot <- grep("_", names(denovo), value = T)
annot <- grep("LINCRNA", annot, value = T, invert = T)
# annot <- annot[!annot %in% c("PROTEIN_CODING__NEAREST_TSS", "PROTEIN_CODING__INTERGENIC")]

lof <- grep("LOF", annot, value = T)
cg <- grep("COPY_GAIN", annot, value = T)
msv <- grep("MSV_EXON_OVR", annot, value = T)
inv_span <- grep("INV_SPAN", annot, value = T)
dup_partial <- grep("DUP_PARTIAL", annot, value = T)
intronic <- grep("INTRONIC", annot, value = T)
utr <- grep("UTR", annot, value = T)
promoter<- grep("PROMOTER", annot, value = T)
nearest_tss<- grep("NEAREST_TSS", annot, value = T)
intergenic<- grep("INTERGENIC", annot, value = T)
coding <- grep("CODING", annot, value = T)
coding <- coding[!coding %in% c('PROTEIN_CODING__NEAREST_TSS', "PROTEIN_CODING__INTERGENIC")]

denovo$is_lof <- apply(denovo[,lof, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_cg <- apply(denovo[,cg, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_msv <- apply(denovo[,msv, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_inv_span <- apply(denovo[,inv_span, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_dup_partial <- apply(denovo[,dup_partial, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_intronic <- apply(denovo[,intronic, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_utr <- apply(denovo[,utr, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_promoter <- apply(denovo[,promoter, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_nearest_tss <- apply(denovo[,nearest_tss, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))
denovo$is_intergenic <- apply(denovo[,intergenic, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == '' | r == FALSE))
denovo$is_coding <- apply(denovo[,coding, drop = FALSE], 1, function(r) !all(is.na(r) | r == 'NA' | r == ''))

keep_cols <- c(grep("SVTYPE", names(denovo)), (ncol(denovo)-9-1):ncol(denovo))
denovo_cons <- denovo[,keep_cols]

ref_cols <- grep("is_", names(denovo_cons))
denovo_cons_type <- cbind(denovo_cons[,'SVTYPE', drop = FALSE], (denovo_cons[,ref_cols])*1)

# denovo_cons_coding <- subset(denovo_cons, is_coding == 1)
# denovo_cons_non_coding <- subset(denovo_cons, is_coding == 0)

sv_cols <- list(
  list(query = elements, 
       params = list("SVTYPE", c("DEL","DUP", "INS", "INV", "CPX", "CTX")), color = del_col, active = T
       , query.name = "DEL"
  ), #red
  list(query = elements, 
       params = list("SVTYPE", c("DUP", "INS", "INV", "CPX", "CTX")), color = dup_col, active = T
       , query.name = "DUP"
  ), #blue
  list(query = elements, 
       params = list("SVTYPE", c("INS", "INV", "CPX", "CTX")), color = ins_col, active = T
       , query.name = "INS"
  ), #orange
  list(query = elements, 
       params = list("SVTYPE", c("INV", "CPX", "CTX")), color = inv_col, active = T
       , query.name = "INV"
  ), #green
  list(query = elements, 
       params = list("SVTYPE", "CPX", "CTX"), color = cpx_col, active = T
       , query.name = "CPX" )
  # ), #purple
  #list(query = elements, 
  #params = list("SVTYPE", "CTX"), color = ctx_col, active = T,
  #query.name = "CTX"
  #)
)


denovo_cons_type %>% 
  upset(sets = grep("is_", names(denovo_cons_type), value = T),
        queries = sv_cols,
        order.by = "freq", 
        set_size.show = TRUE,
        # set_size.scale_max = 15000, 
        text.scale = 2.5, 
        query.legend = "top") -> p_denovo_upset_all

grob_annotation_upset_plot <- as.grob(p_denovo_upset_all)

# denovo_cons_coding %>% 
#   upset(sets = names(denovo_cons_coding),
#         order.by = "freq", 
#         set_size.show = TRUE,
#         # set_size.scale_max = 15000, 
#         text.scale = 2, 
#         query.legend = "top") -> denovo_upset_coding
# 
# denovo_cons_non_coding %>% 
#   upset(sets = names(denovo_cons_non_coding),
#         order.by = "freq", 
#         set_size.show = TRUE,
#         # set_size.scale_max = 15000, 
#         text.scale = 2, 
#         query.legend = "top") -> denovo_upset_non_coding
# denovo_upset_coding


##KNOWN DENOVO PLOTS
##De novo count by size

#subset(denovo, name %in% known_denovo) %>%
  #select(SVLEN_interv, name, SVTYPE) %>%
  #unique() %>%
  #ggplot(aes(x = SVLEN_interv, fill = SVTYPE)) +
  #geom_bar() +
  # scale_y_log10()+
  #labs(x="Allele Frequency", y="Number de novo SVs", title= "Known De Novo SVs by Size") +
  #scale_fill_manual(values = rev(c(del_col, dup_col, ins_col, inv_col, ctx_col))) + 
  #theme_classic() +
  #theme(
    #legend.position = "right",
    #axis.text = element_text(size = 10),
    #axis.title = element_text(size = 11),
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.border = element_blank(),
#     axis.line = element_line(colour = "black")) -> p_size_count_known
# p_size_count_known

# ac_1_freq <- min(subset(denovo, AC == 1)$AF)
# denovo$AF_interv <- cut(denovo$AF, breaks = c(0, ac_1_freq, 0.001, 0.01, 0.03, max(denovo$AF)),
#                         labels = c("AC=1", "AC 1 - AF<=0.001", "AF 0.001-0.01", "AF 0.01-0.03", "AF>0.03"))
# 
# subset(denovo, name %in% known_denovo) %>%
#   select(AF_interv, name, SVTYPE) %>%
#   unique() %>%
#   ggplot(aes(x = AF_interv, fill = SVTYPE)) +
#   geom_bar() +
#   # scale_y_log10()+
#   labs(x="Allele Frequency", y="Number de novo SVs", title= "Known De Novo SVs by Frequency") +
#   scale_fill_manual(values = rev(c(del_col, dup_col, ins_col, inv_col, ctx_col))) + 
#   theme_classic() +
#   theme(
#     legend.position = "right",
#     axis.text = element_text(size = 10),
#     axis.title = element_text(size = 11),
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.border = element_blank(),
#     axis.line = element_line(colour = "black")) -> p_af_count_known
# p_af_count_known
# 
# 
# subset(denovo, name %in% known_denovo) %>%
#   select(name, SVTYPE, EVIDENCE_FIX) %>%
#   unique() %>%
#   ggplot(aes(x = forcats::fct_infreq(EVIDENCE_FIX), fill = SVTYPE)) +
#   geom_bar() +
#   # scale_y_log10()+
#   labs(x="Evidence", y="Number de novo SVs", title = "Number of Known De Novo SVs by Evidence") +
#   scale_fill_manual(values = rev(c(del_col, dup_col, ins_col, inv_col, ctx_col))) +
#   theme_classic() +
#   theme(
#     legend.position = "right",
#     axis.text = element_text(size = 10),
#     axis.title = element_text(size = 11),
#     axis.text.x = element_text(angle = 45, hjust = 1),
#     panel.border = element_blank(),
#     axis.line = element_line(colour = "black")) -> p_evidence_known
# p_evidence_known


#Creating a panel of plots
lay <- rbind(c(1,1,2,2,2,2,2,2),
             c(3,3,NA,4,4,NA,5,5),
             c(6,6,6,6,7,7,7,7),
             c(8,8,8,8,8,8,NA,NA))

ml <- grid.arrange(p, p_chr_count, p_af_count, p_af_count_in_gd, p_af_count_not_in_gd, p_size_count, p_evidence, grob_annotation_upset_plot, layout_matrix = lay, top=textGrob("De Novo SV Data", gp=gpar(fontsize=30)))
ggsave(out_file, ml, width = 35, height = 45)


# 
# 
# lay <- rbind(c(1,1,2,2,3,3))
# ml <- grid.arrange(p_af_count_known, p_size_count_known, p_evidence_known, layout_matrix = lay, top=textGrob("Chr16 Known CMG De Novo SV Data"))
# ggsave(known_out_file, ml, width = 8, height = 5)

















