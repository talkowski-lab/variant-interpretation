library("dplyr")
library(data.table)

args = commandArgs(trailingOnly=TRUE)

sample_list <- args[1]
cram_list <- args[2]
ped_input <- args[3]
out <- args[4]

ped <- fread(ped_input)

ped <- subset(ped, select=c(family_id, subject_id, maternal_id, paternal_id, family_size))

#for this script I am using it will only work with trios
trios <- subset(ped, family_size==3)

sample_list <- as.data.frame(read.table(sample_list, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
cram_list <- as.data.frame(read.table(cram_list, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
#crai_list <- as.data.frame(read.table(crai_list, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
data_table <- data.frame(sample_list, cram_list)



colnames(data_table) <- c('subject_id','bam_or_cram_file')
#Xuefang said we can get it working for if only one parent also (this means there are two kids and one parent)
#this gets only pb in subject_id column
trios <- subset(trios, (paternal_id != 0) | (maternal_id != 0))

#puts parent in both maternal and paternal if only one parent --> BECAUSE THIS MEANS THERE ARE SIBLINGS,
#FAMILY_ID IS NOT UNIQUE
trios$maternal_id <- ifelse(trios$maternal_id == 0, trios$paternal_id, trios$maternal_id)
trios$paternal_id <- ifelse(trios$paternal_id == 0, trios$maternal_id, trios$paternal_id)
  
pb_table <- subset(data_table, data_table$subject_id %in% trios$subject_id)
colnames(pb_table) <- c('subject_id','pb_cram_file')
  

#get pb cram into ped
trios <- merge(trios,pb_table, by="subject_id")
  
  
#set up mo and fa tables to get crams in trio table
mo_table <- subset(data_table, data_table$subject_id %in% trios$maternal_id)
colnames(mo_table) <- c('maternal_id','mo_cram_file')
fa_table <- subset(data_table, data_table$subject_id %in% trios$paternal_id)
colnames(fa_table) <- c('paternal_id','fa_cram_file')
  
  
#add mo and fa crams to table
trios <- merge(trios,mo_table, by="maternal_id", all.x = TRUE)
trios <- merge(trios,fa_table, by="paternal_id", all.x = TRUE)

trios <- subset(trios, select=c("family_id", "subject_id", "maternal_id", "paternal_id", "pb_cram_file","mo_cram_file","fa_cram_file"))
colnames(trios)[1] <- "entity:family_id"

write.table(trios, out, sep="\t", quote=F, row.names=F)

