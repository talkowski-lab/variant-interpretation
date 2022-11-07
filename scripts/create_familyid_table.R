library("dplyr")
library(data.table)

args = commandArgs(trailingOnly=TRUE)

sample_list <- args[1]
cram_list <- args[2]
crai_list <- args[3]
ped_input <- args[4]
out_trios <- args[5]
# out_trios_single_parent <- args[5]
out_singletons <- args[6]

ped <- fread(ped_input)

ped <- subset(ped, select=c(family_id, subject_id, maternal_id, paternal_id, family_size))

#Trios
trio_samples <- subset(ped, family_size==3)

sample_list <- as.data.frame(read.table(sample_list, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
cram_list <- as.data.frame(read.table(cram_list, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
crai_list <- as.data.frame(read.table(crai_list, header = FALSE, sep="\t",stringsAsFactors=FALSE, quote=""))
data_table <- data.frame(sample_list, cram_list, crai_list)



colnames(data_table) <- c('subject_id','bam_or_cram_file')
#Xuefang said we can get it working for if only one parent also (this means there are two kids and one parent)
#this gets only pb in subject_id column
trios <- subset(trio_samples, (paternal_id != 0) & (maternal_id != 0))

pb_table <- subset(data_table, data_table$subject_id %in% trios$subject_id)
colnames(pb_table) <- c('subject_id','pb_cram_file', 'pb_crai_file')


#get pb cram into ped
trios <- merge(trios,pb_table, by="subject_id")


#set up mo and fa tables to get crams in trio table
mo_table <- subset(data_table, data_table$subject_id %in% trios$maternal_id)
colnames(mo_table) <- c('maternal_id','mo_cram_file', 'mo_crai_file')
fa_table <- subset(data_table, data_table$subject_id %in% trios$paternal_id)
colnames(fa_table) <- c('paternal_id','fa_cram_file', 'fa_crai_file')


#add mo and fa crams to table
trios <- merge(trios,mo_table, by="maternal_id", all.x = TRUE)
trios <- merge(trios,fa_table, by="paternal_id", all.x = TRUE)

trios <- subset(trios, select=c("family_id", "subject_id", "maternal_id", "paternal_id", "pb_cram_file","mo_cram_file","fa_cram_file", "pb_crai_file","mo_crai_file","fa_crai_file"))
colnames(trios)[1] <- "entity:trio_id"

write.table(trios, out_trios, sep="\t", quote=F, row.names=F)


# #Trios with single parent, duos, quads
# trio_samples <- subset(ped, family_size==3)
# 
# 
# 
# #FAMILY_ID IS NOT UNIQUE AFTER THIS STEP
# trios_single_parent <- subset(trio_samples, ((paternal_id != 0) & (maternal_id ==0)) | ((maternal_id != 0) & (paternal_id == 0)))
# 
# #puts parent in both maternal and paternal if only one parent --> BECAUSE THIS MEANS THERE ARE SIBLINGS,
# #FAMILY_ID IS NOT UNIQUE
# trios_single_parent$maternal_id <- ifelse(trios_single_parent$maternal_id == 0, trios_single_parent$paternal_id, trios_single_parent$maternal_id)
# trios_single_parent$paternal_id <- ifelse(trios_single_parent$paternal_id == 0, trios_single_parent$maternal_id, trios_single_parent$paternal_id)
#   
# pb_table <- subset(data_table, data_table$subject_id %in% trios_single_parent$subject_id)
# colnames(pb_table) <- c('subject_id','pb_cram_file')
#   
# 
# #get pb cram into ped
# trios_single_parent <- merge(trios_single_parent,pb_table, by="subject_id")
#   
#   
# #set up mo and fa tables to get crams in trio table
# mo_table <- subset(data_table, data_table$subject_id %in% trios_single_parent$maternal_id)
# colnames(mo_table) <- c('maternal_id','mo_cram_file')
# fa_table <- subset(data_table, data_table$subject_id %in% trios_single_parent$paternal_id)
# colnames(fa_table) <- c('paternal_id','fa_cram_file')
#   
#   
# #add mo and fa crams to table
# trios_single_parent <- merge(trios_single_parent,mo_table, by="maternal_id", all.x = TRUE)
# trios_single_parent <- merge(trios_single_parent,fa_table, by="paternal_id", all.x = TRUE)
# 
# trios_single_parent <- subset(trios_single_parent, select=c("family_id", "subject_id", "maternal_id", "paternal_id", "pb_cram_file","mo_cram_file","fa_cram_file"))
# colnames(trios_single_parent)[1] <- "entity:single_parent_trios_id"
# 
# write.table(trios_single_parent, out_trios_single_parent, sep="\t", quote=F, row.names=F)
# 

#Singletons
singletons <- subset(ped, family_size==1)
pb_table <- subset(data_table, data_table$subject_id %in% singletons$subject_id)
colnames(pb_table) <- c('subject_id','pb_cram_file', 'pb_crai_file')


#get pb cram into ped
singletons <- merge(singletons,pb_table, by="subject_id")

singletons <- subset(singletons, select=c("family_id", "subject_id", "maternal_id", "paternal_id", "pb_cram_file", 'pb_crai_file'))
colnames(singletons)[1] <- "entity:singleton_id"
write.table(singletons, out_singletons, sep="\t", quote=F, row.names=F)



