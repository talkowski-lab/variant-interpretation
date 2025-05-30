#!/usr/bin/env Rscript
#Author: Alba Sanchis-Juan

#####################
##Define parameters##
#####################
library("optparse")

option_list = list(
  make_option(c("-c", "--rconfig"), type="character", default=NULL,
              help="Config file", metavar="character"),
  make_option(c("-o", "--out_file"), type="character", default=NULL,
              help="Path to output file", metavar="character"),
  make_option(c("-u", "--rfunctions"), type="character", default=NULL,
              help="Path to r functions file", metavar="character"),
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
        help="Verbosity [default]")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

rconfig <- opt$rconfig
#variants_path <- opt$input_file
#variants_gt_path <- opt$input_gt
#fam <- opt$fam
#manifest_path <- opt$manifest
#gd_path <- opt$gd_path
out_file <- opt$out_file
rfunctions <- opt$rfunctions

if (is.null(rconfig)){
  print_help(opt_parser)
  stop("The config file must be supplied.n", call.=FALSE)
}

#################
##Source config##
#################

source(rconfig)
source(rfunctions)
verbose(paste0("Family ID: ", fam))

###################################
##Read files and define variables##
###################################
verbose("Reading annotation files")
hpodb <- fread(hpodb_path)
gdroi <- fread(gd_path, header = F)
pli <- fread(pli_path, header = F)
prec <- fread(prec_path, header = F)
eo <- fread(eo_path, header = F)
genelist <- fread(genelist_path, header = F)
omim <- fread(mim_path, fill=TRUE)

af_columns <- strsplit(af_columns, ",")[[1]]

verbose("Reading manifest and variants")
manifest <- fread(manifest_path, stringsAsFactors = F, header = T)
fam_info <- subset(manifest, family_id == fam)
if( nrow(subset(fam_info, affected == 2)) == 0 ){
  cat("No affected individuals in this family\n")
  quit()
}
vars <- fread(variants_path, header = T, stringsAsFactors = F)
names(vars)[1] <- "chrom"

vars_gt_command <- paste0("zcat ", variants_gt_path, " | grep -v ^##")
vars_gt <- fread(vars_gt_command, header = T, stringsAsFactors = F)
fam_info <- subset(fam_info, subject_id %in% names(vars_gt))
fam_samples <- fam_info$subject_id
fam_struct <- unique(fam_info$family_structure)
# fam_struct <- ifelse(length(fam_samples) == 1, "singleton", 
#                      ifelse(length(fam_samples) == 2, "duo", 
#                             ifelse(length(fam_samples) == 3, "trio", 
#                                    ifelse(length(fam_samples) == 4, "quad", "other") ) ) )
mother <- unique(subset(fam_info, !maternal_id %in% c("", 0))$maternal_id)
mother <- mother[mother %in% names(vars_gt)]
father <- unique(subset(fam_info, !paternal_id %in% c("", 0))$paternal_id)
father <- father[father %in% names(vars_gt)]

affected <- subset(fam_info, affected == 2)$subject_id
affected_all <- subset(manifest, affected == 2)$subject_id
unaffected <- subset(fam_info, affected == 1)$subject_id
unaffected_all <- subset(manifest, affected == 1)$subject_id

##TO DO: add sanity checks

#####################
##Initial filtering##
#####################
verbose("Family and frequency initial filtering")

#keep only if affected sample is a carrier
vars$keep <- apply(vars,1,function(row){
  row_samples <- strsplit(row['samples'], ",")[[1]]
  any(row_samples %in% affected)
})

vars_aff <- subset(vars, keep == TRUE)

#Add specific affected sample carrying the variant
vars_aff$sample <- apply(vars_aff,1,function(row){
 row_samples <- strsplit(row['samples'], ",")[[1]]
 paste(row_samples[row_samples %in% affected], collapse = ",")
})

#Replace NA for 0 in AF and AC columns
vars_aff %>%
  mutate_at(vars(grep("_AF", names(vars_aff), value = T)), ~replace_na(., 0)) -> vars_aff

#Flag if in genomic disorders
if( exists("gdroi") ){
  vars_aff$IN_GD <- vars_aff$name %in% gdroi$V10
}else{
  vars_aff$IN_GD <- ""
}

#Filter rare variants in the AF columns given as input
vars_aff %>%
  filter_at(vars(af_columns), all_vars(. <= 0.01)) -> vars_aff_rare

#Append variants in GDs that are not rare
if( exists("gdroi") ){
  vars_gd_only <- subset(vars_aff, IN_GD == TRUE & !name %in% vars_aff_rare$name)
  vars_aff_rare_gd <- rbind(vars_aff_rare, vars_gd_only)
}else{
  vars_aff_rare_gd <- vars_aff_rare
}


############
##Annotate##
############
verbose("Defining variables for annotation")

#Genotype information
gt_info <- data.frame(subset(vars_gt, ID %in% vars_aff_rare_gd$name, select = c("ID", fam_samples)))
names(gt_info)[1] <- "name"
names(gt_info) <- gsub("X__", "__", names(gt_info))

gt_info[,fam_samples] <- data.frame(apply(data.frame(gt_info[,fam_samples]), 2, function(c) {
	get_sv_gt(c) 
}))

gt_info_fix <- gt_info

##Fix genotypes before merging with main table
gt_info_fix[names(gt_info_fix)[2:ncol(gt_info_fix)]] <- apply(gt_info_fix[names(gt_info_fix)[2:ncol(gt_info_fix)]], 2, function(column) 
  gsub("1/0", "0/1", column))

gt_info_fix[names(gt_info_fix)[2:ncol(gt_info_fix)]] <- apply(gt_info_fix[names(gt_info_fix)[2:ncol(gt_info_fix)]], 2, function(column) 
  gsub("\\.\\/\\.", "0/0", column))

vars_aff_rare_gt <- merge(vars_aff_rare_gd, gt_info_fix, by = "name", all.x = T, all.y = F)

#Define columns that have gene annotations
gene_cols <- grep("PREDICTED_", names(vars_aff_rare_gt), value = T)
gene_cols <- gene_cols[gene_cols %ni% c("PREDICTED_NONCODING_SPAN", "PREDICTED_NONCODING_BREAKPOINT", "PREDICTED_INTERGENIC")]

verbose("Annotating sample counts")

#Family sample count
vars_aff_rare_gt$SC_FAM <- apply(vars_aff_rare_gt, 1, function(row) 
	get_counts(row, samp = fam_samples))

#Family affected sample count
vars_aff_rare_gt$SC_FAM_AFF <- apply(vars_aff_rare_gt, 1, function(row) 
	get_counts(row, samp = affected))

#Family unaffected sample count
vars_aff_rare_gt$SC_FAM_UNAFF <- apply(vars_aff_rare_gt, 1, function(row)
    get_counts(row, samp = unaffected))

#Cohort affected sample count
vars_aff_rare_gt$SC_ALL_AFF <- apply(vars_aff_rare_gt, 1, function(row)
    get_counts(row, samp = affected_all))

#Cohort unaffected sample count
vars_aff_rare_gt$SC_ALL_UNAFF <- apply(vars_aff_rare_gt, 1, function(row)
    get_counts(row, samp = unaffected_all))

verbose("Finding HPO match")

#HPO information
if("hpo_present" %in% names(vars_aff_rare_gt)){
  vars_aff_rare_gt$HPO_MATCH <- apply(vars_aff_rare_gt, 1, function(row)
    get_hpo_match(row, samp = affected))
}else{
  vars_aff_rare_gt$HPO_MATCH <- FALSE
}

verbose("Finding OMIM annotations")

#HPO information
vars_aff_rare_gt$OMIM <- apply(vars_aff_rare_gt, 1, function(row)
	get_omim(row, samp = affected))

verbose("Finding genelist match")

#Gene in genelist 
vars_aff_rare_gt$GENELIST_MATCH <- apply(vars_aff_rare_gt, 1, function(row)
	get_gene_match(row, genes = genelist$V1))

verbose("Annotating for constraint information")

#Gene in eo list
vars_aff_rare_gt$eo_ANY <- apply(vars_aff_rare_gt, 1, function(row)
    get_gene_match(row, genes = eo$V1))

#Gene has pRec >= 0.9
vars_aff_rare_gt$pRec_ANY <- apply(vars_aff_rare_gt, 1, function(row)
    get_gene_match(row, genes = prec$V1))

#vars_aff_rare_gt$VAR <- paste0(vars_aff_rare_gt$chrom, ":", vars_aff_rare_gt$start, "-", vars_aff_rare_gt$end)

#vars_aff_rare_gt <- vars_aff_rare_gt[order(vars_aff_rare_gt$chrom, vars_aff_rare_gt$start, vars_aff_rare_gt$end),]

#svs_gdroi <- bedr(input = list(a = vars_aff_rare_gt$VAR, 
					#b = paste0(gdroi$V1, ":", gdroi$V2, "-", gdroi$V3)), 
	  			  #method = "intersect", 
        		  #params = "-wa -f 0.1",
        		  #verbose = FALSE)

#vars_aff_rare_gt$IN_GD_ROI <- FALSE

#vars_aff_rare_gt$IN_GD_ROI <- apply(vars_aff_rare_gt, 1, function(row)
    #get_gd(row))

#vars_aff_rare_gt[vars_aff_rare_gt$VAR %in% svs_gdroi &
			#vars_aff_rare_gt$SC_ALL_UNAFF <= 5,]$IN_GD_ROI <- TRUE
	

####################
##Inheritance flag##
####################

#1. Absent in unaffected / de novo
##Absent in any unaffected, filter by frequency, gene list, HPO match and/or constraint

verbose("Absent in unaffected flag")
vars_aff_rare_gt$FILT_ABSENT_UNAFF <- FALSE

#filter stricter gnomAD AF
vars_aff_rare_gt %>%
  filter_at(vars(af_columns), all_vars(. <= 1e-3)) -> tmp_vars_aff_rare_gt
vars_aff_rare_gt <- tmp_vars_aff_rare_gt

vars_aff_rare_gt[vars_aff_rare_gt$SC_FAM_UNAFF == 0 &
                        vars_aff_rare_gt$AC <=10 &
                        vars_aff_rare_gt$SC_ALL_UNAFF <= 5
                        # (is.na(vars_aff_rare_gt$`gnomad_v4.1_sv_AF`) | vars_aff_rare_gt$`gnomad_v4.1_sv_AF` <= 1e-3)
                        # (vars_aff_rare_gt$GENELIST_MATCH | vars_aff_rare_gt$HPO_MATCH | vars_aff_rare_gt$eo_ANY ) &
      ,]$FILT_ABSENT_UNAFF <- TRUE


#2. Compound het SV-SV
##If trio, looks for cosegregation - if not trio, returns genes with multiple hits
if(comp_het == "true"){
    verbose("Compound heterozygous flag")
    tmp_vars_aff_rare_gt <- vars_aff_rare_gt
    if(fam_struct %in% c("trio", "quad")){
        vars_aff_rare_gt <- cbind(tmp_vars_aff_rare_gt,
                            do.call(rbind, apply(tmp_vars_aff_rare_gt, 1, function(row)
                                get_comphet_trio(row, tmp_vars_aff_rare_gt, gene_cols)
                            )))

        #vars_aff_rare_gt$FILT_MULT_HIT <- apply(vars_aff_rare_gt, 1, function(row, vars_aff_rare_gt)
            #get_comphet_trio(row) )
    } else {
        #vars_aff_rare_gt$FILT_MULT_HIT <- apply(vars_aff_rare_gt, 1, function(row)
            #get_comphet_singleton(row) )
        vars_aff_rare_gt <- cbind(tmp_vars_aff_rare_gt,
                            do.call(rbind, apply(tmp_vars_aff_rare_gt, 1, function(row)
                                get_comphet_other(row, tmp_vars_aff_rare_gt, gene_cols)
                            )))
    }

    #vars_aff_rare_gt$FILT_MULT_HIT <- FALSE
    # vars_aff_rare_gt[vars_aff_rare_gt$MULT_HIT &
    # 	(vars_aff_rare_gt$GENELIST_MATCH | vars_aff_rare_gt$HPO_MATCH | vars_aff_rare_gt$pRec_ANY),]$FILT_MULT_HIT <- TRUE
    vars_aff_rare_gt$FILT_MULT_HIT <- vars_aff_rare_gt$MULT_HIT
}

#3. Autosomal recessive
##Any affected are hom variant, none unaffected can be hom
verbose("Recessive flag")

vars_aff_rare_gt$FILT_AR <- FALSE

vars_aff_rare_gt$AFF_AR <- apply(vars_aff_rare_gt, 1, function(row) any(row[affected] == "1/1")) 

if(fam_struct != "singleton"){
	vars_aff_rare_gt$UNAFF_NO_AR <- apply(vars_aff_rare_gt, 1, function(row) all(row[unaffected] != "1/1"))
}else{
	vars_aff_rare_gt$UNAFF_NO_AR <- TRUE
}

vars_aff_rare_gt[vars_aff_rare_gt$AFF_AR &
	vars_aff_rare_gt$UNAFF_NO_AR &
	# (vars_aff_rare_gt$GENELIST_MATCH | vars_aff_rare_gt$HPO_MATCH | vars_aff_rare_gt$pRec_ANY) &
	vars_aff_rare_gt$N_HOMALT <= 10,]$FILT_AR <- TRUE

vars_aff_rare_gt$AFF_AR <- NULL
vars_aff_rare_gt$UNAFF_NO_AR <- NULL

#4. X-linked recessive
##Flag FILT_AR if male and variant in chrX
male_probands <- subset(fam_info, !paternal_id %in% c("", "0") & !maternal_id %in% c("", "0") & affected == 2 & sex == 1)$subject_id

if(length(male_probands) > 0 ){
  
  vars_aff_rare_gt$FILT_XLR <- apply(vars_aff_rare_gt, 1, function(row) 
    any(row[male_probands] %in% c("1/1", "1")) & 
      all(row[unaffected] %ni% c("1/1", "1")) & 
      (row["chrom"] == "chrX" | row["CHR2"] == "chrX") 
      # (row["GENELIST_MATCH"] == TRUE | row["HPO_MATCH"] == TRUE)
    )

}else{
  vars_aff_rare_gt$FILT_XLR <- NA
}


##Inherited: 
vars_aff_rare_gt$IN_GOI <- vars_aff_rare_gt$GENELIST_MATCH | vars_aff_rare_gt$eo_ANY | vars_aff_rare_gt$pRec_ANY

vars_aff_rare_gt$FILT_INHERITED <- FALSE

vars_aff_rare_gt[  vars_aff_rare_gt$IN_GOI &
                   vars_aff_rare_gt$AC <= 10 &
                   vars_aff_rare_gt$SC_FAM_UNAFF > 0 &
                   vars_aff_rare_gt$SC_ALL_UNAFF <= 5,]$FILT_INHERITED <- TRUE


#Identify LOF ones to give later higher priority
##Pending? Classify if LOF, UTR, PROMOTER, Intronic, Copy Gain and
lof_cols <- grep("LOF", names(vars_aff_rare_gt), value = T)
vars_aff_rare_gt$IS_LOF <- apply(vars_aff_rare_gt[,..lof_cols], 1, function(r) !all(is.na(r)))
vars_aff_rare_gt$IS_LOF_GENELIST <- apply(vars_aff_rare_gt[,..lof_cols], 1, function(r){any(r %in% genelist$V1)})

verbose("Writting to output")

##Rformat columns before writting
vars_aff_rare_gt$unaffected <- apply(vars_aff_rare_gt[,..unaffected], 1, function(r) paste(r, collapse=","))
vars_aff_rare_gt$affected <- apply(vars_aff_rare_gt[,..affected], 1, function(r) paste(r, collapse=","))

keep_cols <- names(vars_aff_rare_gt)[names(vars_aff_rare_gt) %ni% c(affected, unaffected)]

if(comp_het == "true"){
    vars_out <- subset(vars_aff_rare_gt, IN_GD | FILT_ABSENT_UNAFF | FILT_MULT_HIT | FILT_AR | FILT_XLR | FILT_INHERITED, select = keep_cols )
}else{
    vars_out <- subset(vars_aff_rare_gt, IN_GD | FILT_ABSENT_UNAFF | FILT_AR | FILT_XLR | FILT_INHERITED, select = keep_cols )
}

vars_out$FAMILY <- fam

#Write output
write.table(vars_out, out_file, sep = '\t', quote = F, row.names = F)
