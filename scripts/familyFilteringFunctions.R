#!/usr/bin/env Rscript


##################
##Load libraries##
##################

library("optparse")
library("plyr")
library("data.table")
library("future.apply")
# library("bedr")
library("dplyr")
library("tidyr")

`%ni%` <- Negate(`%in%`)

eo_path <- '/data/talkowski/an436/resources/gnomad/constrained_gnomad_oe_3.5.genes.list'
prec_path <- '/data/talkowski/an436/resources/gnomad/gnomad.v2.1.1.pRec.0.9.txt'
pli_path <- '/data/talkowski/an436/resources/gnomad/gnomad.v2.1.1.pLI.0.9.txt'
# gdroi_path <- '/data/talkowski/an436/projects/CMG/batch1/data/GD_segments_hg38.sort.merge.bed'
hpodb_path <- '/data/talkowski/an436/resources/hpo/latest/phenotype_to_genes.txt'
omim_path <- '/data/talkowski/an436/resources/omim/20201201/genemap2.ref.txt'

#############
##Functions##
#############

verbose <- function(cmd) {
    if ( opt$verbose ) {
    write(paste0(cmd, "\n"), stderr())
}}

get_counts <- function(row, samp){
  length(samp[samp %in% strsplit(row['samples'], ",")[[1]]])
}

get_hpo_match <- function(row, samp){
    rsamp <- samp[samp %in% strsplit(as.character(row['samples']), ",")[[1]]]
	#if(length(rsamp) > 0){
    hpos <- strsplit(subset(fam_info, subject_id %in% rsamp)$hpo_present, ",")[[1]]

    if(length(hpos) > 0){
   	    hpo_genes <- subset(hpodb, HPO_ID %in% hpos)$GENE_SYMBOL
        rgenes <- unlist(strsplit(unlist(row[gene_cols], use.names = FALSE), ","))
        any(rgenes %in% hpo_genes)
    }else{
        FALSE
    }
	#}else{
		#FALSE
	#}
}

get_omim <- function(row, samp){
    rgenes <- unlist(strsplit(unlist(row[gene_cols], use.names = FALSE), ","))
    paste(subset(omim, symbol %in% rgenes)$phenotypes, collapse = ";")
}

get_gene_match <- function(row, genes){
    #rsamp <- samp[samp %in% strsplit(row['samples'], ",")[[1]]]
    rgenes <- unlist(strsplit(unlist(row[gene_cols], use.names = FALSE), ","))
    any(rgenes %in% genes)
}

get_gd <- function(row){
    grepl(row['svtype'], row['type_GD'])
}

get_sv_gt <- function(col) {
	sapply(col, function(row) {
		strsplit(row, ":")[[1]][1]
	})
}

get_comphet_other <- function(row, variants) {
	rgenes <- unlist(strsplit(unlist(row[gene_cols],use.names = FALSE), ","))
	df <- with(variants, variants[grepl(paste(rgenes[!is.na(rgenes)], collapse="|"), get(gene_cols)),])

    if(nrow(df) > 1){
	#ogenes <- unlist(strsplit(unlist(vars_aff_rare_gt[vars_aff_rare_gt$name != row['name'], ..gene_cols], use.names = FALSE), ","))
	#compHet <- any(rgenes[!is.na(rgenes)] %in% ogenes[!is.na(ogenes)])
		IS_COMPHET <- TRUE
		COMPHET_ID <- paste0(df$name, collapse = ",")
	}else{
		IS_COMPHET <- FALSE
		COMPHET_ID <- NA
	}

	data.table("MULT_HIT" = IS_COMPHET, "MULT_HIT_ID" = COMPHET_ID)

}

get_comphet_trio <- function(row, variants) {
	rgenes <- unlist(strsplit(unlist(row[gene_cols],use.names = FALSE), ","))
	rgenes_all <- paste(rgenes[!is.na(rgenes)], collapse="|")

	if(rgenes_all != ""){

    	df <- with(variants, variants[grepl(paste(rgenes[!is.na(rgenes)], collapse="|"), get(gene_cols)),])

			if(nrow(df) > 1){
				#df$affectedMask <- df[,..affected] == "0/1"
				df$affectedMask <- apply(df, 1, function(r) all(r[affected] == "0/1"))
				#df$unaffectedMask <- df[,unaffected] == "0/1"
				df$fatherMask <- df[,..father] == "0/1"
				df$motherMask <- df[,..mother] == "0/1"

				df$compHetFatherMask <- df$affectedMask & !df$motherMask & df$fatherMask
				df$compHetMotherMask <- df$affectedMask & df$motherMask & !df$fatherMask
				df$compHetDenovoMask <- df$FILT_ABSENT_UNAFF

				IS_COMPHET <- any(df$compHetFatherMask) & any(df$compHetDenovoMask) | 
							  any(df$compHetMotherMask) & any(df$compHetDenovoMask) | 
					 		 any(df$compHetFatherMask) & any(df$compHetMotherMask)

				if(IS_COMPHET) {
					COMPHET_ID <- paste0(df$name, collapse = ",")
				}else{
					COMPHET_ID <- NA
				}

				data.table("MULT_HIT" = IS_COMPHET, "MULT_HIT_ID" = COMPHET_ID)
			}else{

				data.table("MULT_HIT" = FALSE, "MULT_HIT_ID" = NA)
			}
	}else{

	data.table("MULT_HIT" = FALSE, "MULT_HIT_ID" = NA)

	}
}