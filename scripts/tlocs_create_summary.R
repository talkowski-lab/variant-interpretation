#!/usr/bin/env Rscript
#Author: Alba Sanchis-Juan
##This script creates the tlocs summary table

#Load libraries
library(data.table)
library(optparse)

##Arguments
option_list = list(
  make_option(c("-s", "--samples"), type="character", default=NULL,
              help="samples file name", metavar="character"),
  make_option(c("-i", "--info"), type="character", default=NULL,
              help="info name", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=".",
              help="output file name", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#Read files
samples <- fread(opt$samples)
info <- fread(opt$info)
out <- opt$out

#Create summary table
summary_df <- do.call(rbind, apply(samples, 1, function(row){
  sample_ID <- row['V6']
  tloc_ID <- row['V5']

  info_sample <- subset(info, V7 == sample_ID)

  BP1_1_plus <- min(subset(info_sample, V3 == "+")$V2)
  BP1_2_plus <- max(subset(info_sample, V3 == "+")$V2)
  BP1_1_minus <- min(subset(info_sample, V3 == "-")$V2)
  BP1_2_minus <- max(subset(info_sample, V3 == "-")$V2)
  BP1_int_plus <- BP1_1_plus-BP1_2_plus
  BP1_int_minus <- BP1_1_minus-BP1_2_minus
  BP2_1_plus <- min(subset(info_sample, V6 == "+")$V5)
  BP2_2_plus <- max(subset(info_sample, V6 == "+")$V5)
  BP2_1_minus <- min(subset(info_sample, V6 == "-")$V5)
  BP2_2_minus <- max(subset(info_sample, V6 == "-")$V5)
  BP2_int_plus <- BP2_1_plus-BP2_2_plus
  BP2_int_minus <- BP2_1_minus-BP2_2_minus

  First_PE_dir <- info_sample[order(info_sample$V2),]$V3[1]
  PE <- nrow(info_sample)
  Count1 <- nrow(subset(info_sample, V3 == "+" & V6 == "+"))
  Count2 <- nrow(subset(info_sample, V3 == "+" & V6 == "-"))
  Count3 <- nrow(subset(info_sample, V3 == "-" & V6 == "+"))
  Count4 <- nrow(subset(info_sample, V3 == "-" & V6 == "-"))

  if(info_sample[order(info_sample$V2),]$V3[1] == "-"){
    count <- 0
    vector <- info_sample[order(info_sample$V2),]$V3
    for (i in 1:(length(vector) - 1)) {
      if (vector[i] == "-" && vector[i + 1] == "-") {
        count <- count + 1
      } else if (vector[i] == "+") {
        count <- count + 1
        break  # Interrupt loop if "+" is encountered
      }
    }
    # Output the result
    Count_ins_PE <- count

  }else{
    Count_ins_PE <- 0
  }

  return(data.frame(BP1_1_plus, BP1_2_plus, BP1_1_minus, BP1_2_minus, BP1_int_plus, BP1_int_minus,
                    BP2_1_plus, BP2_2_plus, BP2_1_minus, BP2_2_minus, BP2_int_plus, BP2_int_minus,
                    First_PE_dir, PE, Count1, Count2, Count3, Count4, Count_ins_PE))

}))

#Append to main samples table
samples_summary <- cbind(samples, summary_df)

##Replace negative values to absolute
replace_abs <- function(x) {
  if (is.numeric(x)) {
    return(abs(x))
  } else {
    return(x)
  }
}

samples_summary_abs <- as.data.frame(lapply(samples_summary, replace_abs))

##Add pass counts
samples_summary_abs$pass1 <- lapply(1:nrow(samples_summary_abs), function(i) if ((as.numeric(samples_summary_abs$BP1_int_plus[i])) > 50) {"1"} else {"0"})
samples_summary_abs$pass2 <- lapply(1:nrow(samples_summary_abs), function(i) if ((as.numeric(samples_summary_abs$BP1_int_minus[i])) > 50) {"1"} else {"0"})
samples_summary_abs$pass3 <- lapply(1:nrow(samples_summary_abs), function(i) if ((as.numeric(samples_summary_abs$BP2_int_plus[i])) > 50) {"1"} else {"0"})
samples_summary_abs$pass4 <- lapply(1:nrow(samples_summary_abs), function(i) if ((as.numeric(samples_summary_abs$BP2_int_minus[i])) > 50) {"1"} else {"0"})
samples_summary_abs$Total_pass<-  lapply(1:nrow(samples_summary_abs), function(i) ( sum(as.numeric(unlist(samples_summary_abs$pass1[i])), as.numeric(unlist(samples_summary_abs$pass2[i])),
                                                                                as.numeric(unlist(samples_summary_abs$pass3[i])), as.numeric(unlist(samples_summary_abs$pass4[i])))))

#samples_summary_abs <- apply(samples_summary_abs,2,as.character)

#Write output
fwrite(samples_summary_abs, out, sep = "\t", quote = F, row.names = F, col.names = F)
