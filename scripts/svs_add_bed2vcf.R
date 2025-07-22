##################
##Load libraries##
##################
library(dplyr)
library(readr)
library(data.table)
library(optparse)

#######################
##Define input-output##
#######################
# Define Input Arguments
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Final denovo SVs", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Final denovo SVs prep for annotation", metavar="character"),
  make_option(c("-s", "--samples"), type="character", default=NULL,
              help="Optional samples list", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)

#Define input parameteres
bed_file <- opt$input
vcf_file <- opt$output

if(length(opt$samples) >0){
  sample_list <- opt$samples
}

#############
##FUNCTIONS##
#############

# Function to convert BED to VCF
convert_bed_to_vcf <- function(bed_file, vcf_file, vcf_header) {
  # Read the BED file
  bed_data <- fread(bed_file)
  bed_data$SVLEN <- as.numeric(bed_data$end) - as.numeric(bed_data$start)

  # Extract unique samples
  # samples <- unique(bed_data$sample)
  if(length(opt$samples) >0){
    samples <- fread(sample_list, header = FALSE)$V1
  }

  # Create the header line with samples
  if(length(opt$samples) >0){
    header_line <- paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", paste(samples, collapse = "\t"), sep = "\t")
  }else{
    header_line <- paste("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE", sep = "\t")
  }

  # Write VCF header to the file
  writeLines(c(vcf_header, header_line), vcf_file)

  # Process BED data
  sv_types <- unique(bed_data$SVTYPE)

  bed_data_no_cpx <- bed_data %>%
    subset(!SVTYPE %in% c("CPX", "CTX")) %>%
    mutate(
      REF = "N",
      ALT = ifelse(SVTYPE == "MULTIPLE", "<INS>", paste0("<", SVTYPE, ">")),
      QUAL = ".",
      FILTER = "PASS",
      INFO = paste0("END=", end, ";SVTYPE=", SVTYPE, ";ALGORITHMS=", name, ";EVIDENCE=MANUAL;SVLEN=", SVLEN),
      FORMAT = "GT:GQ:RD_CN:RD_GQ:PE_GT:PE_GQ:SR_GT:SR_GQ:EV:CN",
    )
  bed_data <- bed_data_no_cpx

  if("CPX" %in% sv_types){
    bed_data_cpx <- bed_data %>%
      subset(SVTYPE == "CPX") %>%
      mutate(
        REF = "N",
        ALT = ifelse(CPX_TYPE_VCF == "UNRESOLVED", "<BND>", paste0("<", SVTYPE, ">")),
        QUAL = ".",
        FILTER = ifelse(CPX_TYPE_VCF == "UNRESOLVED", "UNRESOLVED", "PASS"),
        INFO = ifelse(CPX_TYPE_VCF == "UNRESOLVED",
          paste0("END=", end, ";SVTYPE=BND;ALGORITHMS=", name, ";EVIDENCE=MANUAL;SVLEN=-1"),
          paste0("END=", end, ";SVTYPE=", SVTYPE, ";ALGORITHMS=", name, ";EVIDENCE=MANUAL;SVLEN=", SVLEN, ";CPX_INTERVALS=", CPX_INTERVALS_VCF, ";CPX_TYPE=", CPX_TYPE_VCF)
        ),
        FORMAT = "GT:GQ:RD_CN:RD_GQ:PE_GT:PE_GQ:SR_GT:SR_GQ:EV:CN",
      )
    bed_data <- rbind(bed_data, bed_data_cpx)
  }

  if("CTX" %in% sv_types){
    bed_data_ctx <- bed_data %>%
      subset(SVTYPE == "CTX") %>%
      mutate(
        REF = "N",
        ALT = paste0("<", SVTYPE, ">"),
        QUAL = ".",
        FILTER = "PASS",
        INFO = paste0("END=", end, ";SVTYPE=CTX;ALGORITHMS=", name, ";EVIDENCE=MANUAL;SVLEN=1;CHR2=", CHR2, ";END2=", END2),
        FORMAT = "GT:GQ:RD_CN:RD_GQ:PE_GT:PE_GQ:SR_GT:SR_GQ:EV:CN",
      )
    bed_data <- rbind(bed_data, bed_data_ctx)
  }

  # bed_data <- rbind(bed_data_no_cpx, bed_data_cpx, bed_data_ctx)
  # bed_data <- rbind(bed_data_no_cpx, bed_data_cpx)

  # Prepare the genotype columns for each sample
  if(length(opt$samples) >0){
    for (sample in samples) {
      bed_data[[sample]] <- ifelse(bed_data$sample == sample, "0/1:.:.:.:.:.:.:.:.:2", "0/0:.:.:.:.:.:.:.:.:2")
    }
  }else{
    bed_data[["SAMPLE"]] <- "0/1:.:.:.:.:.:.:.:.:2"
  }

  # Select the columns to write to the VCF file
  if(length(opt$samples) >0){
    vcf_cols <- c("chrom", "start", "name", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples)
  }else{
    vcf_cols <- c("chrom", "start", "name", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")
  }
  bed_data_write <- bed_data %>% select(all_of(vcf_cols))

  # Write VCF entries to the file
  write_delim(bed_data_write, vcf_file, delim = "\t", append = TRUE, col_names = FALSE)

}

# Define the VCF header
vcf_header <- c(
  "##fileformat=VCFv4.2",
  "##FILTER=<ID=PASS,Description=\"All filters passed\">",
  "##ALT=<ID=CNV,Description=\"Copy Number Polymorphism\">",
  "##ALT=<ID=CPX,Description=\"Complex SV\">",
  "##ALT=<ID=CTX,Description=\"Reciprocal chromosomal translocation\">",
  "##ALT=<ID=INS,Description=\"Insertion\">",
  "##ALT=<ID=INS:ME,Description=\"Mobile element insertion of unspecified ME class\">",
  "##ALT=<ID=INS:ME:ALU,Description=\"Alu element insertion\">",
  "##ALT=<ID=INS:ME:LINE1,Description=\"LINE1 element insertion\">",
  "##ALT=<ID=INS:ME:SVA,Description=\"SVA element insertion\">",
  "##ALT=<ID=INS:UNK,Description=\"Sequence insertion of unspecified origin\">",
  "##CPX_TYPE_INS_iDEL=\"Insertion with deletion at insertion site.\"",
  "##CPX_TYPE_INVdel=\"Complex inversion with 3' flanking deletion.\"",
  "##CPX_TYPE_INVdup=\"Complex inversion with 3' flanking duplication.\"",
  "##CPX_TYPE_dDUP=\"Dispersed duplication.\"",
  "##CPX_TYPE_dDUP_iDEL=\"Dispersed duplication with deletion at insertion site.\"",
  "##CPX_TYPE_delINV=\"Complex inversion with 5' flanking deletion.\"",
  "##CPX_TYPE_delINVdel=\"Complex inversion with 5' and 3' flanking deletions.\"",
  "##CPX_TYPE_delINVdup=\"Complex inversion with 5' flanking deletion and 3' flanking duplication.\"",
  "##CPX_TYPE_dupINV=\"Complex inversion with 5' flanking duplication.\"",
  "##CPX_TYPE_dupINVdel=\"Complex inversion with 5' flanking duplication and 3' flanking deletion.\"",
  "##CPX_TYPE_dupINVdup=\"Complex inversion with 5' and 3' flanking duplications.\"",
  "##CPX_TYPE_piDUP_FR=\"Palindromic inverted tandem duplication, forward-reverse orientation.\"",
  "##CPX_TYPE_piDUP_RF=\"Palindromic inverted tandem duplication, reverse-forward orientation.\"",
  "##FILTER=<ID=BOTHSIDES_SUPPORT,Description=\"Variant has read-level support for both sides of breakpoint\">",
  "##FILTER=<ID=HIGH_SR_BACKGROUND,Description=\"High number of SR splits in background samples indicating messy region\">",
  "##FILTER=<ID=MULTIALLELIC,Description=\"Multiallelic site\">",
  "##FILTER=<ID=PESR_GT_OVERDISPERSION,Description=\"High PESR dispersion count\">",
  "##FILTER=<ID=UNRESOLVED,Description=\"Variant is unresolved\">",
  "##FORMAT=<ID=CN,Number=1,Type=Integer,Description=\"Predicted copy state\">",
  "##FORMAT=<ID=CNQ,Number=1,Type=Integer,Description=\"Read-depth genotype quality\">",
  "##FORMAT=<ID=EV,Number=.,Type=String,Description=\"Classes of evidence supporting final genotype\">",
  "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",
  "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
  "##FORMAT=<ID=PE_GQ,Number=1,Type=Integer,Description=\"Paired-end genotype quality\">",
  "##FORMAT=<ID=PE_GT,Number=1,Type=Integer,Description=\"Paired-end genotype\">",
  "##FORMAT=<ID=RD_CN,Number=1,Type=Integer,Description=\"Predicted copy state\">",
  "##FORMAT=<ID=RD_GQ,Number=1,Type=Integer,Description=\"Read-depth genotype quality\">",
  "##FORMAT=<ID=SR_GQ,Number=1,Type=Integer,Description=\"Split read genotype quality\">",
  "##FORMAT=<ID=SR_GT,Number=1,Type=Integer,Description=\"Split-read genotype\">",
  "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">",
  "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Source algorithms\">",
  "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">",
  "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Second contig\">",
  "##INFO=<ID=CPX_INTERVALS,Number=.,Type=String,Description=\"Genomic intervals constituting complex variant.\">",
  "##INFO=<ID=CPX_TYPE,Number=1,Type=String,Description=\"Class of complex variant.\">",
  "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">",
  "##INFO=<ID=END2,Number=1,Type=Integer,Description=\"Second position\">",
  "##INFO=<ID=EVIDENCE,Number=.,Type=String,Description=\"Classes of random forest support.\">",
  "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Source of inserted sequence.\">",
  "##INFO=<ID=STRANDS,Number=1,Type=String,Description=\"First and second strands\">",
  "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of affected segment on the reference\">",
  "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
  "##INFO=<ID=UNRESOLVED_TYPE,Number=1,Type=String,Description=\"Class of unresolved variant.\">",
  "##contig=<ID=chr1,assembly=38,length=248956422>",
  "##contig=<ID=chr2,assembly=38,length=242193529>",
  "##contig=<ID=chr3,assembly=38,length=198295559>",
  "##contig=<ID=chr4,assembly=38,length=190214555>",
  "##contig=<ID=chr5,assembly=38,length=181538259>",
  "##contig=<ID=chr6,assembly=38,length=170805979>",
  "##contig=<ID=chr7,assembly=38,length=159345973>",
  "##contig=<ID=chr8,assembly=38,length=145138636>",
  "##contig=<ID=chr9,assembly=38,length=138394717>",
  "##contig=<ID=chr10,assembly=38,length=133797422>",
  "##contig=<ID=chr11,assembly=38,length=135086622>",
  "##contig=<ID=chr12,assembly=38,length=133275309>",
  "##contig=<ID=chr13,assembly=38,length=114364328>",
  "##contig=<ID=chr14,assembly=38,length=107043718>",
  "##contig=<ID=chr15,assembly=38,length=101991189>",
  "##contig=<ID=chr16,assembly=38,length=90338345>",
  "##contig=<ID=chr17,assembly=38,length=83257441>",
  "##contig=<ID=chr18,assembly=38,length=80373285>",
  "##contig=<ID=chr19,assembly=38,length=58617616>",
  "##contig=<ID=chr20,assembly=38,length=64444167>",
  "##contig=<ID=chr21,assembly=38,length=46709983>",
  "##contig=<ID=chr22,assembly=38,length=50818468>",
  "##contig=<ID=chrX,assembly=38,length=156040895>",
  "##contig=<ID=chrY,assembly=38,length=57227415>"
)

#########################
##RUN CONVER BED TO VCF##
#########################
convert_bed_to_vcf(bed_file, vcf_file, vcf_header)