#!/usr/bin/env Rscript
################################## 
### Function: SNP ID harmonizer for GWAS results
### date 25-Oct-2023
### version 1.0
### author: EALM (ealopera@gmail.com)
##################################
## Notes
# This script takes GWAS summary statistics and a list of SNPs as reference
# both need to contain these columns: 
# CHR, POS, ALT, REF
# the target file needs to contain additional Pval and optional SNPname columns
# The harmonizer will create SNP ID column harmonized with the reference list
# and output a file ready to be used by plink
# the harmonize ID will be in format CHR_POS_A1_A2
# Req. columns for input: rCHR, rPOS , rA1, rA2, CHR, POS , A1, A2,
# indicated by column names.
## New
##################################
#### set up ####
options(stringsAsFactors=F)
Sys.setlocale("LC_CTYPE", "C.UTF-8")

req_packages <- c("optparse", "data.table","tidyverse")
for (pack in req_packages) {
  if(!require(pack, character.only = TRUE)) {
    #install.packages(pack, repos = "https://cloud.r-project.org")
    install.packages(pack, repos='http://cran.us.r-project.org')
  }
}
## bash parser
option_list <- list(
  make_option(c("-r", "--refinput"), type="character", default="",
              help="full path to reference input file"),
  make_option(c("-i","--input"), type="character", default="",
              help="input file to be harmonized"),
  make_option("--rCHR", type="character", default="rchrom",
              help="colnames for chromosome in the reference"),
  make_option("--rA1", type="character", default="rA1",
              help="colnames for allele1 in the reference"),
  make_option("--rA2", type="character", default="rA2",
              help="colnames for allele2 in the reference"),
  make_option("--rPOS", type="character", default="rpos",
              help="colnames for genome position in the reference"),
  make_option("--CHR", type="character", default="chrom",
              help="colnames for chromosome"),
  make_option("--A1", type="character", default="A1",
              help="colnames for allele1"),
  make_option("--A2", type="character", default="A2",
              help="colnames for allele2"),
  make_option("--POS", type="character", default="pos",
              help="colnames for genome position"),
  make_option("--PVAL", type="character", default="P",
              help="colnames for p value of target file"),
  make_option("--SNP", type="character", default="SNP",
              help="colnames for SNP identifier  (rs code) in the target file"),
  make_option(c("-p", "--prefix"), type="character", default="",
              help="prefix for output files")
)
opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt)
## test & developing data
#  opt<-list()
# opt$CHR<-"CHR"
# opt$POS<-"BP"
# opt$A2<-"A2"
# opt$A1<-"A1"
# opt$rCHR<-"CHROM"
# opt$rPOS<-"POS"
# opt$rA2<-"ALT"
# opt$rA1<-"REF"
# opt$input<-"PGC_BP_5e-04.txt"
# opt$refinput<-"genotype_SNPlist.txt"
# opt$prefix<-"BDpaisa"
# opt$PVAL<-"P"
# opt$SNP<-"SNP"
##################################
#### main ####
##################################
## define variables
chrcol <- opt$CHR
poscol <- opt$POS
A1col <- opt$A1
A2col <- opt$A2
rchrcol <- opt$rCHR
rposcol <- opt$rPOS
rA1col <- opt$rA1
rA2col <- opt$rA2
prefix<-opt$prefix
input<-opt$input
refinput<-opt$refinput
pvalcol<-opt$PVAL
snpcol<-opt$SNP

outfile <- paste0(gsub(".+/([^/]+$)","\\1",prefix),"_SNPharmo")
file_result <- paste0(outfile,".txt")
file_report <- paste0(outfile,"_report.txt")

## import objective file
if( grepl(".gz$",input) | grepl(".bgz$",input) ) {
  dat1_1 = fread(cmd=paste0("gunzip -c ", input), header=T, select=c(snpcol,pvalcol,chrcol, poscol, A1col, A2col))
} else {
  dat1_1 <- fread(input, header=T, select=c(snpcol,pvalcol,chrcol, poscol, A1col, A2col))
}
## import reference file
if( grepl(".gz$",input) | grepl(".bgz$",refinput) ) {
  datref <-fread(cmd=paste0("gunzip -c ", reefinput), header=T, select=c(rchrcol, rposcol, rA1col, rA2col))
} else {
  datref<- fread(refinput, header=T, select=c(rchrcol, rposcol, rA1col, rA2col))
}

print("[INFO] ref and input files imported,  with columns:")
print(names(dat1_1))
print("[INFO] ref and input files imported,  with columns:")
print(names(datref))

## format SNP data

dat1 = dat1_1[, c(snpcol,pvalcol,chrcol, poscol, A1col, A2col), with=FALSE]
setnames(dat1, c(snpcol,pvalcol,chrcol, poscol, A1col, A2col), c("SNPname","P","CHR", "BP", "A1", "A2"))
setnames(datref, c(rchrcol, rposcol, rA1col, rA2col), c("rCHR", "rBP", "rA1", "rA2"))
print("OK1")

dat1$A1 = toupper(dat1$A1)
dat1$A2 = toupper(dat1$A2)
datref$rA1 = toupper(datref$rA1)
datref$rA2 = toupper(datref$rA2)
dat1$CHR <- gsub("chr","",dat1$CHR)
datref$CHR <- gsub("chr","",datref$CHR)
dat1$CHR <- as.numeric(gsub("X","23",dat1$CHR))
#build equal positions column 
dat1$CHR_POS<-paste0(dat1$CHR,"_",dat1$BP)
datref$CHR_POS<-paste0(datref$rCHR,"_",datref$rBP)
#cut reference to the SNPs contained in data only
inref<-which(datref$CHR_POS %in% dat1$CHR_POS )
datref<-datref[inref,]
# bring data from ref
dat1$rA1<-datref$rA1[match(dat1$CHR_POS,datref$CHR_POS)]
dat1$rA2<-datref$rA2[match(dat1$CHR_POS,datref$CHR_POS)]
dat1$SNP<-case_when(dat1$rA1==dat1$A1 & dat1$rA2==dat1$A2~ paste0(dat1$CHR_POS,"_",dat1$A1,"_",dat1$A2),
                        dat1$rA1==dat1$A2 & dat1$rA2==dat1$A1~ paste0(dat1$CHR_POS,"_",dat1$A2,"_",dat1$A1),
                        dat1$rA1==dat1$A1 & dat1$rA2!=dat1$A2~ paste0(dat1$CHR_POS,"_",dat1$A1,"_",dat1$A2,"b"),
                        dat1$rA2==dat1$A1 & dat1$rA1!=dat1$A2~ paste0(dat1$CHR_POS,"_",dat1$A1,"_",dat1$A2,"b"))

dat1$SNPchange<-case_when(dat1$rA1==dat1$A1 & dat1$rA2==dat1$A2~ "no_change",
                        dat1$rA1==dat1$A2 & dat1$rA2==dat1$A1~ "swapped",
                        dat1$rA1==dat1$A1 & dat1$rA2!=dat1$A2~ "possible_triallelic",
                        dat1$rA2==dat1$A1 & dat1$rA1!=dat1$A2~ "possible_triallelic")

report<-table(dat1$SNPchange)
print("finished harmonization, found:")
print(report)
write.table(dat1,file_result,quote=F,sep='\t',row.names = F)
write.table(report,file_report,quote=F,sep='\t',row.names = F)
print("Harmonized file saved in current directory")
##done###
