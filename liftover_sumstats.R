#!/usr/bin/env Rscript
################################## 
### Function: liftover from hg19 to hg38 for GWAS summary stats
### date 27-Oct-2023
### version 1.0
### author: EALM (ealopera@gmail.com)
##################################
## Notes
# This script takes GWAS summary statistics and uses liftover to convert from
# genome builds hg19 to hg38. It requires the columns:
# CHR, POS, ALT, REF
# indicated by column names.
######
## New
####### load packages #############
library(data.table)
library(tidyverse)
library(stringr)
library(optparse)
options("scipen"=100, "digits"=10)
############################

#########################################################################################################
option_list = list(
  make_option(c("-w", "--wkdir"), type="character", default=NULL, 
              help="Path were output directory will be created", 
              metavar="character"),
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="path to the summ-stats file", 
              metavar="character"),
  make_option(c("-t", "--tdir"), type="character", default=F, 
              help="path to the liftover executalbe file",
              metavar="character"),
  	make_option("--SNP", type="character", default="snpid",
                help="colnames for SNP ID"),
  make_option(c("-r", "--refdir"), type="character", default=NULL, 
              help="Path to directory for reference hg19ToHg38.over.chain.gz",
              metavar="character"),
  make_option("--CHR", type="character", default="chrom",
              help="colnames for chromosome"),
  make_option("--A1", type="character", default="A1",
              help="colnames for allele1"),
  make_option("--A2", type="character", default="A2",
              help="colnames for allele2"),
  make_option("--POS", type="character", default="pos",
              help="colnames for genome position"),
  make_option(c("-p", "--prefix"), type="character", default="",
              help="prefix for output files")
); 

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

##step controls
if (is.null(opt$wkdir)){
  print_help(opt_parser)
  stop("At least one argument must be supplied as working dir", call.=FALSE)
}
if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied as working file", call.=FALSE)
}

####test ############
# tdir<-"/u/project/loes/elopera/tools/liftOver"
# inpfile<-"/u/project/loes/sservice/GWAS_SumStats/daner_PGC_BIP32b_mds7a.gz"
# wkdir<-"/u/project/loes/elopera/pubGWAS_sumstats/"
# refdir<-"/u/project/loes/elopera/tools/"
# outdir<-"/u/project/loes/elopera/"
# opt$prefix<-"daner_PGC_BIP32b_mds7a"
### deinfe variables
snpcol<- opt$SNP
chrcol <- opt$CHR
poscol <- opt$POS
A1col <- opt$A1
A2col <- opt$A2
wkdir<-opt$wkdir
inpfile<-opt$file
lift<-opt$lift
refdir<-opt$refdir
tdir<-opt$tdir
phename<-opt$prefix
#################################### main ######################################
### load data
input_file<-file.path(wkdir,inpfile)
if( grepl(".gz$",inpfile) | grepl(".bgz$",inpfile) ) {
  sumstats = fread(cmd=paste0("gunzip -c ", inpfile), header=T,fill=TRUE)
} else {
  sumstats <- fread(inpfile, header=T,fill=TRUE)
}
setnames(sumstats, c(snpcol,chrcol, poscol, A1col, A2col), c("SNP","CHR", "BP", "A1", "A2"))
### convert X chromosome snps name
sumstats[sumstats$CHR=="23",1]<-"X"
##name output dir
outfile<-file.path(wkdir,paste0(phename,"_lifted.txt"))
bedfile<-file.path(wkdir,paste0(phename,"_snps_pos.bed"))
bimfile<-file.path(wkdir,paste0(phename,"_hg38remapped.bim"))
nomapfile<-file.path(wkdir,paste0(phename,".unmapped_hg38"))
dir.create(wkdir,recursive = T)

#### liftover usage
  ## create bed file to convert with liftover
  snpcoord<-data.frame(select(sumstats,c( "SNP","CHR",  "BP", "A1","A2")))
  snpcoord$CHR<-paste0("chr",snpcoord$`CHR`)
  snpcoord$BP<-as.numeric(  snpcoord$BP)
  ### handle insertions
  sep1<-str_split_fixed(snpcoord$A1,pattern="",n=2)
  sep2<-str_split_fixed(snpcoord$A2,pattern="",n=2)
  snpcoord$endpos1<-as.numeric(sep1[,2])
  snpcoord$endpos2<-as.numeric(sep2[,2])
  snpcoord$endposdel<-ifelse(is.na(snpcoord$endpos1),snpcoord$endpos2,snpcoord$endpos1)
  snpcoord$BP2<-ifelse(is.na(snpcoord$endposdel),snpcoord$BP+1,snpcoord$BP+snpcoord$endposdel)
 
  snpcoord<-snpcoord[,c("CHR","BP","BP2","SNP","A1","A2","BP")]
  names(snpcoord)<-c("chr","start","end","name","A1","A1","hg37pos")
  
  write.table(snpcoord,bedfile,
              quote = F,row.names = F,sep="\t",col.names = F)
  
  ### do the remmaping using lifover
  lift.call<- paste0( tdir, " -bedPlus=4 ",
                      bedfile," ",
                     refdir,"/hg19ToHg38.over.chain.gz ",
                     bimfile," ",
                      nomapfile)
  system(lift.call)

  ### bring the results of liftover together with the original results 
  remaped<-fread(bimfile,data.table = F )
  head(remaped)
  names(sumstats)
  setnames(sumstats, c("BP","A1","A2"),c("hg19_POS","hg19_A1","hg19_A2"))
  sumstats$POS<-remaped$V2[match(sumstats$SNP,remaped$V4)]
  sumstats$A1<-remaped$V5[match(sumstats$SNP,remaped$V4)]
  sumstats$A2<-remaped$V6[match(sumstats$SNP,remaped$V4)]
unlink(bedfile)
print("[INFO] temporary befile was erased")
write.table(sumstats,file.path(wkdir,paste0(phename,"_lifted.txt")),row.names = F,quote = F)
print("[INFO] Lifted summary stats written in working directory ")
###done##
