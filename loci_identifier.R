#!/usr/bin/env Rscript
################################## 
### Function: Loci identifier for GWAS summary statistics
### date 10-Oct-2023
### version 1.0
### author: EALM (ealopera@gmail.com)
##################################
## Notes
# This script takes GWAS summary statistics and returns a summary of: 
# the hits (hits: pval<5e-4), 
# the loci (loci: regions > [RegionFlank] base pairs),
# and the significant hits (pval < [sigThreshold]) found per loci. 
# Req. columns for input: CHRomosome, position , ALLELE1, ALLELE2, PVAL
# indicated by column name.
# The algorithm spans a region around each [sigThreshold] hit, merging together
# the significant hits around [RegionFlank] base pairs

## New
## 16-Oct-2023 
# future iterations might accept rscodes and include an additional LD filter
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
	make_option(c("-i", "--input"), type="character", default="",
                help="full path to input file"),
	make_option(c("-p", "--prefix"), type="character", default="",
                help="prefix for output files"),
	make_option("--pheno", type="character", default="",
                help="phenotype"),
	make_option("--SNP", type="character", default="snpid",
                help="colnames for SNP ID"),
	make_option("--CHR", type="character", default="chrom",
		help="colnames for chromosome"),
	make_option("--ALLELE1", type="character", default="Allele1",
                help="colnames for allele1"),
	make_option("--ALLELE2", type="character", default="Allele2",
                help="colnames for allele2"),
	make_option("--POS", type="character", default="pos",
                help="colnames for genome position"),
	make_option("--PVAL", type="character", default="",
                help="p values"),
	make_option("--sigThreshold", type="numeric", default=5E-8,
                help="significance threshold for p values [default=5E-8]"),
#	make_option("--knownGWASList", type="character", default="",
#	help="known GWAS list provided by the user"),
	make_option("--RegionFlank", type="numeric", default=200000,
                help="known region flank [default=200000]")
)
opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
print(opt)
## test & developing data
    # opt<-list()
    # opt$sigThreshold<-5E-8
    # opt$CHR<-"chromosome"
    # opt$POS<-"base_pair_location"
    # opt$ALLELE2<-"effect_allele"
    # opt$ALLELE1<-"other_allele"
    # opt$pheno<-"atshma"
    # opt$PVAL<-"p_value"
    # opt$RegionFlank<-200000
    # opt$pheno<-"asthma"
    # opt$input<-"29273806-GCST006862-EFO_0000270.h.tsv.gz"
    # opt$prefix<-paste0(date_label,"_",opt$pheno)
    # #snprs<-"variant_id"
##################################
#### main ####
##################################
## define variables
sigThreshold <- opt$sigThreshold
yLine <- -log10(sigThreshold)
chrcol <- opt$CHR
poscol <- opt$POS
A1col <- opt$ALLELE1
A2col <- opt$ALLELE2
pheno <- opt$pheno
pvalcol <- opt$PVAL
flank <- as.numeric(opt$RegionFlank)
prefix<-opt$prefix
input<-opt$input
date_label<- paste(strsplit(as.character(date()) 
                            ,' ',fixed=TRUE)[[1]][c(3,2,5)],collapse="")
outfile <- paste0(gsub(".+/([^/]+$)","\\1",prefix),"_",date_label,"_",pheno)
file_summarytable <- paste0(outfile,"_loci.txt")
file_tophits <- paste0(outfile,"_hitspE-4.txt")

## import summary statistics
if( grepl(".gz$",input) | grepl(".bgz$",input) ) {
        resIn_1 = fread(cmd=paste0("gunzip -c ", input), header=T, select=c(chrcol, poscol, A1col, A2col,pvalcol))
} else {
        resIn_1 <- fread(input, header=T, select=c(chrcol, poscol, A1col, A2col,pvalcol))
}
print("[INFO] summary statistics file imported, column to be used:")
print( names(resIn_1))

## format summary stats
resIn = resIn_1[, c(chrcol, poscol, A1col, A2col,pvalcol), with=FALSE]
setnames(resIn, c(chrcol, poscol, A1col, A2col,pvalcol), c("CHR", "BP", "ALLELE1", "ALLELE2","PVAL"))
print("OK1")
resIn$ALLELE1 = toupper(resIn$ALLELE1)
resIn$ALLELE2 = toupper(resIn$ALLELE2)
resIn$SNP = paste0(resIn$CHR,"_", resIn$BP, "_", resIn$ALLELE1,"/",resIn$ALLELE2)
resIn$CHR <- gsub("chr","",resIn$CHR)
resIn$CHR <- as.numeric(gsub("X","23",resIn$CHR))

## needed second position to create data.table interval
resIn$BP2 <- resIn$BP
Nmarkers = nrow(resIn)
resIn$PVAL = as.numeric(resIn$PVAL)

## control for extreme pvalues
resIn$PVAL[which(resIn$PVAL < 1*10^-300)] = 1*10^-300
print(min(resIn$PVAL))
print(max(resIn$PVAL))
resIn$log10P = -log10(resIn$PVAL)
print("OK2")
table(resIn$CHR)
which(is.na(resIn$CHR))

## find the significant hits
	tophits <- which(resIn$log10P > -log10(sigThreshold))
	if(length(tophits)>0){
		tophits <- resIn[tophits,]
        	tophits$numCHR <- as.numeric(gsub("X","23",tophits$CHR))
        	x <- as.numeric(tophits$BP)
        	y <- tophits$numCHR
## define rows of start for loci
## the starting positions of loci defined by
## the first row then the next region if chromosome changes or 
## if distance betwen hits gets bigger than the flank variable ($flank)
        	start = c(1, which(diff(y) != 0 | abs(diff(x)) >= flank) + 1) 
## all end of regions defined by $start and the last row    	
        	end = c(start - 1, length(x))	
    
## snippet to use when known gwas regions are added
		highlightRegions <- data.frame(
			'Nearest_Gene'=NA,
			'CHROM'=tophits$CHR[start],
			'START'=tophits$BP[start] - flank ,
			'END'=tophits$BP[end] + flank,
			'STATUS'="potentially_novel",
			'TOP_SNP'=NA,
			'minPVALUE'=NA,
			'Reported_GWAS_Catalog'=NA,
			'Reported_closeGen'=NA
		)
		print(highlightRegions)

## Extract top SNPs from each hit region
		hits <- data.table(
			CHROM=highlightRegions$CHROM,
			BEGIN=highlightRegions$START,
			END=highlightRegions$END, key = c("CHROM", "BEGIN", "END"))
		print(hits)
		## scan the region for all the hits and build summaries
		regionHITS = NULL
		Locisummary = NULL
		for(a in 1:dim(highlightRegions)[1]){
		 ## bring the rows with SNPs within the region
			withinKnownRegion <- !is.na(foverlaps(resIn, hits[a,], 
			by.x=c("CHR", "BP", "BP2"),
			by.y=c("CHROM", "BEGIN", "END"),
			type="any", which=TRUE, mult="first"))
			regionhits <- resIn[which(withinKnownRegion & resIn$log10P > -log10(1E-4)),]
			regionhits$CHROMPOS = paste0(regionhits$CHR, ":", regionhits$BP)
			locusid<-paste0(pheno,"_Locus",	str_pad(a, 4, pad = "0"))
			regionhits$LOCUSID = locusid
			regionHITS = rbind(regionHITS, regionhits)
			##summary of the locus
	    nvariants<-length(which(withinKnownRegion))
	    nregion_hits<-nrow(regionhits)
	    nregion_gwas<-length(which (regionhits$PVAL<sigThreshold))
	    max_hit<-regionhits$SNP[which.max(regionhits$log10P)]
	    min_pval<-min(regionhits$PVAL)
	    locusrow<-cbind(locusid,nvariants,nregion_hits,nregion_gwas,max_hit,min_pval)
	    Locisummary<-rbind(Locisummary,locusrow)
		}
	###write loci summary	
	Locisummary<-data.frame(Locisummary)
	Locisummary<-data.frame(cbind(Locisummary[,1],hits,Locisummary[,-1]))
	colnames(Locisummary)<-c("LOCUS",names(hits),"N_VARIANTS","N_HITS","N_SIGNIF","MOST_SIG_SNP","MIN_PVAL")
	write.table(Locisummary, file_summarytable, col.names=T, row.names=F, quote=F, sep="\t")
	###write summary of hits by loci
	nrow(regionHITS)
	regionHitSumm<-regionHITS[,c("CHR","BP","ALLELE1","ALLELE2","PVAL","SNP","log10P","CHROMPOS","LOCUSID")]
	write.table(regionHitSumm, file_tophits, col.names=T, row.names=F, quote=F, sep="\t")
	print(paste0("[INFO] The files ",file_summarytable," and ", file_tophits, "were written in the current directory") )
	print("[INFO] done")
	} else {
	  print(paste0("[INFO] No significant hits were found at the pval < ",sigThreshold," threshold" ))
	}
### done
