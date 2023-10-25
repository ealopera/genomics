##This script will generate regions ${outprefix}.regions.txt which contains merged regions with ${RegionFlank} bp up- and down-stream top hits

### stablish positional variables
## input: variable 1: full path to input file
## pheno: variable 2 phenotype name
## SNP var3 colname for SNP ID
## CHR var4 colname for chromosome
## PVAL: var5 colname for p value
## POS: var6 colname for chromosome position
## ALLELE1: var7 colname for allele1
## ALLELE2: var8 colname for allele2
## RegionFlank  var9  known region flank to stablish loci [default=200000]
## sigThreshold var10 significance threshold for p values [default=5E-8]
## prefix: var11 prefix for output files
input=$1
pheno=$2
SNP=$3
CHR=$4	
PVAL=$5	
POS=$6	
ALLELE2=$7
ALLELE1=$8
RegionFlank=$9	
sigThreshold=${10}
prefix=${11}
### failsafe and default variables
if [ -z "$1" ]
then
      echo "input variable is empty"
      exit 1
fi

if [ -z "$2" ]
then
      echo "pheno variable is empty"
      exit 1
fi

if  [ -z "$9" ]
then
      RegionFlank=200000
fi

if  [ -z "{$10}" ]
then
      RegionFlank=5E-8
fi
#load R
ml RPlus
############################# main #############################################

Rscript genomics/loci_identifier.R	\
--input=$input	\
--PVAL=$PVAL	\
--RegionFlank=$RegionFlank	\
--sigThreshold=$sigThreshold	\
--pheno=$pheno \
--CHR=$CHR	\
--POS=$POS	\
--ALLELE1=$ALLELE1	\
--ALLELE2=$ALLELE2 \
--prefix=$prefix \
--SNP=$SNP

######test run
#Rscript  genomics/loci_identifier.R	\
#--input="29273806-GCST006862-EFO_0000270.h.tsv.gz"	\
#--PVAL="p_value"	\
#--RegionFlank=200000	\
#--sigThreshold=5E-8	\
#--pheno="asthma" \
#--CHR="chromosome"	\
#--POS="base_pair_location"	\
#--ALLELE1="other_allele"	\
#--ALLELE2="effect_allele"
