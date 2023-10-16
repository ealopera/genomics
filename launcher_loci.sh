##This script will generate regions ${outprefix}.regions.txt which contains merged regions with ${RegionFlank} bp up- and down-stream top hits

Rscript loci_identifier.R	\
--input="29273806-GCST006862-EFO_0000270.h.tsv.gz"	\
--PVAL="p_value"	\
--RegionFlank=500000	\
--prefix=${outprefix}	\
--sigThreshold=5E-8	\
--pheno="asthma" \
--CHR="chromosome"	\
--POS="base_pair_location"	\
--ALLELE1="other_allele"	\
--ALLELE2="effect_allele"
