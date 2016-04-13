### Manhattan plots of Outliers from LFMM and Arlequin
###    using this link: http://www.gettinggeneticsdone.com/2014/05/qqman-r-package-for-qq-and-manhattan-plots-for-gwas-results.html
### Carolyn Tarpey | April 2016
### ---------------------------------------

#install.packages("qqman")
library(qqman)
vignette("qqman")

setwd('G:/Analysis/Mapping/Outliers_Manhattan_Plot/Arlequin_LFMM')


outliers = read.table('Outliers_Map.txt', header = TRUE)
as.data.frame(outliers)

head(outliers)
table(outliers$CHR)



#pdf("plots/LFMMresults_Manhattan_plots.pdf", width = 8, height = 5)


#	odd_0.01	even_0.01	overlap_0.01	overlap_PT1	overlap_PT2	overlap_LO	overlap_LA

###############################################
#LFMM
# odd pvalues: o_p_PT1	o_p_PT2	o_p_LO	o_p_LA	
# outlisers: odd_PT1	odd_PT2	odd_LO	odd_LA	odd_any	

plot = manhattan(outliers, p = "o_p_PT1", chr = "CHR", bp = "BP", snp = "SNP", main= "Odd Lineage PT1 Outliers")
plot = manhattan(outliers, p = "o_p_PT2", chr = "CHR", bp = "BP", snp = "SNP", main= "Odd Lineage PT2 Outliers")
plot = manhattan(outliers, p = "o_p_LO", chr = "CHR", bp = "BP", snp = "SNP", main= "Odd Lineage Longitude Outliers")
plot = manhattan(outliers, p = "o_p_LA", chr = "CHR", bp = "BP", snp = "SNP", main= "Odd Lineage Latitude Outliers")

#even pvalues: e_p_PT1	e_p_PT2	e_p_Le	e_p_LA	
# outliers: even_PT1	even_PT2	even_LO	even_LA	even_any

plot = manhattan(outliers, p = "LFMM_ASE_K2_PT1", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASE_K2_PT1", col = c("blue4", "orange3"))
plot = manhattan(outliers, p = "LFMM_ASE_K2_LO", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASE_K2_LO", col = c("blue4", "orange3"))
plot = manhattan(outliers, p = "LFMM_ASE_K2_LA", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASE_K2_LA", col = c("blue4", "orange3"))

#overlap
#pvalues:
# outliers: overlap_PT1	overlap_PT2	overlap_LO	overlap_LA

plot = manhattan(outliers, p = "LFMM_ASE_K2_LO", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASE_K2_LO", col = c("blue4", "orange3"))

#odd vs even
#pvalues:
# outliers: 
plot = manhattan(outliers, p = "LFMM_ASE_K2_LO", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASE_K2_LO", col = c("blue4", "orange3"))


###############################################
#Arlequin

#pvalues: o_p_Fst_arl	e_p_Fst_arl	
#ouliers: odd_0.01	even_0.01

plot = manhattan(outliers, p = "LFMM_ASE_K2_PT2", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASE_K2_PT2", col = c("blue4", "orange3"))
plot = manhattan(outliers, p = "LFMM_ASE_K2_PT1", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASE_K2_PT1", col = c("blue4", "orange3"))

#overlap
#pvalues:
# outliers: 
plot = manhattan(outliers, p = "LFMM_ASE_K2_LO", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASE_K2_LO", col = c("blue4", "orange3"))

#odd vs even
#pvalues:
# outliers: 

plot = manhattan(outliers, p = "LFMM_ASE_K2_LA", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASE_K2_LA", col = c("blue4", "orange3"))


#dev.off()
