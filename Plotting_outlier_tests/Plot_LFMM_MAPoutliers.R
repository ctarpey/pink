### Manhattan plots of LFMM results
###    using this link: http://www.gettinggeneticsdone.com/2014/05/qqman-r-package-for-qq-and-manhattan-plots-for-gwas-results.html
### Carolyn Tarpey | March 2106
### ---------------------------------------

#install.packages("qqman")
library(qqman)


setwd('G:/Analysis/Mapping/Outliers_Manhattan_Plot')


vignette("qqman")

LFMMresults = read.table('LFMM_MAPoutliers.txt', header = TRUE)
as.data.frame(LFMMresults)

head(LFMMresults)
table(LFMMresults$CHR)

pdf("plots/LFMMresults_Manhattan_plots.pdf", width = 8, height = 5)
plot = manhattan(LFMMresults, p = "LFMM_ASE_K2_PT2", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASE_K2_PT2", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_ASE_K2_PT1", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASE_K2_PT1", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_ASE_K2_LO", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASE_K2_LO", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_ASE_K2_LA", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASE_K2_LA", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_ASO_K2_PT2", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASO_K2_PT2", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_ASO_K2_PT1", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASO_K2_PT1", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_ASO_K2_LO", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASO_K2_LO", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_ASO_K2_LA", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASO_K2_LA", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_NAO_K3_PT2", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_NAO_K3_PT2", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_NAO_K3_PT1", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_NAO_K3_PT1", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_NAO_K3_LO", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_NAO_K3_LO", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_NAO_K3_LA", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_NAO_K3_LA", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_NAE_K2_PT2", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_NAE_K2_PT2", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_NAE_K2_PT1", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_NAE_K2_PT1", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_NAE_K2_LO", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_NAE_K2_LO", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_NAE_K2_LA", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_NAE_K2_LA", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_ALL_K4_PT2", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ALL_K4_PT2", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_ALL_K4_PT1", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ALL_K4_PT1", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_ALL_K4_LO", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ALL_K4_LO", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_ALL_K4_LA", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ALL_K4_LA", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_ALL_K7_PT2", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ALL_K7_PT2", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_ALL_K7_PT1", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ALL_K7_PT1", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_ALL_K7_LO", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ALL_K7_LO", col = c("blue4", "orange3"))
plot = manhattan(LFMMresults, p = "LFMM_ALL_K7_LA", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ALL_K7_LA", col = c("blue4", "orange3"))


dev.off()

-log10(0.05)
