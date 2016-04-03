### Manhattan plots of Outliers from LFMM and Arlequin
###    using this link: http://www.gettinggeneticsdone.com/2014/05/qqman-r-package-for-qq-and-manhattan-plots-for-gwas-results.html
### Carolyn Tarpey | April 2106
### ---------------------------------------

#install.packages("qqman")
library(qqman)
vignette("qqman")

setwd('G:/Analysis/Mapping/Outliers_Manhattan_Plot/Arlequin_LFMM')


outliers = read.table('Outliers_Map.txt', header = TRUE)
as.data.frame(outliers)

head(outliers)
table(outliers$CHR)


#SNP CHR     BP Paralog    o_p_PT1     o_p_PT2      o_p_LO     o_p_LA    e_p_PT1    e_p_PT2      e_p_Le     e_p_LA 
#o_p_Fst_arl e_p_Fst_arl odd_PT1 odd_PT2 odd_LO odd_LA odd_any even_PT1 even_PT2 even_LO even_LA even_any 
#odd_0.01 even_0.01 overlap_0.01 overlap_PT1 overlap_PT2 overlap_LO   overlap_LA

#Position is defined as position on chromosome. Need to add largest position from last marker on previous chr so that markers are ordered correctly along x-axis

chrNum=26

for (i in 1:chrNum){ ndx <- which(outliers[, 2]==i)
lstMrk <- max(outliers[ndx, 3])
if (i < chrNum) ndx2 <- which(outliers[, 2]==i+1)
if (i < chrNum) outliers[ndx2, 3] <- outliers[ndx2, 3] + lstMrk
}

bpMidVec <- vector(length=chrNum)

for (i in 1:chrNum){ndx <- which(outliers[, 2]==i)
posSub <- outliers[ndx, 3]
bpMidVec[i] <- ((max(posSub) - min(posSub))/2) + min(posSub)
}


p <- ggplot(outliers) + geom_point(aes(x=BP, y=o_p_PT1, size=3.5, colour=as.factor(CHR)), alpha=1/3) + 
  scale_color_manual(values=rep(c('black', 'dark green'), 6)) + theme_bw(base_size=15) + theme(legend.position='none') + 
  scale_x_continuous(labels=as.character(1:chrNum), breaks=bpMidVec) +  
  ggtitle('Odd Lineage Outliers PT1') + xlab('Linkage Group ') + ylab('-log10(P)')


#geom_hline(y=4.08, linetype=1, col='red', lwd=1.5) +







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
