### Manhattan plots of LFMM results
###    using this link: http://www.gettinggeneticsdone.com/2014/05/qqman-r-package-for-qq-and-manhattan-plots-for-gwas-results.html
### Carolyn Tarpey | March 2106
### ---------------------------------------

#install.packages("qqman")
library(qqman)

setwd('G:/Analysis/Mapping/Outliers_Manhattan_Plot')

#vignette("qqman")

Bayenv2results = read.table('Bayenv2_MAPoutliers.txt', header = TRUE)
as.data.frame(Bayenv2results)

#Bayenv2results$logPT1 = log10(Bayenv2results$PT1)
head(Bayenv2results)
table(Bayenv2results$CHR)

pdf("plots/Bayenv2results_Manhattan_plots2.pdf", width = 8, height = 5)
# plot = manhattan(Bayenv2results,logp = FALSE,, p = "LOG10PT1", chr = "CHR", bp = "BP", snp = "SNP", main= "LOG10PT1", col = c("blue4", "orange3"))
# plot = manhattan(Bayenv2results, p = "LOG10PT2", chr = "CHR", bp = "BP", snp = "SNP", main= "LOG10PT2", col = c("blue4", "orange3"))
# plot = manhattan(Bayenv2results,  p = "LOG10LONG", chr = "CHR", bp = "BP", snp = "SNP", main= "LOG10LONG", col = c("blue4", "orange3"))
# plot = manhattan(Bayenv2results, p = "LOG10LAT", chr = "CHR", bp = "BP", snp = "SNP", main= "LOG10LAT", col = c("blue4", "orange3"))

plot = manhattan(Bayenv2results, logp = FALSE, p = "PT1", chr = "CHR", bp = "BP", snp = "SNP", main= "log10 BF PT1", col = c("blue4", "orange3"))
plot = manhattan(Bayenv2results, logp = FALSE, p = "LONG", chr = "CHR", bp = "BP", snp = "SNP", main= "log10 BF LONG", col = c("blue4", "orange3"))
plot = manhattan(Bayenv2results, logp = FALSE, p = "LAT", chr = "CHR", bp = "BP", snp = "SNP", main= "log10 BF LAT", col = c("blue4", "orange3"))
plot = manhattan(Bayenv2results, logp = FALSE, p = "PT2", chr = "CHR", bp = "BP", snp = "SNP", main= "log10 BF PT2", col = c("blue4", "orange3"))


dev.off()

-log10(0.05)
