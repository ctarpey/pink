### Plotting a PCA for each pop with Adegenet
###   Data not completely filtered
### Carolyn Tarpey | November 2017
### ---------------------------------------


library("RColorBrewer")
library(ggplot2)
library(colorspace)
library(plyr)
library(colorRamps)
#install.packages("colorspace")
#install.packages("colorRamps")
###############DATA

pop_key <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/Summary_Stats/POPINFO.txt", header = TRUE, sep = '\t')

All_Genepop <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/PCA/filter_inds_pca.eigenvec")
head(ALL_eigenvec_table)
ALL_tt = merge(ALL_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')
ALL_geo <- ALL_tt[order(ALL_tt$Order_geo),] #sort by a geographical order, then odd then even 
head(ALL_tt)
