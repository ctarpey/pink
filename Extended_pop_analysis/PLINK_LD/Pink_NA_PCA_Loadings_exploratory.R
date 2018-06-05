###  North American PCA loadings Investigation
###   Data exploration to see if there are any loci 
###     that are differentiating the populations. 
###
### Carolyn Tarpey | June 2018
### ---------------------------------------

library(RColorBrewer)
library(ggplot2)
library(colorspace)
library(plyr)
library(colorRamps)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)
library(reshape2)
library(plotly)
library(gridExtra)
library(scales) 
library(grid)


###############Import the pop info for the populations (IT IS UPDATED!!)
#pop_key <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/POPINFO.txt", header = TRUE, sep = '\t')

###############Import the eigenvalues that were run in PLINK and merge them with the Pop info
NA_eigen_val_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAsofNA/na_all_pca.eigenvec.var")
NA_eigen_val_table$V1 <- NULL
names(NA_eigen_val_table)<- c("Locus", "Dim1","Dim2","Dim3","Dim4","Dim5","Dim6","Dim7","Dim8","Dim9","Dim10","Dim11","Dim12","Dim13","Dim14","Dim15","Dim16","Dim17","Dim18","Dim19","Dim20")
head(NA_eigen_val_table)

EVEN_NA_eigen_val_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAsofNA/even_na_pca.eigenvec.var")
EVEN_NA_eigen_val_table$V1 <- NULL
names(EVEN_NA_eigen_val_table)<- c("Locus", "Dim1","Dim2","Dim3","Dim4","Dim5","Dim6","Dim7","Dim8","Dim9","Dim10","Dim11","Dim12","Dim13","Dim14","Dim15","Dim16","Dim17","Dim18","Dim19","Dim20")
head(EVEN_NA_eigen_val_table)

ODD_NA_eigen_val_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAsofNA/odd_na_pca.eigenvec.var")
ODD_NA_eigen_val_table$V1 <- NULL
names(ODD_NA_eigen_val_table)<- c("Locus", "Dim1","Dim2","Dim3","Dim4","Dim5","Dim6","Dim7","Dim8","Dim9","Dim10","Dim11","Dim12","Dim13","Dim14","Dim15","Dim16","Dim17","Dim18","Dim19","Dim20")
head(ODD_NA_eigen_val_table)


###plot the eigenvalues to see what the ranges are for each of the top 3 dimmensions foreach grouping. 

ggplot(NA_eigen_val_table) + geom_histogram(aes(x=NA_eigen_val_table$Dim1), binwidth = 0.1) +theme_bw() +
  xlab("Dimension 1 Loadings") + ylab("Count of loadings per bin") + ggtitle("Histogram of the Loadings for PCA North American populations Dimesion 1") 

ggplot(NA_eigen_val_table) + geom_histogram(aes(x=NA_eigen_val_table$Dim2), binwidth = 0.1) +theme_bw() +
  xlab("Dimension 2 Loadings") + ylab("Count of loadings per bin") + ggtitle("Histogram of the Loadings for PCA North American populations Dimesion 2") 

ggplot(NA_eigen_val_table) + geom_histogram(aes(x=NA_eigen_val_table$Dim3), binwidth = 0.1) +theme_bw() +
  xlab("Dimension 3 Loadings") + ylab("Count of loadings per bin") + ggtitle("Histogram of the Loadings for PCA North American populations Dimesion 3") 



