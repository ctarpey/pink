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
pop_key <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/POPINFO_LS.txt", header = TRUE, sep = '\t')

###############Import the eigenvalues that were run in PLINK and merge them with the Pop info
NA_Nome_eigen_vec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/na_nome_pops_pca.eigenvec.var")
NA_Nome_eigen_vec_table$V1 <- NULL
names(NA_Nome_eigen_vec_table)<- c("Locus", "Dim1","Dim2","Dim3","Dim4","Dim5","Dim6","Dim7","Dim8","Dim9","Dim10","Dim11","Dim12","Dim13","Dim14","Dim15","Dim16","Dim17","Dim18","Dim19","Dim20")
head(NA_Nome_eigen_vec_table)

Even_NA_Nome_eigen_vec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/e_na_nome_pca.eigenvec.var")
Even_NA_Nome_eigen_vec_table$V1 <- NULL
names(Even_NA_Nome_eigen_vec_table)<- c("Locus", "Dim1","Dim2","Dim3","Dim4","Dim5","Dim6","Dim7","Dim8","Dim9","Dim10","Dim11","Dim12","Dim13","Dim14","Dim15","Dim16","Dim17","Dim18","Dim19","Dim20")
head(Even_NA_Nome_eigen_vec_table)

Odd_NA_Nome_eigen_vec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/o_na_nome_pca.eigenvec.var")
Odd_NA_Nome_eigen_vec_table$V1 <- NULL
names(Odd_NA_Nome_eigen_vec_table)<- c("Locus", "Dim1","Dim2","Dim3","Dim4","Dim5","Dim6","Dim7","Dim8","Dim9","Dim10","Dim11","Dim12","Dim13","Dim14","Dim15","Dim16","Dim17","Dim18","Dim19","Dim20")
head(Odd_NA_Nome_eigen_vec_table)

####Plot the eigenvalues  to see which dimmensions appear to be the  most influential for each group

NA_Nome_eigenval <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/na_nome_pops_pca.eigenval")
NA_Nome_eigenval$Dim <- c(1:20)

Even_NA_Nome_eigenval <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/e_na_nome_pca.eigenval")
Even_NA_Nome_eigenval$Dim <- c(1:20)

Odd_NA_Nome_eigenval <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/o_na_nome_pca.eigenval")
Odd_NA_Nome_eigenval$Dim <- c(1:20)

##plot 
ggplot() + geom_bar(data=NA_Nome_eigenval, aes(x=Dim, y=V1),stat= "identity") + theme_bw() +
  xlab("Dimension") + ylab("Eigenvalue") + ggtitle("Eigenvalues the PCA of North American w/ Nome") 

ggplot() + geom_bar(data=Even_NA_Nome_eigenval, aes(x=Dim, y=V1),stat= "identity") + theme_bw() +
  xlab("Dimension") + ylab("Eigenvalue") + ggtitle("Eigenvalues the PCA of Even lineage North American w/ Nome") 

ggplot() + geom_bar(data=Odd_NA_Nome_eigenval, aes(x=Dim, y=V1),stat= "identity") + theme_bw() +
  xlab("Dimension") + ylab("Eigenvalue") + ggtitle("Eigenvalues the PCA of Odd lineage North American w/ Nome") 


###plot the eigenvalues to see what the ranges are for each of the top 3 dimmensions foreach grouping. 
###################################

ggplot(NA_Nome_eigen_vec_table) + geom_histogram(aes(x=NA_Nome_eigen_vec_table$Dim1), binwidth = 0.1) +theme_bw() +
  xlab("Dimension 1 Loadings") + ylab("Count of loadings per bin") + ggtitle("Histogram of the Loadings for PCA North American w/Nome Dimesion 1") 

ggplot(NA_Nome_eigen_vec_table) + geom_histogram(aes(x=NA_Nome_eigen_vec_table$Dim2), binwidth = 0.1) +theme_bw() +
  xlab("Dimension 2 Loadings") + ylab("Count of loadings per bin") + ggtitle("Histogram of the Loadings for PCA North American w/Nome Dimesion 2") 


ggplot(Even_NA_Nome_eigen_vec_table) + geom_histogram(aes(x=Even_NA_Nome_eigen_vec_table$Dim2), binwidth = 0.1) +theme_bw() +
  xlab("Dimension 1 Loadings") + ylab("Count of loadings per bin") + ggtitle("Histogram of the Loadings for PCA of Even North American w/Nome Dimesion 1") 

ggplot(Even_NA_Nome_eigen_vec_table) + geom_histogram(aes(x=Even_NA_Nome_eigen_vec_table$Dim2), binwidth = 0.1) +theme_bw() +
  xlab("Dimension 2 Loadings") + ylab("Count of loadings per bin") + ggtitle("Histogram of the Loadings for PCA of Even North American w/Nome Dimesion 2") 

ggplot(Even_NA_Nome_eigen_vec_table) + geom_histogram(aes(x=Even_NA_Nome_eigen_vec_table$Dim2), binwidth = 0.1) +theme_bw() +
  xlab("Dimension 3 Loadings") + ylab("Count of loadings per bin") + ggtitle("Histogram of the Loadings for PCA of Even North American  w/Nome Dimesion 3") 


ggplot(Odd_NA_Nome_eigen_vec_table) + geom_histogram(aes(x=Odd_NA_Nome_eigen_vec_table$Dim3), binwidth = 0.1) +theme_bw() +
  xlab("Dimension 1 Loadings") + ylab("Count of loadings per bin") + ggtitle("Histogram of the Loadings for PCA of Odd North American  w/Nome  Dimesion 1") 

ggplot(Odd_NA_Nome_eigen_vec_table) + geom_histogram(aes(x=Odd_NA_Nome_eigen_vec_table$Dim3), binwidth = 0.1) +theme_bw() +
  xlab("Dimension 2 Loadings") + ylab("Count of loadings per bin") + ggtitle("Histogram of the Loadings for PCA of Odd North American  w/Nome  Dimesion 2") 

ggplot(Odd_NA_Nome_eigen_vec_table) + geom_histogram(aes(x=Odd_NA_Nome_eigen_vec_table$Dim3), binwidth = 0.1) +theme_bw() +
  xlab("Dimension 3 Loadings") + ylab("Count of loadings per bin") + ggtitle("Histogram of the Loadings for PCA of Odd North American  w/Nome  Dimesion 3") 


