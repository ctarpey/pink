### Population Genotype QC based on initial PCA
###   This code looks at the summary stats for the new pops compared to the old to 
###    better understand the spread of the populations in PCAs
###  Carolyn Tarpey | November 2017
### ---------------------------------------


library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)

#load genepop Fst output
New_Old_Fst<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Genepop/QC_NewPops_63758/Old_New_fst.txt")
colnames(New_Old_Fst)<- c("Locus","Fst")
New_Old_Fst<-New_Old_Fst[-1,]
head(New_Old_Fst)

New_Old_Fst_ordered<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Genepop/QC_NewPops_63758/Fst_sorted.txt")
colnames(New_Old_Fst_ordered)<- c("Locus","Fst","Rank")
New_Old_Fst_ordered<-New_Old_Fst_ordered[-1,]
head(New_Old_Fst_ordered)

#need to convert to integers
ggplot()+geom_histogram(data=New_Old_Fst_ordered[New_Old_Fst_ordered$Fst>=0.01],aes(x=New_Old_Fst_ordered$Fst),binwidth=0.01)+ggtitle("Histogram of Fst")
#plot genotype rate per locus relative to MAF for MAF 0.02 filtered loci
ggplot(data=New_Old_Fst_ordered[New_Old_Fst_ordered$Fst>=0.01],aes(x=New_Old_Fst_ordered$Rank, y=New_Old_Fst_ordered$Fst))
