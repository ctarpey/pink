### Plotting the LD r2 from PLINK 
###     REGIONAL GROUPINGS 
###    for the full dataset of One SNP per Tag- 23759 loci. 
###    THE CODE FOR THE PLOTS OF THE ACTUAL LD /CHROM IS IN _2 
###  Carolyn Tarpey | May 2018
### ----------------------------------

#The full 13 population groups takes up too much memory, so use the workspace Pink_PlottingLD_mainHierarchies.R.RDATA for these 7 groups
# Load the rest of the groups in the _mainHierarchies R scripts, and the NA_subset for the breakout of the North American populations. 

#This R code was originally written by Garrett McKinney to plot the output of PLINK LD r2 for each chromosome. 
#It requires the LD output from PLINK. We used the alignment of the pinks to Chinook to get the Chromosome 
#and position assignments here, so there are 34 Chromosomes, and chromosome 35 are those that aligned to the unassigned 

#TO MAKE THE PLOTS OF THE LD PER CHROMOSOME, USE THIS SCRIPT TO IMPORT THE DATA AND THE Pink_PlottingLD_2.R TO MAKE THE ACTUAL PLOTS

library(ggplot2)
library(vcfR)
library(stringr)
library(plyr)
library(dplyr)
library(gdata)
library(reshape2)
library(plotly)
library(gridExtra)
library(scales) 
library(grid)


#load PLINK LD output files

#### Even Asia
#load PLINK LD output files
LD_even_A_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr1.ld",sep="")
LD_even_A_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr2.ld",sep="")
LD_even_A_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr3.ld",sep="")
LD_even_A_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr4.ld",sep="")
LD_even_A_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr5.ld",sep="")
LD_even_A_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr6.ld",sep="")
LD_even_A_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr7.ld",sep="")
LD_even_A_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr8.ld",sep="")
LD_even_A_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr9.ld",sep="")
LD_even_A_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr10.ld",sep="")
LD_even_A_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr11.ld",sep="")
LD_even_A_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr12.ld",sep="")
LD_even_A_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr13.ld",sep="")
LD_even_A_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr14.ld",sep="")
LD_even_A_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr15.ld",sep="")
LD_even_A_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr16.ld",sep="")
LD_even_A_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr17.ld",sep="")
LD_even_A_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr18.ld",sep="")
LD_even_A_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr19.ld",sep="")
LD_even_A_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr20.ld",sep="")
LD_even_A_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr21.ld",sep="")
LD_even_A_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr22.ld",sep="")
LD_even_A_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr23.ld",sep="")
LD_even_A_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr24.ld",sep="")
LD_even_A_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr25.ld",sep="")
LD_even_A_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr26.ld",sep="")
LD_even_A_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr27.ld",sep="")
LD_even_A_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr28.ld",sep="")
LD_even_A_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr29.ld",sep="")
LD_even_A_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr30.ld",sep="")
LD_even_A_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr31.ld",sep="")
LD_even_A_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr32.ld",sep="")
LD_even_A_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr33.ld",sep="")
LD_even_A_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr34.ld",sep="")
LD_even_A_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EA/PLINK_out_EA_chr35.ld",sep="")


LD_data_even_A <- rbind( LD_even_A_chr1,  LD_even_A_chr2, LD_even_A_chr3, LD_even_A_chr4, LD_even_A_chr5, LD_even_A_chr6,  LD_even_A_chr7, LD_even_A_chr8, LD_even_A_chr9,
                         LD_even_A_chr10, LD_even_A_chr11,  LD_even_A_chr12,  LD_even_A_chr13, LD_even_A_chr14,  LD_even_A_chr15,  LD_even_A_chr16,  LD_even_A_chr17,  LD_even_A_chr18,  LD_even_A_chr19,
                         LD_even_A_chr20,  LD_even_A_chr21,  LD_even_A_chr22, LD_even_A_chr23,  LD_even_A_chr24, LD_even_A_chr25,  LD_even_A_chr26,  LD_even_A_chr27,  LD_even_A_chr28,  LD_even_A_chr29,  LD_even_A_chr30,  LD_even_A_chr31,  LD_even_A_chr32,  LD_even_A_chr33, LD_even_A_chr34, LD_even_A_chr0)
dim(LD_data_even_A)
rm( LD_even_A_chr1,  LD_even_A_chr2, LD_even_A_chr3, LD_even_A_chr4, LD_even_A_chr5, LD_even_A_chr6,  LD_even_A_chr7, LD_even_A_chr8, LD_even_A_chr9,
    LD_even_A_chr10, LD_even_A_chr11,  LD_even_A_chr12,  LD_even_A_chr13, LD_even_A_chr14,  LD_even_A_chr15,  LD_even_A_chr16,  LD_even_A_chr17,  LD_even_A_chr18,  LD_even_A_chr19,
    LD_even_A_chr20,  LD_even_A_chr21,  LD_even_A_chr22, LD_even_A_chr23,  LD_even_A_chr24, LD_even_A_chr25,  LD_even_A_chr26,  LD_even_A_chr27,  LD_even_A_chr28,  LD_even_A_chr29,  LD_even_A_chr30,  LD_even_A_chr31,  LD_even_A_chr32,  LD_even_A_chr33, LD_even_A_chr34, LD_even_A_chr0)



#### Even North America
#load PLINK LD output files
LD_even_NA_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr1.ld",sep="")
LD_even_NA_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr2.ld",sep="")
LD_even_NA_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr3.ld",sep="")
LD_even_NA_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr4.ld",sep="")
LD_even_NA_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr5.ld",sep="")
LD_even_NA_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr6.ld",sep="")
LD_even_NA_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr7.ld",sep="")
LD_even_NA_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr8.ld",sep="")
LD_even_NA_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr9.ld",sep="")
LD_even_NA_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr10.ld",sep="")
LD_even_NA_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr11.ld",sep="")
LD_even_NA_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr12.ld",sep="")
LD_even_NA_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr13.ld",sep="")
LD_even_NA_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr14.ld",sep="")
LD_even_NA_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr15.ld",sep="")
LD_even_NA_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr16.ld",sep="")
LD_even_NA_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr17.ld",sep="")
LD_even_NA_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr18.ld",sep="")
LD_even_NA_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr19.ld",sep="")
LD_even_NA_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr20.ld",sep="")
LD_even_NA_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr21.ld",sep="")
LD_even_NA_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr22.ld",sep="")
LD_even_NA_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr23.ld",sep="")
LD_even_NA_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr24.ld",sep="")
LD_even_NA_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr25.ld",sep="")
LD_even_NA_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr26.ld",sep="")
LD_even_NA_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr27.ld",sep="")
LD_even_NA_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr28.ld",sep="")
LD_even_NA_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr29.ld",sep="")
LD_even_NA_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr30.ld",sep="")
LD_even_NA_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr31.ld",sep="")
LD_even_NA_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr32.ld",sep="")
LD_even_NA_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr33.ld",sep="")
LD_even_NA_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr34.ld",sep="")
LD_even_NA_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA/PLINK_out_ENA_chr35.ld",sep="")



LD_data_even_NA <- rbind( LD_even_NA_chr1,  LD_even_NA_chr2, LD_even_NA_chr3, LD_even_NA_chr4, LD_even_NA_chr5, LD_even_NA_chr6,  LD_even_NA_chr7, LD_even_NA_chr8, LD_even_NA_chr9,
                          LD_even_NA_chr10, LD_even_NA_chr11,  LD_even_NA_chr12,  LD_even_NA_chr13, LD_even_NA_chr14,  LD_even_NA_chr15,  LD_even_NA_chr16,  LD_even_NA_chr17,  LD_even_NA_chr18,  LD_even_NA_chr19,
                          LD_even_NA_chr20,  LD_even_NA_chr21,  LD_even_NA_chr22, LD_even_NA_chr23,  LD_even_NA_chr24, LD_even_NA_chr25,  LD_even_NA_chr26,  LD_even_NA_chr27,  LD_even_NA_chr28,  LD_even_NA_chr29,  LD_even_NA_chr30,  LD_even_NA_chr31,  LD_even_NA_chr32,  LD_even_NA_chr33, LD_even_NA_chr34, LD_even_NA_chr0)
dim(LD_data_even_NA)
rm( LD_even_NA_chr1,  LD_even_NA_chr2, LD_even_NA_chr3, LD_even_NA_chr4, LD_even_NA_chr5, LD_even_NA_chr6,  LD_even_NA_chr7, LD_even_NA_chr8, LD_even_NA_chr9,
    LD_even_NA_chr10, LD_even_NA_chr11,  LD_even_NA_chr12,  LD_even_NA_chr13, LD_even_NA_chr14,  LD_even_NA_chr15,  LD_even_NA_chr16,  LD_even_NA_chr17,  LD_even_NA_chr18,  LD_even_NA_chr19,
    LD_even_NA_chr20,  LD_even_NA_chr21,  LD_even_NA_chr22, LD_even_NA_chr23,  LD_even_NA_chr24, LD_even_NA_chr25,  LD_even_NA_chr26,  LD_even_NA_chr27,  LD_even_NA_chr28,  LD_even_NA_chr29,  LD_even_NA_chr30,  LD_even_NA_chr31,  LD_even_NA_chr32,  LD_even_NA_chr33, LD_even_NA_chr34, LD_even_NA_chr0)



#### odd Asia
#load PLINK LD output files
LD_odd_a_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr1.ld",sep="")
LD_odd_a_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr2.ld",sep="")
LD_odd_a_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr3.ld",sep="")
LD_odd_a_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr4.ld",sep="")
LD_odd_a_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr5.ld",sep="")
LD_odd_a_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr6.ld",sep="")
LD_odd_a_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr7.ld",sep="")
LD_odd_a_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr8.ld",sep="")
LD_odd_a_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr9.ld",sep="")
LD_odd_a_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr10.ld",sep="")
LD_odd_a_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr11.ld",sep="")
LD_odd_a_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr12.ld",sep="")
LD_odd_a_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr13.ld",sep="")
LD_odd_a_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr14.ld",sep="")
LD_odd_a_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr15.ld",sep="")
LD_odd_a_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr16.ld",sep="")
LD_odd_a_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr17.ld",sep="")
LD_odd_a_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr18.ld",sep="")
LD_odd_a_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr19.ld",sep="")
LD_odd_a_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr20.ld",sep="")
LD_odd_a_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr21.ld",sep="")
LD_odd_a_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr22.ld",sep="")
LD_odd_a_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr23.ld",sep="")
LD_odd_a_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr24.ld",sep="")
LD_odd_a_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr25.ld",sep="")
LD_odd_a_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr26.ld",sep="")
LD_odd_a_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr27.ld",sep="")
LD_odd_a_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr28.ld",sep="")
LD_odd_a_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr29.ld",sep="")
LD_odd_a_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr30.ld",sep="")
LD_odd_a_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr31.ld",sep="")
LD_odd_a_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr32.ld",sep="")
LD_odd_a_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr33.ld",sep="")
LD_odd_a_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr34.ld",sep="")
LD_odd_a_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_OA/PLINK_out_OA_chr35.ld",sep="")



LD_data_odd_a <- rbind( LD_odd_a_chr1,  LD_odd_a_chr2, LD_odd_a_chr3, LD_odd_a_chr4, LD_odd_a_chr5, LD_odd_a_chr6,  LD_odd_a_chr7, LD_odd_a_chr8, LD_odd_a_chr9,
                        LD_odd_a_chr10, LD_odd_a_chr11,  LD_odd_a_chr12,  LD_odd_a_chr13, LD_odd_a_chr14,  LD_odd_a_chr15,  LD_odd_a_chr16,  LD_odd_a_chr17,  LD_odd_a_chr18,  LD_odd_a_chr19,
                        LD_odd_a_chr20,  LD_odd_a_chr21,  LD_odd_a_chr22, LD_odd_a_chr23,  LD_odd_a_chr24, LD_odd_a_chr25,  LD_odd_a_chr26,  LD_odd_a_chr27,  LD_odd_a_chr28,  LD_odd_a_chr29,  LD_odd_a_chr30,  LD_odd_a_chr31,  LD_odd_a_chr32,  LD_odd_a_chr33, LD_odd_a_chr34, LD_odd_a_chr0)
dim(LD_data_odd_a)
rm( LD_odd_a_chr1,  LD_odd_a_chr2, LD_odd_a_chr3, LD_odd_a_chr4, LD_odd_a_chr5, LD_odd_a_chr6,  LD_odd_a_chr7, LD_odd_a_chr8, LD_odd_a_chr9,
    LD_odd_a_chr10, LD_odd_a_chr11,  LD_odd_a_chr12,  LD_odd_a_chr13, LD_odd_a_chr14,  LD_odd_a_chr15,  LD_odd_a_chr16,  LD_odd_a_chr17,  LD_odd_a_chr18,  LD_odd_a_chr19,
    LD_odd_a_chr20,  LD_odd_a_chr21,  LD_odd_a_chr22, LD_odd_a_chr23,  LD_odd_a_chr24, LD_odd_a_chr25,  LD_odd_a_chr26,  LD_odd_a_chr27,  LD_odd_a_chr28,  LD_odd_a_chr29,  LD_odd_a_chr30,  LD_odd_a_chr31,  LD_odd_a_chr32,  LD_odd_a_chr33, LD_odd_a_chr34, LD_odd_a_chr0)


#### Odd North America
#load PLINK LD output files
LD_odd_na_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr1.ld",sep="")
LD_odd_na_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr2.ld",sep="")
LD_odd_na_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr3.ld",sep="")
LD_odd_na_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr4.ld",sep="")
LD_odd_na_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr5.ld",sep="")
LD_odd_na_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr6.ld",sep="")
LD_odd_na_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr7.ld",sep="")
LD_odd_na_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr8.ld",sep="")
LD_odd_na_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr9.ld",sep="")
LD_odd_na_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr10.ld",sep="")
LD_odd_na_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr11.ld",sep="")
LD_odd_na_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr12.ld",sep="")
LD_odd_na_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr13.ld",sep="")
LD_odd_na_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr14.ld",sep="")
LD_odd_na_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr15.ld",sep="")
LD_odd_na_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr16.ld",sep="")
LD_odd_na_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr17.ld",sep="")
LD_odd_na_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr18.ld",sep="")
LD_odd_na_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr19.ld",sep="")
LD_odd_na_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr20.ld",sep="")
LD_odd_na_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr21.ld",sep="")
LD_odd_na_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr22.ld",sep="")
LD_odd_na_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr23.ld",sep="")
LD_odd_na_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr24.ld",sep="")
LD_odd_na_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr25.ld",sep="")
LD_odd_na_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr26.ld",sep="")
LD_odd_na_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr27.ld",sep="")
LD_odd_na_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr28.ld",sep="")
LD_odd_na_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr29.ld",sep="")
LD_odd_na_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr30.ld",sep="")
LD_odd_na_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr31.ld",sep="")
LD_odd_na_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr32.ld",sep="")
LD_odd_na_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr33.ld",sep="")
LD_odd_na_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr34.ld",sep="")
LD_odd_na_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA/PLINK_out_ONA_chr35.ld",sep="")


LD_data_odd_na <- rbind( LD_odd_na_chr1,  LD_odd_na_chr2, LD_odd_na_chr3, LD_odd_na_chr4, LD_odd_na_chr5, LD_odd_na_chr6,  LD_odd_na_chr7, LD_odd_na_chr8, LD_odd_na_chr9,
                         LD_odd_na_chr10, LD_odd_na_chr11,  LD_odd_na_chr12,  LD_odd_na_chr13, LD_odd_na_chr14,  LD_odd_na_chr15,  LD_odd_na_chr16,  LD_odd_na_chr17,  LD_odd_na_chr18,  LD_odd_na_chr19,
                         LD_odd_na_chr20,  LD_odd_na_chr21,  LD_odd_na_chr22, LD_odd_na_chr23,  LD_odd_na_chr24, LD_odd_na_chr25,  LD_odd_na_chr26,  LD_odd_na_chr27,  LD_odd_na_chr28,  LD_odd_na_chr29,  LD_odd_na_chr30,  LD_odd_na_chr31,  LD_odd_na_chr32,  LD_odd_na_chr33, LD_odd_na_chr34, LD_odd_na_chr0)
dim(LD_data_odd_na)
rm( LD_odd_na_chr1,  LD_odd_na_chr2, LD_odd_na_chr3, LD_odd_na_chr4, LD_odd_na_chr5, LD_odd_na_chr6,  LD_odd_na_chr7, LD_odd_na_chr8, LD_odd_na_chr9,
    LD_odd_na_chr10, LD_odd_na_chr11,  LD_odd_na_chr12,  LD_odd_na_chr13, LD_odd_na_chr14,  LD_odd_na_chr15,  LD_odd_na_chr16,  LD_odd_na_chr17,  LD_odd_na_chr18,  LD_odd_na_chr19,
    LD_odd_na_chr20,  LD_odd_na_chr21,  LD_odd_na_chr22, LD_odd_na_chr23,  LD_odd_na_chr24, LD_odd_na_chr25,  LD_odd_na_chr26,  LD_odd_na_chr27,  LD_odd_na_chr28,  LD_odd_na_chr29,  LD_odd_na_chr30,  LD_odd_na_chr31,  LD_odd_na_chr32,  LD_odd_na_chr33, LD_odd_na_chr34, LD_odd_na_chr0)


########################THESE HAVE NO SUSITNA POPULATIONS INCLUDED :

#### Odd North America NO SUSITNA
#load PLINK LD output files
LD_odd_na_ns_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr1.ld",sep="")
LD_odd_na_ns_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr2.ld",sep="")
LD_odd_na_ns_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr3.ld",sep="")
LD_odd_na_ns_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr4.ld",sep="")
LD_odd_na_ns_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr5.ld",sep="")
LD_odd_na_ns_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr6.ld",sep="")
LD_odd_na_ns_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr7.ld",sep="")
LD_odd_na_ns_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr8.ld",sep="")
LD_odd_na_ns_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr9.ld",sep="")
LD_odd_na_ns_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr10.ld",sep="")
LD_odd_na_ns_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr11.ld",sep="")
LD_odd_na_ns_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr12.ld",sep="")
LD_odd_na_ns_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr13.ld",sep="")
LD_odd_na_ns_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr14.ld",sep="")
LD_odd_na_ns_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr15.ld",sep="")
LD_odd_na_ns_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr16.ld",sep="")
LD_odd_na_ns_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr17.ld",sep="")
LD_odd_na_ns_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr18.ld",sep="")
LD_odd_na_ns_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr19.ld",sep="")
LD_odd_na_ns_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr20.ld",sep="")
LD_odd_na_ns_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr21.ld",sep="")
LD_odd_na_ns_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr22.ld",sep="")
LD_odd_na_ns_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr23.ld",sep="")
LD_odd_na_ns_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr24.ld",sep="")
LD_odd_na_ns_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr25.ld",sep="")
LD_odd_na_ns_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr26.ld",sep="")
LD_odd_na_ns_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr27.ld",sep="")
LD_odd_na_ns_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr28.ld",sep="")
LD_odd_na_ns_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr29.ld",sep="")
LD_odd_na_ns_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr30.ld",sep="")
LD_odd_na_ns_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr31.ld",sep="")
LD_odd_na_ns_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr32.ld",sep="")
LD_odd_na_ns_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr33.ld",sep="")
LD_odd_na_ns_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr34.ld",sep="")
LD_odd_na_ns_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ONA_NS/PLINK_out_ONA_NS_chr35.ld",sep="")


LD_data_odd_na_ns <- rbind( LD_odd_na_ns_chr1,  LD_odd_na_ns_chr2, LD_odd_na_ns_chr3, LD_odd_na_ns_chr4, LD_odd_na_ns_chr5, LD_odd_na_ns_chr6,  LD_odd_na_ns_chr7, LD_odd_na_ns_chr8, LD_odd_na_ns_chr9,
                            LD_odd_na_ns_chr10, LD_odd_na_ns_chr11,  LD_odd_na_ns_chr12,  LD_odd_na_ns_chr13, LD_odd_na_ns_chr14,  LD_odd_na_ns_chr15,  LD_odd_na_ns_chr16,  LD_odd_na_ns_chr17,  LD_odd_na_ns_chr18,  LD_odd_na_ns_chr19,
                            LD_odd_na_ns_chr20,  LD_odd_na_ns_chr21,  LD_odd_na_ns_chr22, LD_odd_na_ns_chr23,  LD_odd_na_ns_chr24, LD_odd_na_ns_chr25,  LD_odd_na_ns_chr26,  LD_odd_na_ns_chr27,  LD_odd_na_ns_chr28,  LD_odd_na_ns_chr29,  LD_odd_na_ns_chr30,  LD_odd_na_ns_chr31,  LD_odd_na_ns_chr32,  LD_odd_na_ns_chr33, LD_odd_na_ns_chr34, LD_odd_na_ns_chr0)
dim(LD_data_odd_na_ns)
rm( LD_odd_na_ns_chr1,  LD_odd_na_ns_chr2, LD_odd_na_ns_chr3, LD_odd_na_ns_chr4, LD_odd_na_ns_chr5, LD_odd_na_ns_chr6,  LD_odd_na_ns_chr7, LD_odd_na_ns_chr8, LD_odd_na_ns_chr9,
    LD_odd_na_ns_chr10, LD_odd_na_ns_chr11,  LD_odd_na_ns_chr12,  LD_odd_na_ns_chr13, LD_odd_na_ns_chr14,  LD_odd_na_ns_chr15,  LD_odd_na_ns_chr16,  LD_odd_na_ns_chr17,  LD_odd_na_ns_chr18,  LD_odd_na_ns_chr19,
    LD_odd_na_ns_chr20,  LD_odd_na_ns_chr21,  LD_odd_na_ns_chr22, LD_odd_na_ns_chr23,  LD_odd_na_ns_chr24, LD_odd_na_ns_chr25,  LD_odd_na_ns_chr26,  LD_odd_na_ns_chr27,  LD_odd_na_ns_chr28,  LD_odd_na_ns_chr29,  LD_odd_na_ns_chr30,  LD_odd_na_ns_chr31,  LD_odd_na_ns_chr32,  LD_odd_na_ns_chr33, LD_odd_na_ns_chr34, LD_odd_na_ns_chr0)


#### EVEN North America NO SUSITNA
#load PLINK LD output files
LD_even_na_ns_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr1.ld",sep="")
LD_even_na_ns_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr2.ld",sep="")
LD_even_na_ns_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr3.ld",sep="")
LD_even_na_ns_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr4.ld",sep="")
LD_even_na_ns_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr5.ld",sep="")
LD_even_na_ns_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr6.ld",sep="")
LD_even_na_ns_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr7.ld",sep="")
LD_even_na_ns_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr8.ld",sep="")
LD_even_na_ns_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr9.ld",sep="")
LD_even_na_ns_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr10.ld",sep="")
LD_even_na_ns_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr11.ld",sep="")
LD_even_na_ns_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr12.ld",sep="")
LD_even_na_ns_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr13.ld",sep="")
LD_even_na_ns_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr14.ld",sep="")
LD_even_na_ns_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr15.ld",sep="")
LD_even_na_ns_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr16.ld",sep="")
LD_even_na_ns_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr17.ld",sep="")
LD_even_na_ns_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr18.ld",sep="")
LD_even_na_ns_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr19.ld",sep="")
LD_even_na_ns_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr20.ld",sep="")
LD_even_na_ns_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr21.ld",sep="")
LD_even_na_ns_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr22.ld",sep="")
LD_even_na_ns_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr23.ld",sep="")
LD_even_na_ns_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr24.ld",sep="")
LD_even_na_ns_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr25.ld",sep="")
LD_even_na_ns_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr26.ld",sep="")
LD_even_na_ns_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr27.ld",sep="")
LD_even_na_ns_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr28.ld",sep="")
LD_even_na_ns_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr29.ld",sep="")
LD_even_na_ns_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr30.ld",sep="")
LD_even_na_ns_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr31.ld",sep="")
LD_even_na_ns_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr32.ld",sep="")
LD_even_na_ns_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr33.ld",sep="")
LD_even_na_ns_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr34.ld",sep="")
LD_even_na_ns_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ENA_NS/PLINK_out_ENA_NS_chr35.ld",sep="")


LD_data_even_na_ns <- rbind( LD_even_na_ns_chr1,  LD_even_na_ns_chr2, LD_even_na_ns_chr3, LD_even_na_ns_chr4, LD_even_na_ns_chr5, LD_even_na_ns_chr6,  LD_even_na_ns_chr7, LD_even_na_ns_chr8, LD_even_na_ns_chr9,
                             LD_even_na_ns_chr10, LD_even_na_ns_chr11,  LD_even_na_ns_chr12,  LD_even_na_ns_chr13, LD_even_na_ns_chr14,  LD_even_na_ns_chr15,  LD_even_na_ns_chr16,  LD_even_na_ns_chr17,  LD_even_na_ns_chr18,  LD_even_na_ns_chr19,
                             LD_even_na_ns_chr20,  LD_even_na_ns_chr21,  LD_even_na_ns_chr22, LD_even_na_ns_chr23,  LD_even_na_ns_chr24, LD_even_na_ns_chr25,  LD_even_na_ns_chr26,  LD_even_na_ns_chr27,  LD_even_na_ns_chr28,  LD_even_na_ns_chr29,  LD_even_na_ns_chr30,  LD_even_na_ns_chr31,  LD_even_na_ns_chr32,  LD_even_na_ns_chr33, LD_even_na_ns_chr34, LD_even_na_ns_chr0)
dim(LD_data_even_na_ns)
rm( LD_even_na_ns_chr1,  LD_even_na_ns_chr2, LD_even_na_ns_chr3, LD_even_na_ns_chr4, LD_even_na_ns_chr5, LD_even_na_ns_chr6,  LD_even_na_ns_chr7, LD_even_na_ns_chr8, LD_even_na_ns_chr9,
    LD_even_na_ns_chr10, LD_even_na_ns_chr11,  LD_even_na_ns_chr12,  LD_even_na_ns_chr13, LD_even_na_ns_chr14,  LD_even_na_ns_chr15,  LD_even_na_ns_chr16,  LD_even_na_ns_chr17,  LD_even_na_ns_chr18,  LD_even_na_ns_chr19,
    LD_even_na_ns_chr20,  LD_even_na_ns_chr21,  LD_even_na_ns_chr22, LD_even_na_ns_chr23,  LD_even_na_ns_chr24, LD_even_na_ns_chr25,  LD_even_na_ns_chr26,  LD_even_na_ns_chr27,  LD_even_na_ns_chr28,  LD_even_na_ns_chr29,  LD_even_na_ns_chr30,  LD_even_na_ns_chr31,  LD_even_na_ns_chr32,  LD_even_na_ns_chr33, LD_even_na_ns_chr34, LD_even_na_ns_chr0)

### to release the memory after deleting all those individual chromosome files.

gc()  

#######################Calculate the rank of each of the LDs per chromosome for all  groups.  

#Add a column to each dataframe that is for the rank of the R2 per chromosome
LD_data_even_A$Rank <- NA
LD_data_even_NA$Rank <-NA
LD_data_odd_a$Rank <- NA
LD_data_odd_na$Rank <- NA
LD_data_even_na_ns$Rank <- NA
LD_data_odd_na_ns$Rank <- NA

#this uses dyplr to rank the r2 per each chromosome, then orders the rank within each chromosome

#LD_data_even_A 
LD_data_even_A <- LD_data_even_A  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LD_data_even_A <- as.data.frame(LD_data_even_A)
LD_data_even_A <- LD_data_even_A[order(LD_data_even_A[,1],LD_data_even_A[,8]), ]
head(LD_data_even_A)

#LD_data_even_NA
LD_data_even_NA <- LD_data_even_NA  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LD_data_even_NA <- as.data.frame(LD_data_even_NA)
LD_data_even_NA <- LD_data_even_NA[order(LD_data_even_NA[,1],LD_data_even_NA[,8]), ]
head(LD_data_even_NA)

#LD_data_odd_a 
LD_data_odd_a <- LD_data_odd_a  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LD_data_odd_a <- as.data.frame(LD_data_odd_a)
LD_data_odd_a <- LD_data_odd_a[order(LD_data_odd_a[,1],LD_data_odd_a[,8]), ]
head(LD_data_odd_a)

#LD_data_odd_na
LD_data_odd_na <- LD_data_odd_na  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LD_data_odd_na <- as.data.frame(LD_data_odd_na)
LD_data_odd_na <- LD_data_odd_na[order(LD_data_odd_na[,1],LD_data_odd_na[,8]), ]
head(LD_data_odd_na)

#LD_data_even_na_ns
LD_data_even_na_ns <- LD_data_even_na_ns  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LD_data_even_na_ns <- as.data.frame(LD_data_even_na_ns)
LD_data_even_na_ns <- LD_data_even_na_ns[order(LD_data_even_na_ns[,1],LD_data_even_na_ns[,8]), ]
head(LD_data_even_na_ns)

#LD_data_odd_na_ns
LD_data_odd_na_ns <- LD_data_odd_na_ns  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LD_data_odd_na_ns <- as.data.frame(LD_data_odd_na_ns)
LD_data_odd_na_ns <- LD_data_odd_na_ns[order(LD_data_odd_na_ns[,1],LD_data_odd_na_ns[,8]), ]
head(LD_data_odd_na_ns)

#Calculate the physical distance between each of the paired loci
#Add a column to each dataframe that is the physical distance between each of the loci in the pairs
LD_data_even_A$Distance <- abs(LD_data_even_A$BP_A - LD_data_even_A$BP_B)
LD_data_even_NA$Distance <- abs(LD_data_even_NA$BP_A - LD_data_even_NA$BP_B)
LD_data_odd_a$Distance <- abs(LD_data_odd_a$BP_A - LD_data_odd_a$BP_B)
LD_data_odd_na$Distance <-  abs(LD_data_odd_na$BP_A - LD_data_odd_na$BP_B)
LD_data_even_na_ns$Distance <- abs(LD_data_even_na_ns$BP_A - LD_data_even_na_ns$BP_B)
LD_data_odd_na_ns$Distance <-  abs(LD_data_odd_na_ns$BP_A - LD_data_odd_na_ns$BP_B)



##########################################################
#For running the plots without having to import all the original data 
#these datasets were exported at the end of the code originally

LD_data_even_A <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_even_A.txt",  header=TRUE, stringsAsFactors = FALSE, na.strings = "-" )
LD_data_even_NA <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_even_NA.txt",  header=TRUE, stringsAsFactors = FALSE, na.strings = "-" )
LD_data_odd_a <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_odd_a.txt",  header=TRUE, stringsAsFactors = FALSE, na.strings = "-" )
LD_data_odd_na <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_odd_na.txt",  header=TRUE, stringsAsFactors = FALSE, na.strings = "-" )
LD_data_even_na_ns <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_even_na_ns.txt",  header=TRUE, stringsAsFactors = FALSE, na.strings = "-" )
LD_data_odd_na_ns <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_odd_na_ns.txt",  header=TRUE, stringsAsFactors = FALSE, na.strings = "-" )

###################################################
#Facet wrap ggplot the rank of r2 per each chromosome for each of the 6 groups

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Ranked_LD_regional_groupings.pdf", width = 9, height = 7)

# plot Individual Genotype Rate Post_Filtered_data, Ignore the Unassigned LD
ggplot(LD_data_even_A[(LD_data_even_A$CHR_A != "35") & (LD_data_even_A$R2 >=0.05),], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_even_A[(LD_data_even_A$CHR_A != "35") & (LD_data_even_A$R2 >=0.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale= "free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value >0.05 by chromosome for Even Asian samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LD_data_even_NA[(LD_data_even_NA$CHR_A != "35") & (LD_data_even_NA$R2 >=0.05),], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_even_NA[(LD_data_even_NA$CHR_A != "35") & (LD_data_even_NA$R2 >=0.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale= "free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value >0.05 by chromosome for Even North American samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LD_data_odd_a[(LD_data_odd_a$CHR_A != "35") & (LD_data_odd_a$R2 >=0.05),], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_odd_a[(LD_data_odd_a$CHR_A != "35") & (LD_data_odd_a$R2 >=0.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale= "free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value >0.05 by chromosome for Odd Asian samples") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LD_data_odd_na[(LD_data_odd_na$CHR_A != "35") & (LD_data_odd_na$R2 >=0.05),], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_odd_na[(LD_data_odd_na$CHR_A != "35") & (LD_data_odd_na$R2 >=0.05),,], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale= "free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value >0.05 by chromosome for Odd North American samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LD_data_even_na_ns[(LD_data_even_na_ns$CHR_A != "35") & (LD_data_even_na_ns$R2 >=0.05),], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_even_na_ns[(LD_data_even_na_ns$CHR_A != "35") & (LD_data_even_na_ns$R2 >=0.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale= "free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value >0.05 by chromosome for Even North American No Susitna samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LD_data_odd_na_ns[(LD_data_odd_na_ns$CHR_A != "35") & (LD_data_odd_na_ns$R2 >=0.05),], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_odd_na_ns[(LD_data_odd_na_ns$CHR_A != "35") & (LD_data_odd_na_ns$R2 >=0.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale= "free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value >0.05 by chromosome for Odd North American No Susitna samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

dev.off()

###############################################################
# plot R2 rank, Ignore the Unassigned LD
###THESE have to be individually saved as jpegs because they have too many points to be rendered in pdf each time. 

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Ranked_LD_reginal_LD_data_even_A.jpg", ggplot(LD_data_even_A[LD_data_even_A$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
         geom_point(data= LD_data_even_A[LD_data_even_A$CHR_A != "35",], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale = "free") + theme_bw() + guides(fill= FALSE) +
         ggtitle("Ranked LD R2 value by chromosome for Even Asian samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none"), width = 9,  height = 7,  dpi = 1200)

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Ranked_LD_reginal_LD_data_even_NA.jpg", ggplot(LD_data_even_NA[LD_data_even_NA$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
         geom_point(data= LD_data_even_NA[LD_data_even_NA$CHR_A != "35",], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale = "free") + theme_bw() + guides(fill= FALSE) +
         ggtitle("Ranked LD R2 value by chromosome for Even North American Samples ")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none"), width = 9,  height = 7,  dpi = 1200)

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Ranked_LD_reginal_LD_data_odd_a.jpg", ggplot(LD_data_odd_a[LD_data_odd_a$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
         geom_point(data= LD_data_odd_a[LD_data_odd_a$CHR_A != "35",], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale = "free") + theme_bw() + guides(fill= FALSE) +
         ggtitle("Ranked LD R2 value by chromosome for Odd Asian Samples") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none"), width = 9,  height = 7,  dpi = 1200)

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Ranked_LD_reginal_LD_data_odd_na.jpg", ggplot(LD_data_odd_na[LD_data_odd_na$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
         geom_point(data= LD_data_odd_na[LD_data_odd_na$CHR_A != "35",], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale = "free") + theme_bw() + guides(fill= FALSE) +
         ggtitle("Ranked LD R2 value by chromosome for Odd North American samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none"), width = 9,  height = 7,  dpi = 1200)

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Ranked_LD_reginal_LD_data_even_na_ns.jpg", ggplot(LD_data_even_na_ns[LD_data_even_na_ns$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
         geom_point(data= LD_data_even_na_ns[LD_data_even_na_ns$CHR_A != "35",], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale = "free") + theme_bw() + guides(fill= FALSE) +
         ggtitle("Ranked LD R2 value by chromosome for Even North American No Susitna samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none"), width = 9,  height = 7,  dpi = 1200)

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Ranked_LD_reginal_LD_data_odd_na_ns.jpg", ggplot(LD_data_odd_na_ns[LD_data_odd_na_ns$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
         geom_point(data= LD_data_odd_na_ns[LD_data_odd_na_ns$CHR_A != "35",], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale = "free") + theme_bw() + guides(fill= FALSE) +
         ggtitle("Ranked LD R2 value by chromosome for Odd North American No Susitna samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none"), width = 9,  height = 7,  dpi = 1200)

#################################################
# Plot LD decay plots for each of the loci, R2 > 0.05 only, and ignore the unassigned chromosome. 
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_Decay_by_CHR_regional_groupings.pdf", width = 9, height = 7)

ggplot(LD_data_even_A[(LD_data_even_A$CHR_A != "35") & (LD_data_even_A$R2 >=.05),], aes(x=Distance, y= R2 )) + 
  geom_point(data= LD_data_even_A[(LD_data_even_A$CHR_A != "35") & (LD_data_even_A$R2 >=.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale= "free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 value >0.05 by physical distance per chromosome for Even Asian samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LD_data_even_NA[(LD_data_even_NA$CHR_A != "35") & (LD_data_even_NA$R2 >=.05),], aes(x=Distance, y= R2 )) + 
  geom_point(data= LD_data_even_NA[(LD_data_even_NA$CHR_A != "35") & (LD_data_even_NA$R2 >=.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale= "free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 value >0.05 by physical distance per chromosome for Even North American lineage ")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LD_data_odd_a[(LD_data_odd_a$CHR_A != "35") & (LD_data_odd_a$R2 >=.05),], aes(x=Distance, y= R2 )) + 
  geom_point(data= LD_data_odd_a[(LD_data_odd_a$CHR_A != "35") & (LD_data_odd_a$R2 >=.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale= "free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 value >0.05 by physical distance per chromosome for Odd Asian lineage")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LD_data_odd_na[(LD_data_odd_na$CHR_A != "35") & (LD_data_odd_na$R2 >=.05),], aes(x=Distance, y= R2 )) + 
  geom_point(data= LD_data_odd_na[(LD_data_odd_na$CHR_A != "35") & (LD_data_odd_na$R2 >=.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale= "free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 value >0.05 by physical distance per chromosome for Odd North American samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LD_data_even_na_ns[(LD_data_even_na_ns$CHR_A != "35") & (LD_data_even_na_ns$R2 >=.05),], aes(x=Distance, y= R2 )) + 
  geom_point(data= LD_data_even_na_ns[(LD_data_even_na_ns$CHR_A != "35") & (LD_data_even_na_ns$R2 >=.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale= "free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 value >0.05 by physical distance per chromosome for Even North American No Susitna samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LD_data_odd_na_ns[(LD_data_odd_na_ns$CHR_A != "35") & (LD_data_odd_na_ns$R2 >=.05),], aes(x=Distance, y= R2 )) + 
  geom_point(data= LD_data_odd_na_ns[(LD_data_odd_na_ns$CHR_A != "35") & (LD_data_odd_na_ns$R2 >=.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale= "free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 value >0.05 by physical distance per chromosome for Odd North American No Susitna samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

dev.off()

#clear out any memory that we dont need after making plots. 
rm(LD_data_even_A, LD_data_even_NA, LD_data_odd_a , LD_data_odd_na , LD_data_even_na_ns, LD_data_odd_na_ns)

gc()

#######################
#save the dataframes that were made here so that they can be imported separately. 
####Write the combined data file for each LD R2 group and a filtered version that is easier for plotting LD  to files to be used in other R scripts

#LD_data_even_A
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_even_A.txt", "wb")
write.table(LD_data_even_A,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

LD_data_even_A_.2_BPA_BPB_R2 <- LD_data_even_A[LD_data_even_A$R2 >=0.2, c("CHR_A","BP_A", "BP_B", "R2")]
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_even_A_.2_BPA_BPB_R2.txt", "wb")
write.table(LD_data_even_A_.2_BPA_BPB_R2,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

rm(LD_data_even_A, LD_data_even_A_.2_BPA_BPB_R2)
gc()

#LD_data_even_NA
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_even_NA.txt", "wb")
write.table(LD_data_even_NA,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

LD_data_even_NA_.2_BPA_BPB_R2 <- LD_data_even_NA[LD_data_even_NA$R2 >=0.2, c("CHR_A","BP_A", "BP_B", "R2")]
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_even_NA_.2_BPA_BPB_R2.txt", "wb")
write.table(LD_data_even_NA_.2_BPA_BPB_R2,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

rm(LD_data_even_NA, LD_data_even_NA_.2_BPA_BPB_R2)
gc()

#LD_data_odd_a
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_odd_a.txt", "wb")
write.table(LD_data_odd_a,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

LD_data_odd_a_.2_BPA_BPB_R2 <- LD_data_odd_a[LD_data_odd_a$R2 >=0.2, c("CHR_A","BP_A", "BP_B", "R2")]
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_odd_a_.2_BPA_BPB_R2.txt", "wb")
write.table(LD_data_odd_a_.2_BPA_BPB_R2,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

rm(LD_data_odd_a, LD_data_odd_a_.2_BPA_BPB_R2)
gc()


#LD_data_odd_na
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_odd_na.txt", "wb")
write.table(LD_data_odd_na,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

LD_data_odd_na_.2_BPA_BPB_R2 <- LD_data_odd_na[LD_data_odd_na$R2 >=0.2, c("CHR_A","BP_A", "BP_B", "R2")]
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_odd_na_.2_BPA_BPB_R2.txt", "wb")
write.table(LD_data_odd_na_.2_BPA_BPB_R2,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

rm(LD_data_odd_na, LD_data_odd_na_.2_BPA_BPB_R2)
gc()


#LD_data_even_na_ns
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_even_na_ns.txt", "wb")
write.table(LD_data_even_na_ns,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

LD_data_even_na_ns_.2_BPA_BPB_R2 <- LD_data_even_na_ns[LD_data_even_na_ns$R2 >=0.2, c("CHR_A","BP_A", "BP_B", "R2")]
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_even_na_ns_.2_BPA_BPB_R2.txt", "wb")
write.table(LD_data_even_na_ns_.2_BPA_BPB_R2,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

rm(LD_data_even_na_ns, LD_data_even_na_ns_.2_BPA_BPB_R2)
gc()

#LD_data_odd_na_ns
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_odd_na_ns.txt", "wb")
write.table(LD_data_odd_na_ns,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

LD_data_odd_na_ns_.2_BPA_BPB_R2 <- LD_data_odd_na_ns[LD_data_odd_na_ns$R2 >=0.2, c("CHR_A","BP_A", "BP_B", "R2")]
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_odd_na_ns_.2_BPA_BPB_R2.txt", "wb")
write.table(LD_data_odd_na_ns_.2_BPA_BPB_R2,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

rm(LD_data_odd_na_ns, LD_data_odd_na_ns_.2_BPA_BPB_R2)
gc()

#save this workspace independently !!






