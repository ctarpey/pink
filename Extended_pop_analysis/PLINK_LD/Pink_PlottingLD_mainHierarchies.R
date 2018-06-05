### Plotting the LD r2 from PLINK 
###    for the full dataset of One SNP per Tag- 23759 loci. 
###    THE CODE FOR THE PLOTS OF THE ACTUAL LD /CHROM IS IN _2 
###   Carolyn Tarpey | May 2018
### ---------------------------------------

#The full 13 population groups takes up too much memory, so use the workspace Pink_PlottingLD_mainHierarchies.R.RDATA for these 7 groups
# Load the rest of the groups in the _Regional_groups R scripts, and the NA_subset for the breakout of the North American populations. 

#This R code was originally written by Garrett McKinney to plot the output of PLINK LD r2 for each chromosome. 
#It requires the LD output from PLINK. We used the alignment of the pinks to Chinook to get the Chromosome 
#and position assignments here, so there are 34 Chromosomes, and chromosome 35 are those that aligned to the unassigned 

#TO MAKE THE PLOTS OF THE LD PER CHROMOSOME, USE THIS SCRIPT TO IMPORT THE DATA AND THE Pink_PlottingLD_2_mainHierarchies.R 
## and Pink_PlottingLD_2_thumbs_mainHierarchies.R TO MAKE ADDITIONAL PLOTS

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
#### ALL 
LDdata_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr1.ld",sep="")
LDdata_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr2.ld",sep="")
LDdata_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr3.ld",sep="")
LDdata_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr4.ld",sep="")
LDdata_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr5.ld",sep="")
LDdata_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr6.ld",sep="")
LDdata_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr7.ld",sep="")
LDdata_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr8.ld",sep="")
LDdata_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr9.ld",sep="")
LDdata_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr10.ld",sep="")
LDdata_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr11.ld",sep="")
LDdata_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr12.ld",sep="")
LDdata_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr13.ld",sep="")
LDdata_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr14.ld",sep="")
LDdata_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr15.ld",sep="")
LDdata_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr16.ld",sep="")
LDdata_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr17.ld",sep="")
LDdata_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr18.ld",sep="")
LDdata_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr19.ld",sep="")
LDdata_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr20.ld",sep="")
LDdata_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr21.ld",sep="")
LDdata_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr22.ld",sep="")
LDdata_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr23.ld",sep="")
LDdata_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr24.ld",sep="")
LDdata_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr25.ld",sep="")
LDdata_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr26.ld",sep="")
LDdata_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr27.ld",sep="")
LDdata_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr28.ld",sep="")
LDdata_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr29.ld",sep="")
LDdata_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr30.ld",sep="")
LDdata_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr31.ld",sep="")
LDdata_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr32.ld",sep="")
LDdata_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr33.ld",sep="")
LDdata_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr34.ld",sep="")
LDdata_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ALL/PLINK_out_chr35.ld",sep="")

#combine each of the chromosome specific LDs into one big table, then deletes the idividual 
LDdata_all <- rbind( LDdata_chr1,  LDdata_chr2, LDdata_chr3, LDdata_chr4,LDdata_chr5,LDdata_chr6,  LDdata_chr7,  LDdata_chr8,  LDdata_chr9,
                     LDdata_chr10, LDdata_chr11,  LDdata_chr12,  LDdata_chr13, LDdata_chr14,  LDdata_chr15,  LDdata_chr16,  LDdata_chr17,  LDdata_chr18,  LDdata_chr19,
                     LDdata_chr20,  LDdata_chr21,  LDdata_chr22, LDdata_chr23,  LDdata_chr24, LDdata_chr25,  LDdata_chr26,  LDdata_chr27,  LDdata_chr28,  LDdata_chr29,  LDdata_chr30,  LDdata_chr31,  LDdata_chr32,  LDdata_chr33, LDdata_chr34, LDdata_chr0)
dim(LDdata_all)
head(LDdata_all)
tail(LDdata_all)

rm(LDdata_chr1,  LDdata_chr2, LDdata_chr3, LDdata_chr4,LDdata_chr5,LDdata_chr6,  LDdata_chr7,  LDdata_chr8,  LDdata_chr9,
   LDdata_chr10, LDdata_chr11,  LDdata_chr12,  LDdata_chr13, LDdata_chr14,  LDdata_chr15,  LDdata_chr16,  LDdata_chr17,  LDdata_chr18,  LDdata_chr19,
   LDdata_chr20,  LDdata_chr21,  LDdata_chr22, LDdata_chr23,  LDdata_chr24, LDdata_chr25,  LDdata_chr26,  LDdata_chr27,  LDdata_chr28,  LDdata_chr29,  LDdata_chr30,  LDdata_chr31,  LDdata_chr32,  LDdata_chr33, LDdata_chr34, LDdata_chr0)


#### Even
#load PLINK LD output files
LD_even_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr1.ld",sep="")
LD_even_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr2.ld",sep="")
LD_even_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr3.ld",sep="")
LD_even_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr4.ld",sep="")
LD_even_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr5.ld",sep="")
LD_even_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr6.ld",sep="")
LD_even_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr7.ld",sep="")
LD_even_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr8.ld",sep="")
LD_even_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr9.ld",sep="")
LD_even_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr10.ld",sep="")
LD_even_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr11.ld",sep="")
LD_even_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr12.ld",sep="")
LD_even_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr13.ld",sep="")
LD_even_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr14.ld",sep="")
LD_even_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr15.ld",sep="")
LD_even_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr16.ld",sep="")
LD_even_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr17.ld",sep="")
LD_even_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr18.ld",sep="")
LD_even_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr19.ld",sep="")
LD_even_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr20.ld",sep="")
LD_even_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr21.ld",sep="")
LD_even_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr22.ld",sep="")
LD_even_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr23.ld",sep="")
LD_even_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr24.ld",sep="")
LD_even_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr25.ld",sep="")
LD_even_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr26.ld",sep="")
LD_even_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr27.ld",sep="")
LD_even_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr28.ld",sep="")
LD_even_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr29.ld",sep="")
LD_even_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr30.ld",sep="")
LD_even_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr31.ld",sep="")
LD_even_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr32.ld",sep="")
LD_even_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr33.ld",sep="")
LD_even_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr34.ld",sep="")
LD_even_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN/PLINK_out_EVEN_chr35.ld",sep="")

LDdata_even <- rbind( LD_even_chr1,  LD_even_chr2, LD_even_chr3, LD_even_chr4, LD_even_chr5, LD_even_chr6,  LD_even_chr7, LD_even_chr8, LD_even_chr9,
                      LD_even_chr10, LD_even_chr11,  LD_even_chr12,  LD_even_chr13, LD_even_chr14,  LD_even_chr15,  LD_even_chr16,  LD_even_chr17,  LD_even_chr18,  LD_even_chr19,
                      LD_even_chr20,  LD_even_chr21,  LD_even_chr22, LD_even_chr23,  LD_even_chr24, LD_even_chr25,  LD_even_chr26,  LD_even_chr27,  LD_even_chr28,  LD_even_chr29,  LD_even_chr30,  LD_even_chr31,  LD_even_chr32,  LD_even_chr33, LD_even_chr34, LD_even_chr0)
dim(LDdata_even)
rm( LD_even_chr1,  LD_even_chr2, LD_even_chr3, LD_even_chr4, LD_even_chr5, LD_even_chr6,  LD_even_chr7, LD_even_chr8, LD_even_chr9,
    LD_even_chr10, LD_even_chr11,  LD_even_chr12,  LD_even_chr13, LD_even_chr14,  LD_even_chr15,  LD_even_chr16,  LD_even_chr17,  LD_even_chr18,  LD_even_chr19,
    LD_even_chr20,  LD_even_chr21,  LD_even_chr22, LD_even_chr23,  LD_even_chr24, LD_even_chr25,  LD_even_chr26,  LD_even_chr27,  LD_even_chr28,  LD_even_chr29,  LD_even_chr30,  LD_even_chr31,  LD_even_chr32,  LD_even_chr33, LD_even_chr34, LD_even_chr0)

#### Odd
#load PLINK LD output files
LD_odd_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr1.ld",sep="")
LD_odd_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr2.ld",sep="")
LD_odd_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr3.ld",sep="")
LD_odd_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr4.ld",sep="")
LD_odd_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr5.ld",sep="")
LD_odd_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr6.ld",sep="")
LD_odd_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr7.ld",sep="")
LD_odd_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr8.ld",sep="")
LD_odd_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr9.ld",sep="")
LD_odd_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr10.ld",sep="")
LD_odd_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr11.ld",sep="")
LD_odd_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr12.ld",sep="")
LD_odd_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr13.ld",sep="")
LD_odd_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr14.ld",sep="")
LD_odd_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr15.ld",sep="")
LD_odd_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr16.ld",sep="")
LD_odd_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr17.ld",sep="")
LD_odd_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr18.ld",sep="")
LD_odd_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr19.ld",sep="")
LD_odd_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr20.ld",sep="")
LD_odd_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr21.ld",sep="")
LD_odd_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr22.ld",sep="")
LD_odd_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr23.ld",sep="")
LD_odd_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr24.ld",sep="")
LD_odd_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr25.ld",sep="")
LD_odd_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr26.ld",sep="")
LD_odd_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr27.ld",sep="")
LD_odd_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr28.ld",sep="")
LD_odd_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr29.ld",sep="")
LD_odd_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr30.ld",sep="")
LD_odd_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr31.ld",sep="")
LD_odd_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr32.ld",sep="")
LD_odd_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr33.ld",sep="")
LD_odd_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr34.ld",sep="")
LD_odd_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD/PLINK_out_ODD_chr35.ld",sep="")


LDdata_odd <- rbind( LD_odd_chr1,  LD_odd_chr2, LD_odd_chr3, LD_odd_chr4,LD_odd_chr5,LD_odd_chr6,  LD_odd_chr7,  LD_odd_chr8,  LD_odd_chr9,
                     LD_odd_chr10, LD_odd_chr11,  LD_odd_chr12,  LD_odd_chr13, LD_odd_chr14,  LD_odd_chr15,  LD_odd_chr16,  LD_odd_chr17,  LD_odd_chr18,  LD_odd_chr19,
                     LD_odd_chr20,  LD_odd_chr21,  LD_odd_chr22, LD_odd_chr23,  LD_odd_chr24, LD_odd_chr25,  LD_odd_chr26,  LD_odd_chr27,  LD_odd_chr28,  LD_odd_chr29,  LD_odd_chr30,  LD_odd_chr31,  LD_odd_chr32,  LD_odd_chr33, LD_odd_chr34, LD_odd_chr0)
dim(LDdata_odd)
rm( LD_odd_chr1,  LD_odd_chr2, LD_odd_chr3, LD_odd_chr4,LD_odd_chr5,LD_odd_chr6,  LD_odd_chr7,  LD_odd_chr8,  LD_odd_chr9,
    LD_odd_chr10, LD_odd_chr11,  LD_odd_chr12,  LD_odd_chr13, LD_odd_chr14,  LD_odd_chr15,  LD_odd_chr16,  LD_odd_chr17,  LD_odd_chr18,  LD_odd_chr19,
    LD_odd_chr20,  LD_odd_chr21,  LD_odd_chr22, LD_odd_chr23,  LD_odd_chr24, LD_odd_chr25,  LD_odd_chr26,  LD_odd_chr27,  LD_odd_chr28,  LD_odd_chr29,  LD_odd_chr30,  LD_odd_chr31,  LD_odd_chr32,  LD_odd_chr33, LD_odd_chr34, LD_odd_chr0)


########################THESE HAVE NO SUSITNA POPULATIONS INCLUDED :
#### Odd NO SUSITNA
#load PLINK LD output files
LD_odd_ns_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr1.ld",sep="")
LD_odd_ns_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr2.ld",sep="")
LD_odd_ns_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr3.ld",sep="")
LD_odd_ns_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr4.ld",sep="")
LD_odd_ns_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr5.ld",sep="")
LD_odd_ns_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr6.ld",sep="")
LD_odd_ns_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr7.ld",sep="")
LD_odd_ns_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr8.ld",sep="")
LD_odd_ns_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr9.ld",sep="")
LD_odd_ns_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr10.ld",sep="")
LD_odd_ns_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr11.ld",sep="")
LD_odd_ns_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr12.ld",sep="")
LD_odd_ns_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr13.ld",sep="")
LD_odd_ns_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr14.ld",sep="")
LD_odd_ns_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr15.ld",sep="")
LD_odd_ns_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr16.ld",sep="")
LD_odd_ns_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr17.ld",sep="")
LD_odd_ns_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr18.ld",sep="")
LD_odd_ns_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr19.ld",sep="")
LD_odd_ns_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr20.ld",sep="")
LD_odd_ns_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr21.ld",sep="")
LD_odd_ns_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr22.ld",sep="")
LD_odd_ns_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr23.ld",sep="")
LD_odd_ns_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr24.ld",sep="")
LD_odd_ns_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr25.ld",sep="")
LD_odd_ns_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr26.ld",sep="")
LD_odd_ns_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr27.ld",sep="")
LD_odd_ns_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr28.ld",sep="")
LD_odd_ns_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr29.ld",sep="")
LD_odd_ns_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr30.ld",sep="")
LD_odd_ns_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr31.ld",sep="")
LD_odd_ns_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr32.ld",sep="")
LD_odd_ns_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr33.ld",sep="")
LD_odd_ns_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr34.ld",sep="")
LD_odd_ns_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ODD_NS/PLINK_out_ODD_NS_chr35.ld",sep="")


LDdata_odd_ns <- rbind( LD_odd_ns_chr1,  LD_odd_ns_chr2, LD_odd_ns_chr3, LD_odd_ns_chr4,LD_odd_ns_chr5,LD_odd_ns_chr6,  LD_odd_ns_chr7,  LD_odd_ns_chr8,  LD_odd_ns_chr9,
                        LD_odd_ns_chr10, LD_odd_ns_chr11,  LD_odd_ns_chr12,  LD_odd_ns_chr13, LD_odd_ns_chr14,  LD_odd_ns_chr15,  LD_odd_ns_chr16,  LD_odd_ns_chr17,  LD_odd_ns_chr18,  LD_odd_ns_chr19,
                        LD_odd_ns_chr20,  LD_odd_ns_chr21,  LD_odd_ns_chr22, LD_odd_ns_chr23,  LD_odd_ns_chr24, LD_odd_ns_chr25,  LD_odd_ns_chr26,  LD_odd_ns_chr27,  LD_odd_ns_chr28,  LD_odd_ns_chr29,  LD_odd_ns_chr30,  LD_odd_ns_chr31,  LD_odd_ns_chr32,  LD_odd_ns_chr33, LD_odd_ns_chr34, LD_odd_ns_chr0)
dim(LDdata_odd_ns)
rm( LD_odd_ns_chr1,  LD_odd_ns_chr2, LD_odd_ns_chr3, LD_odd_ns_chr4,LD_odd_ns_chr5,LD_odd_ns_chr6,  LD_odd_ns_chr7,  LD_odd_ns_chr8,  LD_odd_ns_chr9,
    LD_odd_ns_chr10, LD_odd_ns_chr11,  LD_odd_ns_chr12,  LD_odd_ns_chr13, LD_odd_ns_chr14,  LD_odd_ns_chr15,  LD_odd_ns_chr16,  LD_odd_ns_chr17,  LD_odd_ns_chr18,  LD_odd_ns_chr19,
    LD_odd_ns_chr20,  LD_odd_ns_chr21,  LD_odd_ns_chr22, LD_odd_ns_chr23,  LD_odd_ns_chr24, LD_odd_ns_chr25,  LD_odd_ns_chr26,  LD_odd_ns_chr27,  LD_odd_ns_chr28,  LD_odd_ns_chr29,  LD_odd_ns_chr30,  LD_odd_ns_chr31,  LD_odd_ns_chr32,  LD_odd_ns_chr33, LD_odd_ns_chr34, LD_odd_ns_chr0)


#### EVEN NO SUSITNA
#load PLINK LD output files
LD_even_ns_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr1.ld",sep="")
LD_even_ns_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr2.ld",sep="")
LD_even_ns_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr3.ld",sep="")
LD_even_ns_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr4.ld",sep="")
LD_even_ns_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr5.ld",sep="")
LD_even_ns_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr6.ld",sep="")
LD_even_ns_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr7.ld",sep="")
LD_even_ns_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr8.ld",sep="")
LD_even_ns_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr9.ld",sep="")
LD_even_ns_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr10.ld",sep="")
LD_even_ns_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr11.ld",sep="")
LD_even_ns_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr12.ld",sep="")
LD_even_ns_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr13.ld",sep="")
LD_even_ns_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr14.ld",sep="")
LD_even_ns_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr15.ld",sep="")
LD_even_ns_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr16.ld",sep="")
LD_even_ns_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr17.ld",sep="")
LD_even_ns_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr18.ld",sep="")
LD_even_ns_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr19.ld",sep="")
LD_even_ns_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr20.ld",sep="")
LD_even_ns_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr21.ld",sep="")
LD_even_ns_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr22.ld",sep="")
LD_even_ns_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr23.ld",sep="")
LD_even_ns_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr24.ld",sep="")
LD_even_ns_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr25.ld",sep="")
LD_even_ns_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr26.ld",sep="")
LD_even_ns_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr27.ld",sep="")
LD_even_ns_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr28.ld",sep="")
LD_even_ns_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr29.ld",sep="")
LD_even_ns_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr30.ld",sep="")
LD_even_ns_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr31.ld",sep="")
LD_even_ns_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr32.ld",sep="")
LD_even_ns_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr33.ld",sep="")
LD_even_ns_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr34.ld",sep="")
LD_even_ns_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_EVEN_NS/PLINK_out_EVEN_NS_chr35.ld",sep="")


LDdata_even_ns <- rbind( LD_even_ns_chr1,  LD_even_ns_chr2, LD_even_ns_chr3, LD_even_ns_chr4, LD_even_ns_chr5, LD_even_ns_chr6,  LD_even_ns_chr7, LD_even_ns_chr8, LD_even_ns_chr9,
                         LD_even_ns_chr10, LD_even_ns_chr11,  LD_even_ns_chr12,  LD_even_ns_chr13, LD_even_ns_chr14,  LD_even_ns_chr15,  LD_even_ns_chr16,  LD_even_ns_chr17,  LD_even_ns_chr18,  LD_even_ns_chr19,
                         LD_even_ns_chr20,  LD_even_ns_chr21,  LD_even_ns_chr22, LD_even_ns_chr23,  LD_even_ns_chr24, LD_even_ns_chr25,  LD_even_ns_chr26,  LD_even_ns_chr27,  LD_even_ns_chr28,  LD_even_ns_chr29,  LD_even_ns_chr30,  LD_even_ns_chr31,  LD_even_ns_chr32,  LD_even_ns_chr33, LD_even_ns_chr34, LD_even_ns_chr0)
dim(LDdata_even_ns)
rm( LD_even_ns_chr1,  LD_even_ns_chr2, LD_even_ns_chr3, LD_even_ns_chr4, LD_even_ns_chr5, LD_even_ns_chr6,  LD_even_ns_chr7, LD_even_ns_chr8, LD_even_ns_chr9,
    LD_even_ns_chr10, LD_even_ns_chr11,  LD_even_ns_chr12,  LD_even_ns_chr13, LD_even_ns_chr14,  LD_even_ns_chr15,  LD_even_ns_chr16,  LD_even_ns_chr17,  LD_even_ns_chr18,  LD_even_ns_chr19,
    LD_even_ns_chr20,  LD_even_ns_chr21,  LD_even_ns_chr22, LD_even_ns_chr23,  LD_even_ns_chr24, LD_even_ns_chr25,  LD_even_ns_chr26,  LD_even_ns_chr27,  LD_even_ns_chr28,  LD_even_ns_chr29,  LD_even_ns_chr30,  LD_even_ns_chr31,  LD_even_ns_chr32,  LD_even_ns_chr33, LD_even_ns_chr34, LD_even_ns_chr0)

#### NORTH AMERICA ALL 
#load PLINK LD output files
LD_na_all_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr1.ld",sep="")
LD_na_all_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr2.ld",sep="")
LD_na_all_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr3.ld",sep="")
LD_na_all_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr4.ld",sep="")
LD_na_all_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr5.ld",sep="")
LD_na_all_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr6.ld",sep="")
LD_na_all_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr7.ld",sep="")
LD_na_all_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr8.ld",sep="")
LD_na_all_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr9.ld",sep="")
LD_na_all_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr10.ld",sep="")
LD_na_all_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr11.ld",sep="")
LD_na_all_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr12.ld",sep="")
LD_na_all_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr13.ld",sep="")
LD_na_all_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr14.ld",sep="")
LD_na_all_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr15.ld",sep="")
LD_na_all_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr16.ld",sep="")
LD_na_all_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr17.ld",sep="")
LD_na_all_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr18.ld",sep="")
LD_na_all_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr19.ld",sep="")
LD_na_all_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr20.ld",sep="")
LD_na_all_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr21.ld",sep="")
LD_na_all_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr22.ld",sep="")
LD_na_all_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr23.ld",sep="")
LD_na_all_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr24.ld",sep="")
LD_na_all_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr25.ld",sep="")
LD_na_all_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr26.ld",sep="")
LD_na_all_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr27.ld",sep="")
LD_na_all_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr28.ld",sep="")
LD_na_all_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr29.ld",sep="")
LD_na_all_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr30.ld",sep="")
LD_na_all_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr31.ld",sep="")
LD_na_all_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr32.ld",sep="")
LD_na_all_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr33.ld",sep="")
LD_na_all_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr34.ld",sep="")
LD_na_all_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_NA_ALL/PLINK_NA_out_chr35.ld",sep="")

LD_data_na_all <- rbind( LD_na_all_chr1,  LD_na_all_chr2, LD_na_all_chr3, LD_na_all_chr4, LD_na_all_chr5, LD_na_all_chr6,  LD_na_all_chr7, LD_na_all_chr8, LD_na_all_chr9,
                         LD_na_all_chr10, LD_na_all_chr11,  LD_na_all_chr12,  LD_na_all_chr13, LD_na_all_chr14,  LD_na_all_chr15,  LD_na_all_chr16,  LD_na_all_chr17,  LD_na_all_chr18,  LD_na_all_chr19,
                         LD_na_all_chr20,  LD_na_all_chr21,  LD_na_all_chr22, LD_na_all_chr23,  LD_na_all_chr24, LD_na_all_chr25,  LD_na_all_chr26,  LD_na_all_chr27,  LD_na_all_chr28,  LD_na_all_chr29,  LD_na_all_chr30,  LD_na_all_chr31,  LD_na_all_chr32,  LD_na_all_chr33, LD_na_all_chr34, LD_na_all_chr0)
dim(LD_data_na_all)
rm( LD_na_all_chr1,  LD_na_all_chr2, LD_na_all_chr3, LD_na_all_chr4, LD_na_all_chr5, LD_na_all_chr6,  LD_na_all_chr7, LD_na_all_chr8, LD_na_all_chr9,
    LD_na_all_chr10, LD_na_all_chr11,  LD_na_all_chr12,  LD_na_all_chr13, LD_na_all_chr14,  LD_na_all_chr15,  LD_na_all_chr16,  LD_na_all_chr17,  LD_na_all_chr18,  LD_na_all_chr19,
    LD_na_all_chr20,  LD_na_all_chr21,  LD_na_all_chr22, LD_na_all_chr23,  LD_na_all_chr24, LD_na_all_chr25,  LD_na_all_chr26,  LD_na_all_chr27,  LD_na_all_chr28,  LD_na_all_chr29,  LD_na_all_chr30,  LD_na_all_chr31,  LD_na_all_chr32,  LD_na_all_chr33, LD_na_all_chr34, LD_na_all_chr0)

#### EVEN North America NO SUSITNA
#load PLINK LD output files
LD_a_all_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr1.ld",sep="")
LD_a_all_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr2.ld",sep="")
LD_a_all_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr3.ld",sep="")
LD_a_all_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr4.ld",sep="")
LD_a_all_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr5.ld",sep="")
LD_a_all_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr6.ld",sep="")
LD_a_all_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr7.ld",sep="")
LD_a_all_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr8.ld",sep="")
LD_a_all_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr9.ld",sep="")
LD_a_all_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr10.ld",sep="")
LD_a_all_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr11.ld",sep="")
LD_a_all_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr12.ld",sep="")
LD_a_all_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr13.ld",sep="")
LD_a_all_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr14.ld",sep="")
LD_a_all_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr15.ld",sep="")
LD_a_all_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr16.ld",sep="")
LD_a_all_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr17.ld",sep="")
LD_a_all_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr18.ld",sep="")
LD_a_all_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr19.ld",sep="")
LD_a_all_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr20.ld",sep="")
LD_a_all_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr21.ld",sep="")
LD_a_all_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr22.ld",sep="")
LD_a_all_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr23.ld",sep="")
LD_a_all_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr24.ld",sep="")
LD_a_all_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr25.ld",sep="")
LD_a_all_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr26.ld",sep="")
LD_a_all_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr27.ld",sep="")
LD_a_all_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr28.ld",sep="")
LD_a_all_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr29.ld",sep="")
LD_a_all_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr30.ld",sep="")
LD_a_all_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr31.ld",sep="")
LD_a_all_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr32.ld",sep="")
LD_a_all_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr33.ld",sep="")
LD_a_all_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr34.ld",sep="")
LD_a_all_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_23759_465_ASIA_ALL/PLINK_A_out_chr35.ld",sep="")

LD_data_a_all <- rbind(LD_a_all_chr1, LD_a_all_chr2,LD_a_all_chr3,LD_a_all_chr4,LD_a_all_chr5,LD_a_all_chr6, LD_a_all_chr7,LD_a_all_chr8,LD_a_all_chr9,
                       LD_a_all_chr10,LD_a_all_chr11, LD_a_all_chr12, LD_a_all_chr13,LD_a_all_chr14, LD_a_all_chr15, LD_a_all_chr16, LD_a_all_chr17, LD_a_all_chr18, LD_a_all_chr19,
                       LD_a_all_chr20, LD_a_all_chr21, LD_a_all_chr22,LD_a_all_chr23, LD_a_all_chr24,LD_a_all_chr25, LD_a_all_chr26, LD_a_all_chr27, LD_a_all_chr28, LD_a_all_chr29, LD_a_all_chr30, LD_a_all_chr31, LD_a_all_chr32, LD_a_all_chr33,LD_a_all_chr34,LD_a_all_chr0)
dim(LD_data_a_all)
rm(LD_a_all_chr1, LD_a_all_chr2,LD_a_all_chr3,LD_a_all_chr4,LD_a_all_chr5,LD_a_all_chr6, LD_a_all_chr7,LD_a_all_chr8,LD_a_all_chr9,
   LD_a_all_chr10,LD_a_all_chr11, LD_a_all_chr12, LD_a_all_chr13,LD_a_all_chr14, LD_a_all_chr15, LD_a_all_chr16, LD_a_all_chr17, LD_a_all_chr18, LD_a_all_chr19,
   LD_a_all_chr20, LD_a_all_chr21, LD_a_all_chr22,LD_a_all_chr23, LD_a_all_chr24,LD_a_all_chr25, LD_a_all_chr26, LD_a_all_chr27, LD_a_all_chr28, LD_a_all_chr29, LD_a_all_chr30, LD_a_all_chr31, LD_a_all_chr32, LD_a_all_chr33,LD_a_all_chr34,LD_a_all_chr0)

### to release the memory after deleting all those individual chromosome files.

gc()  

#######################Calculate the rank of each of the LDs per chromosome for all  groups. 

#Add a column to each dataframe that is for the rank of the R2 per chromosome
LDdata_all$Rank <- NA
LDdata_even$Rank <-NA
LDdata_odd$Rank <- NA
LDdata_even_ns$Rank <- NA
LDdata_odd_ns$Rank <- NA
LD_data_na_all$Rank <- NA
LD_data_a_all$Rank <- NA

#this uses dyplr to rank the r2 per each chromosome, then orders the rank within each chromosome
LDdata_all <- LDdata_all  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LDdata_all <- as.data.frame(LDdata_all)
LDdata_all <- LDdata_all[order(LDdata_all[,1],LDdata_all[,8]), ]
head(LDdata_all)

#LDdata_even 
LDdata_even <- LDdata_even  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LDdata_even <- as.data.frame(LDdata_even)
LDdata_even <- LDdata_even[order(LDdata_even[,1],LDdata_even[,8]), ]
head(LDdata_even)

#LDdata_odd 
LDdata_odd <- LDdata_odd  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LDdata_odd <- as.data.frame(LDdata_odd)
LDdata_odd <- LDdata_odd[order(LDdata_odd[,1],LDdata_odd[,8]), ]
head(LDdata_odd)

#LDdata_even_ns
LDdata_even_ns <- LDdata_even_ns  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LDdata_even_ns <- as.data.frame(LDdata_even_ns)
LDdata_even_ns <- LDdata_even_ns[order(LDdata_even_ns[,1],LDdata_even_ns[,8]), ]
head(LDdata_even_ns)

#LDdata_odd_ns
LDdata_odd_ns <- LDdata_odd_ns  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LDdata_odd_ns <- as.data.frame(LDdata_odd_ns)
LDdata_odd_ns <- LDdata_odd_ns[order(LDdata_odd_ns[,1],LDdata_odd_ns[,8]), ]
head(LDdata_odd_ns)

#LD_data_na_all
LD_data_na_all <- LD_data_na_all  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LD_data_na_all <- as.data.frame(LD_data_na_all)
LD_data_na_all <- LD_data_na_all[order(LD_data_na_all[,1],LD_data_na_all[,8]), ]
head(LD_data_na_all)

#LD_data_a_all
LD_data_a_all <- LD_data_a_all  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LD_data_a_all <- as.data.frame(LD_data_a_all)
LD_data_a_all <- LD_data_a_all[order(LD_data_a_all[,1],LD_data_a_all[,8]), ]
head(LD_data_a_all)

##Calculate the physical distance between each of the paired loci
#Add a column to each dataframe that is the physical distance between each of the loci in the pairs
LDdata_all$Distance <- abs(LDdata_all$BP_A - LDdata_all$BP_B)
LDdata_even$Distance <- abs(LDdata_even$BP_A - LDdata_even$BP_B)
LDdata_odd$Distance <- abs(LDdata_odd$BP_A - LDdata_odd$BP_B)
LDdata_even_ns$Distance <- abs(LDdata_even_ns$BP_A - LDdata_even_ns$BP_B)
LDdata_odd_ns$Distance <- abs(LDdata_odd_ns$BP_A - LDdata_odd_ns$BP_B)
LD_data_na_all$Distance <- abs(LD_data_na_all$BP_A - LD_data_na_all$BP_B)
LD_data_a_all$Distance <- abs(LD_data_a_all$BP_A - LD_data_a_all$BP_B)

# ################################################
# ##### TESTS TO SEE HOW TO USE THE FULL SET OF DATA- HOW TO FILTER R2? 
# LDdata_all_test <- LDdata_all  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=TRUE)))
# LDdata_all_test <- as.data.frame(LDdata_all_test)
# LDdata_all_test <- LDdata_all_test[order(LDdata_all_test[,1],LDdata_all_test[,8]), ]
# head(LDdata_all_test)
# tail(LDdata_all_test)
# 
# Subset <- LDdata_all_test[(LDdata_all_test$CHR_A == 1) & (LDdata_all_test$R2 >=.1),]
# dim(Subset)
# head(Subset)
# tail(Subset)
# 
# ggplot(data=Subset, aes(R2)) + geom_histogram(binwidth= .01)+ theme_bw() + ggtitle("Counts of the LD R2 values greater than 0.05 in All populations data set ")



##########################################################
#For running the plots without having to import all the original data 
#these datasets were exported at the end of the code originally

LDdata_all <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_all.txt",  header=TRUE, stringsAsFactors = FALSE, na.strings = "-" )
LDdata_even <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_even.txt",  header=TRUE, stringsAsFactors = FALSE, na.strings = "-" )
LDdata_odd <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_odd.txt",  header=TRUE, stringsAsFactors = FALSE, na.strings = "-" )
LDdata_even_ns <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_even_ns.txt",  header=TRUE, stringsAsFactors = FALSE, na.strings = "-" )
LDdata_odd_ns <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_odd_ns.txt",  header=TRUE, stringsAsFactors = FALSE, na.strings = "-" )
LD_data_na_all <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_na_all.txt",  header=TRUE, stringsAsFactors = FALSE, na.strings = "-" )
LD_data_a_all <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_a_all.txt",  header=TRUE, stringsAsFactors = FALSE, na.strings = "-" )


###################################################

#Facet wrap ggplot the rank of r2 per each chromosome for each of the 7 groups

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Ranked_LD_main_hierarchy.pdf", width = 9, height = 7)

# plot Individual Genotype Rate Post_Filtered_data, Ignore the Unassigned LD
ggplot(LDdata_all[(LDdata_all$CHR_A != "35") & (LDdata_all$R2 >=0.05),], aes(x=Rank, y= R2 )) + 
  geom_point(data= LDdata_all[(LDdata_all$CHR_A != "35") & (LDdata_all$R2 >=0.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value >0.05 by chromosome for ALL samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LDdata_even[(LDdata_even$CHR_A != "35") & (LDdata_even$R2 >=0.05),], aes(x=Rank, y= R2 )) + 
  geom_point(data= LDdata_even[(LDdata_even$CHR_A != "35") & (LDdata_even$R2 >=0.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value >0.05 by chromosome for Even lineage ")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LDdata_odd[(LDdata_odd$CHR_A != "35") & (LDdata_odd$R2 >=0.05),], aes(x=Rank, y= R2 )) + 
  geom_point(data= LDdata_odd[(LDdata_odd$CHR_A != "35") & (LDdata_odd$R2 >=0.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value >0.05 by chromosome for Odd lineage") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LDdata_even_ns[(LDdata_even_ns$CHR_A != "35") & (LDdata_even_ns$R2 >=0.05),], aes(x=Rank, y= R2 )) + 
  geom_point(data= LDdata_even_ns[(LDdata_even_ns$CHR_A != "35") & (LDdata_even_ns$R2 >=0.05),,], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value >0.05 by chromosome for Even No Susitna samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LDdata_odd_ns[(LDdata_odd_ns$CHR_A != "35") & (LDdata_odd_ns$R2 >=0.05),], aes(x=Rank, y= R2 )) + 
  geom_point(data= LDdata_odd_ns[(LDdata_odd_ns$CHR_A != "35") & (LDdata_odd_ns$R2 >=0.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value >0.05 by chromosome for Odd No Susitna samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LD_data_na_all[(LD_data_na_all$CHR_A != "35") & (LD_data_na_all$R2 >=0.05),], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_na_all[(LD_data_na_all$CHR_A != "35") & (LD_data_na_all$R2 >=0.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value >0.05 by chromosome for North American samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LD_data_a_all[(LD_data_a_all$CHR_A != "35") & (LD_data_a_all$R2 >=0.05),], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_a_all[(LD_data_a_all$CHR_A != "35") & (LD_data_a_all$R2 >=0.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value >0.05 by chromosome for Asian samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

dev.off()


#####################################################
# plot Individual Genotype Rate Post_Filtered_data, Ignore the Unassigned LD
###THESE have to be individually saved as jpegs because they have too many points to be rendered in pdf each time. 

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Ranked_LD_main_hierarchy_LDdata_all.jpg", ggplot(LDdata_all[LDdata_all$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LDdata_all[LDdata_all$CHR_A != "35",], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for ALL samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none"), width = 9,  height = 7,  dpi = 1200)

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Ranked_LD_main_hierarchy_LDdata_even.jpg", ggplot(LDdata_even[LDdata_even$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LDdata_even[LDdata_even$CHR_A != "35",], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for Even lineage ")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none"), width = 9,  height = 7,  dpi = 1200)

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Ranked_LD_main_hierarchy_LDdata_odd.jpg", ggplot(LDdata_odd[LDdata_odd$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LDdata_odd[LDdata_odd$CHR_A != "35",], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for Odd lineage") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none"), width = 9,  height = 7,  dpi = 1200)

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Ranked_LD_main_hierarchy_LDdata_even_ns.jpg", ggplot(LDdata_even_ns[LDdata_even_ns$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LDdata_even_ns[LDdata_even_ns$CHR_A != "35",], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for Even No Susitna samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none"), width = 9,  height = 7,  dpi = 1200)

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Ranked_LD_main_hierarchy_LDdata_odd_ns.jpg", ggplot(LDdata_odd_ns[LDdata_odd_ns$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LDdata_odd_ns[LDdata_odd_ns$CHR_A != "35",], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for Odd No Susitna samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none"), width = 9,  height = 7,  dpi = 1200)

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Ranked_LD_main_hierarchy_LD_data_na_all.jpg", ggplot(LD_data_na_all[LD_data_na_all$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_na_all[LD_data_na_all$CHR_A != "35",], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for North American samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none"), width = 9,  height = 7,  dpi = 1200)

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Ranked_LD_main_hierarchy_LD_data_a_all.jpg", ggplot(LD_data_a_all[LD_data_a_all$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_a_all[LD_data_a_all$CHR_A != "35",], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for Asian samples")+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none"), width = 9,  height = 7,  dpi = 1200)

##############################################
## Plot LD decay plots for each of the loci. 
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_Decay_by_CHR_main_hierarchy.pdf", width = 9, height = 7)

# plot Individual Genotype Rate Post_Filtered_data, Ignore the Unassigned LD
ggplot(LDdata_all[(LDdata_all$CHR_A != "35") & (LDdata_all$R2 >=.05),], aes(x=Distance, y= R2 )) + 
  geom_point(data= LDdata_all[(LDdata_all$CHR_A != "35") & (LDdata_all$R2 >=.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 value >0.05 by physical distance per chromosome for ALL samples")+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LDdata_even[(LDdata_even$CHR_A != "35") & (LDdata_even$R2 >=.05),], aes(x=Distance, y= R2 )) + 
  geom_point(data= LDdata_even[(LDdata_even$CHR_A != "35") & (LDdata_even$R2 >=.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 value >0.05 by physical distance per chromosome for Even lineage ")+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LDdata_odd[(LDdata_odd$CHR_A != "35") & (LDdata_odd$R2 >=.05),], aes(x=Distance, y= R2 )) + 
  geom_point(data= LDdata_odd[(LDdata_odd$CHR_A != "35") & (LDdata_odd$R2 >=.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 value >0.05 by physical distance per chromosome for Odd lineage")+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LDdata_even_ns[(LDdata_even_ns$CHR_A != "35") & (LDdata_even_ns$R2 >=.05),], aes(x=Distance, y= R2 )) + 
  geom_point(data= LDdata_even_ns[(LDdata_even_ns$CHR_A != "35") & (LDdata_even_ns$R2 >=.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 value >0.05 by physical distance per chromosome for Even No Susitna samples")+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LDdata_odd_ns[(LDdata_odd_ns$CHR_A != "35") & (LDdata_odd_ns$R2 >=.05),], aes(x=Distance, y= R2 )) + 
  geom_point(data= LDdata_odd_ns[(LDdata_odd_ns$CHR_A != "35") & (LDdata_odd_ns$R2 >=.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 value >0.05 by physical distance per chromosome for Odd No Susitna samples")+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LD_data_na_all[(LD_data_na_all$CHR_A != "35") & (LD_data_na_all$R2 >=.05),], aes(x=Distance, y= R2 )) + 
  geom_point(data= LD_data_na_all[(LD_data_na_all$CHR_A != "35") & (LD_data_na_all$R2 >=.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 value >0.05 by physical distance per chromosome for North American samples")+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

ggplot(LD_data_a_all[(LD_data_a_all$CHR_A != "35") & (LD_data_a_all$R2 >=.05),], aes(x=Distance, y= R2 )) + 
  geom_point(data= LD_data_a_all[(LD_data_a_all$CHR_A != "35") & (LD_data_a_all$R2 >=.05),], colour = "black", size =.5 ) + facet_wrap(~CHR_A, scale="free") + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 value >0.05 by physical distance per chromosome for Asian samples")+theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

dev.off()

#clear out any memory that we dont need after making plots
#rm( LDdata_all, LDdata_even, LDdata_odd, LDdata_even_ns, LDdata_odd_ns, LD_data_na_all, LD_data_a_all)

gc()


#save the dataframes that were made here so that they can be imported separately. 

####Write the combined data file for each LD R2 group and a filtered version that is easier for plotting LD  to files to be used in other R scripts
##LDdata_all
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_all.txt", "wb")
write.table(LDdata_all,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

LDdata_all_.2_BPA_BPB_R2 <- LDdata_all[LDdata_all$R2 >=0.2, c("CHR_A","BP_A", "BP_B", "R2")]
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_all_.2_BPA_BPB_R2.txt", "wb")
write.table(LDdata_all_.2_BPA_BPB_R2,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

rm(LDdata_all, LDdata_all_.2_BPA_BPB_R2)
gc()

##LDdata_even
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_even.txt", "wb")
write.table(LDdata_even,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

LDdata_even_.2_BPA_BPB_R2 <- LDdata_even[LDdata_even$R2 >=0.2, c("CHR_A","BP_A", "BP_B", "R2")]
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_even_.2_BPA_BPB_R2.txt", "wb")
write.table(LDdata_even_.2_BPA_BPB_R2,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

rm(LDdata_even, LDdata_even_.2_BPA_BPB_R2)
gc()

#LDdata_odd
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_odd.txt", "wb")
write.table(LDdata_odd,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

LDdata_odd_.2_BPA_BPB_R2 <- LDdata_odd[LDdata_odd$R2 >=0.2, c("CHR_A","BP_A", "BP_B", "R2")]
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_odd_.2_BPA_BPB_R2.txt", "wb")
write.table(LDdata_odd_.2_BPA_BPB_R2,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

rm(LDdata_odd, LDdata_odd_.2_BPA_BPB_R2)
gc()

#LDdata_even_ns
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_even_ns.txt", "wb")
write.table(LDdata_even_ns,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

LDdata_even_ns_.2_BPA_BPB_R2 <- LDdata_even_ns[LDdata_even_ns$R2 >=0.2, c("CHR_A","BP_A", "BP_B", "R2")]
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_even_ns_.2_BPA_BPB_R2.txt", "wb")
write.table(LDdata_even_ns_.2_BPA_BPB_R2,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

rm(LDdata_even_ns, LDdata_even_ns_.2_BPA_BPB_R2)
gc()

#LDdata_odd_ns
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_odd_ns.txt", "wb")
write.table(LDdata_odd_ns,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

LDdata_odd_ns_.2_BPA_BPB_R2 <- LDdata_odd_ns[LDdata_odd_ns$R2 >=0.2, c("CHR_A","BP_A", "BP_B", "R2")]
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_odd_ns_.2_BPA_BPB_R2.txt", "wb")
write.table(LDdata_odd_ns_.2_BPA_BPB_R2,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

rm(LDdata_odd_ns, LDdata_odd_ns_.2_BPA_BPB_R2)
gc()

#LD_data_na_all
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_na_all.txt", "wb")
write.table(LD_data_na_all,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

LD_data_na_all_.2_BPA_BPB_R2 <- LD_data_na_all[LD_data_na_all$R2 >=0.2, c("CHR_A","BP_A", "BP_B", "R2")]
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_na_all_.2_BPA_BPB_R2.txt", "wb")
write.table(LD_data_na_all_.2_BPA_BPB_R2,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

rm(LD_data_na_all, LD_data_na_all_.2_BPA_BPB_R2)
gc()

#LD_data_a_all
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_a_all.txt", "wb")
write.table(LD_data_a_all,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

LD_data_a_all_.2_BPA_BPB_R2 <- LD_data_a_all[LD_data_a_all$R2 >=0.2, c("CHR_A","BP_A", "BP_B", "R2")]
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_a_all_.2_BPA_BPB_R2.txt", "wb")
write.table(LD_data_a_all_.2_BPA_BPB_R2,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

rm(LD_data_a_all, LD_data_a_all_.2_BPA_BPB_R2)
gc()

#save this workspace independently !!


