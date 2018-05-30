### Plotting the LD r2 from PLINK 
###    for the pink panel of 2312 markers
###    THE CODE FOR THE PLOTS OF THE ACTUAL LD /CHROM IS IN _TWO
### Garrett McKinney and Carolyn Tarpey | May 2018
### ---------------------------------------

#This R code was originally written by Garrett McKinney to plot the output of PLINK LD r2 for each chromosome. 
#It requires the LD output from PLINK. We used the alignment of the pinks to Chinook to get the Chromosome 
#and position assignments here, so there are 34 Chromosomes, and chromosome 35 are those that aligned to the unassigned 

#TO MAKE THE PLOTS OF THE LD PER CHROMOSOME, USE THIS SCRIPT TO IMPORT THE DATA AND THE PinkMarkerSelectionPlottingLD_TWO.R TO MAKE THE ACTUAL PLOTS

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
#### 2312 ALL 
LDdata_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr1.ld",sep="")
LDdata_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr2.ld",sep="")
LDdata_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr3.ld",sep="")
LDdata_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr4.ld",sep="")
LDdata_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr5.ld",sep="")
LDdata_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr6.ld",sep="")
LDdata_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr7.ld",sep="")
LDdata_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr8.ld",sep="")
LDdata_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr9.ld",sep="")
LDdata_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr10.ld",sep="")
LDdata_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr11.ld",sep="")
LDdata_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr12.ld",sep="")
LDdata_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr13.ld",sep="")
LDdata_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr14.ld",sep="")
LDdata_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr15.ld",sep="")
LDdata_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr16.ld",sep="")
LDdata_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr17.ld",sep="")
LDdata_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr18.ld",sep="")
LDdata_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr19.ld",sep="")
LDdata_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr20.ld",sep="")
LDdata_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr21.ld",sep="")
LDdata_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr22.ld",sep="")
LDdata_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr23.ld",sep="")
LDdata_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr24.ld",sep="")
LDdata_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr25.ld",sep="")
LDdata_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr26.ld",sep="")
LDdata_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr27.ld",sep="")
LDdata_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr28.ld",sep="")
LDdata_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr29.ld",sep="")
LDdata_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr30.ld",sep="")
LDdata_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr31.ld",sep="")
LDdata_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr32.ld",sep="")
LDdata_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr33.ld",sep="")
LDdata_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr34.ld",sep="")
LDdata_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_ALL_out/Panel_2312_PLINK_out_chr35.ld",sep="")


#### 2312 Even
#load PLINK LD output files
LD_even_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr1.ld",sep="")
LD_even_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr2.ld",sep="")
LD_even_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr3.ld",sep="")
LD_even_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr4.ld",sep="")
LD_even_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr5.ld",sep="")
LD_even_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr6.ld",sep="")
LD_even_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr7.ld",sep="")
LD_even_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr8.ld",sep="")
LD_even_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr9.ld",sep="")
LD_even_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr10.ld",sep="")
LD_even_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr11.ld",sep="")
LD_even_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr12.ld",sep="")
LD_even_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr13.ld",sep="")
LD_even_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr14.ld",sep="")
LD_even_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr15.ld",sep="")
LD_even_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr16.ld",sep="")
LD_even_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr17.ld",sep="")
LD_even_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr18.ld",sep="")
LD_even_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr19.ld",sep="")
LD_even_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr20.ld",sep="")
LD_even_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr21.ld",sep="")
LD_even_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr22.ld",sep="")
LD_even_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr23.ld",sep="")
LD_even_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr24.ld",sep="")
LD_even_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr25.ld",sep="")
LD_even_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr26.ld",sep="")
LD_even_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr27.ld",sep="")
LD_even_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr28.ld",sep="")
LD_even_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr29.ld",sep="")
LD_even_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr30.ld",sep="")
LD_even_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr31.ld",sep="")
LD_even_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr32.ld",sep="")
LD_even_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr33.ld",sep="")
LD_even_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr34.ld",sep="")
LD_even_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_out/Panel_2312_Even_out_chr35.ld",sep="")

#### 2312 Odd
#load PLINK LD output files
LD_odd_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr1.ld",sep="")
LD_odd_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr2.ld",sep="")
LD_odd_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr3.ld",sep="")
LD_odd_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr4.ld",sep="")
LD_odd_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr5.ld",sep="")
LD_odd_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr6.ld",sep="")
LD_odd_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr7.ld",sep="")
LD_odd_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr8.ld",sep="")
LD_odd_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr9.ld",sep="")
LD_odd_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr10.ld",sep="")
LD_odd_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr11.ld",sep="")
LD_odd_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr12.ld",sep="")
LD_odd_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr13.ld",sep="")
LD_odd_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr14.ld",sep="")
LD_odd_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr15.ld",sep="")
LD_odd_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr16.ld",sep="")
LD_odd_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr17.ld",sep="")
LD_odd_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr18.ld",sep="")
LD_odd_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr19.ld",sep="")
LD_odd_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr20.ld",sep="")
LD_odd_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr21.ld",sep="")
LD_odd_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr22.ld",sep="")
LD_odd_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr23.ld",sep="")
LD_odd_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr24.ld",sep="")
LD_odd_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr25.ld",sep="")
LD_odd_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr26.ld",sep="")
LD_odd_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr27.ld",sep="")
LD_odd_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr28.ld",sep="")
LD_odd_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr29.ld",sep="")
LD_odd_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr30.ld",sep="")
LD_odd_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr31.ld",sep="")
LD_odd_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr32.ld",sep="")
LD_odd_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr33.ld",sep="")
LD_odd_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr34.ld",sep="")
LD_odd_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_out/Panel_2312_Odd_out_chr35.ld",sep="")

#### 2312 Even Asia
#load PLINK LD output files
LD_even_A_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr1.ld",sep="")
LD_even_A_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr2.ld",sep="")
LD_even_A_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr3.ld",sep="")
LD_even_A_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr4.ld",sep="")
LD_even_A_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr5.ld",sep="")
LD_even_A_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr6.ld",sep="")
LD_even_A_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr7.ld",sep="")
LD_even_A_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr8.ld",sep="")
LD_even_A_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr9.ld",sep="")
LD_even_A_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr10.ld",sep="")
LD_even_A_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr11.ld",sep="")
LD_even_A_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr12.ld",sep="")
LD_even_A_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr13.ld",sep="")
LD_even_A_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr14.ld",sep="")
LD_even_A_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr15.ld",sep="")
LD_even_A_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr16.ld",sep="")
LD_even_A_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr17.ld",sep="")
LD_even_A_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr18.ld",sep="")
LD_even_A_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr19.ld",sep="")
LD_even_A_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr20.ld",sep="")
LD_even_A_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr21.ld",sep="")
LD_even_A_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr22.ld",sep="")
LD_even_A_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr23.ld",sep="")
LD_even_A_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr24.ld",sep="")
LD_even_A_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr25.ld",sep="")
LD_even_A_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr26.ld",sep="")
LD_even_A_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr27.ld",sep="")
LD_even_A_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr28.ld",sep="")
LD_even_A_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr29.ld",sep="")
LD_even_A_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr30.ld",sep="")
LD_even_A_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr31.ld",sep="")
LD_even_A_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr32.ld",sep="")
LD_even_A_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr33.ld",sep="")
LD_even_A_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr34.ld",sep="")
LD_even_A_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_A_out/Panel_2312_Even_A_out_chr35.ld",sep="")

#### 2312 Even North America
#load PLINK LD output files
LD_even_NA_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr1.ld",sep="")
LD_even_NA_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr2.ld",sep="")
LD_even_NA_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr3.ld",sep="")
LD_even_NA_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr4.ld",sep="")
LD_even_NA_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr5.ld",sep="")
LD_even_NA_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr6.ld",sep="")
LD_even_NA_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr7.ld",sep="")
LD_even_NA_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr8.ld",sep="")
LD_even_NA_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr9.ld",sep="")
LD_even_NA_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr10.ld",sep="")
LD_even_NA_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr11.ld",sep="")
LD_even_NA_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr12.ld",sep="")
LD_even_NA_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr13.ld",sep="")
LD_even_NA_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr14.ld",sep="")
LD_even_NA_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr15.ld",sep="")
LD_even_NA_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr16.ld",sep="")
LD_even_NA_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr17.ld",sep="")
LD_even_NA_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr18.ld",sep="")
LD_even_NA_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr19.ld",sep="")
LD_even_NA_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr20.ld",sep="")
LD_even_NA_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr21.ld",sep="")
LD_even_NA_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr22.ld",sep="")
LD_even_NA_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr23.ld",sep="")
LD_even_NA_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr24.ld",sep="")
LD_even_NA_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr25.ld",sep="")
LD_even_NA_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr26.ld",sep="")
LD_even_NA_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr27.ld",sep="")
LD_even_NA_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr28.ld",sep="")
LD_even_NA_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr29.ld",sep="")
LD_even_NA_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr30.ld",sep="")
LD_even_NA_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr31.ld",sep="")
LD_even_NA_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr32.ld",sep="")
LD_even_NA_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr33.ld",sep="")
LD_even_NA_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr34.ld",sep="")
LD_even_NA_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Even_NA_out/Panel_2312_Even_NA_out_chr35.ld",sep="")

#### 2312 odd Asia
#load PLINK LD output files
LD_odd_a_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr1.ld",sep="")
LD_odd_a_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr2.ld",sep="")
LD_odd_a_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr3.ld",sep="")
LD_odd_a_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr4.ld",sep="")
LD_odd_a_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr5.ld",sep="")
LD_odd_a_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr6.ld",sep="")
LD_odd_a_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr7.ld",sep="")
LD_odd_a_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr8.ld",sep="")
LD_odd_a_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr9.ld",sep="")
LD_odd_a_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr10.ld",sep="")
LD_odd_a_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr11.ld",sep="")
LD_odd_a_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr12.ld",sep="")
LD_odd_a_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr13.ld",sep="")
LD_odd_a_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr14.ld",sep="")
LD_odd_a_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr15.ld",sep="")
LD_odd_a_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr16.ld",sep="")
LD_odd_a_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr17.ld",sep="")
LD_odd_a_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr18.ld",sep="")
LD_odd_a_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr19.ld",sep="")
LD_odd_a_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr20.ld",sep="")
LD_odd_a_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr21.ld",sep="")
LD_odd_a_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr22.ld",sep="")
LD_odd_a_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr23.ld",sep="")
LD_odd_a_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr24.ld",sep="")
LD_odd_a_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr25.ld",sep="")
LD_odd_a_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr26.ld",sep="")
LD_odd_a_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr27.ld",sep="")
LD_odd_a_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr28.ld",sep="")
LD_odd_a_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr29.ld",sep="")
LD_odd_a_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr30.ld",sep="")
LD_odd_a_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr31.ld",sep="")
LD_odd_a_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr32.ld",sep="")
LD_odd_a_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr33.ld",sep="")
LD_odd_a_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr34.ld",sep="")
LD_odd_a_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_A_out/Panel_2312_Odd_A_out_chr35.ld",sep="")

#### 2312 Odd North America
#load PLINK LD output files
LD_odd_na_chr1 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr1.ld",sep="")
LD_odd_na_chr2 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr2.ld",sep="")
LD_odd_na_chr3 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr3.ld",sep="")
LD_odd_na_chr4 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr4.ld",sep="")
LD_odd_na_chr5 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr5.ld",sep="")
LD_odd_na_chr6 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr6.ld",sep="")
LD_odd_na_chr7 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr7.ld",sep="")
LD_odd_na_chr8 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr8.ld",sep="")
LD_odd_na_chr9 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr9.ld",sep="")
LD_odd_na_chr10<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr10.ld",sep="")
LD_odd_na_chr11<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr11.ld",sep="")
LD_odd_na_chr12<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr12.ld",sep="")
LD_odd_na_chr13<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr13.ld",sep="")
LD_odd_na_chr14<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr14.ld",sep="")
LD_odd_na_chr15<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr15.ld",sep="")
LD_odd_na_chr16<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr16.ld",sep="")
LD_odd_na_chr17<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr17.ld",sep="")
LD_odd_na_chr18<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr18.ld",sep="")
LD_odd_na_chr19<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr19.ld",sep="")
LD_odd_na_chr20<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr20.ld",sep="")
LD_odd_na_chr21<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr21.ld",sep="")
LD_odd_na_chr22<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr22.ld",sep="")
LD_odd_na_chr23<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr23.ld",sep="")
LD_odd_na_chr24<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr24.ld",sep="")
LD_odd_na_chr25<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr25.ld",sep="")
LD_odd_na_chr26<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr26.ld",sep="")
LD_odd_na_chr27<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr27.ld",sep="")
LD_odd_na_chr28<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr28.ld",sep="")
LD_odd_na_chr29<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr29.ld",sep="")
LD_odd_na_chr30<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr30.ld",sep="")
LD_odd_na_chr31<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr31.ld",sep="")
LD_odd_na_chr32<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr32.ld",sep="")
LD_odd_na_chr33<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr33.ld",sep="")
LD_odd_na_chr34<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr34.ld",sep="")
LD_odd_na_chr0 <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/2312_Odd_NA_out/Panel_2312_Odd_NA_out_chr35.ld",sep="")


#combine each of the chromosome specific LDs into one big table:
LDdata_all <- rbind( LDdata_chr1,  LDdata_chr2, LDdata_chr3, LDdata_chr4,LDdata_chr5,LDdata_chr6,  LDdata_chr7,  LDdata_chr8,  LDdata_chr9,
                     LDdata_chr10, LDdata_chr11,  LDdata_chr12,  LDdata_chr13, LDdata_chr14,  LDdata_chr15,  LDdata_chr16,  LDdata_chr17,  LDdata_chr18,  LDdata_chr19,
                     LDdata_chr20,  LDdata_chr21,  LDdata_chr22, LDdata_chr23,  LDdata_chr24, LDdata_chr25,  LDdata_chr26,  LDdata_chr27,  LDdata_chr28,  LDdata_chr29,  LDdata_chr30,  LDdata_chr31,  LDdata_chr32,  LDdata_chr33, LDdata_chr34, LDdata_chr0)

LDdata_even <- rbind( LD_even_chr1,  LD_even_chr2, LD_even_chr3, LD_even_chr4, LD_even_chr5, LD_even_chr6,  LD_even_chr7, LD_even_chr8, LD_even_chr9,
                      LD_even_chr10, LD_even_chr11,  LD_even_chr12,  LD_even_chr13, LD_even_chr14,  LD_even_chr15,  LD_even_chr16,  LD_even_chr17,  LD_even_chr18,  LD_even_chr19,
                      LD_even_chr20,  LD_even_chr21,  LD_even_chr22, LD_even_chr23,  LD_even_chr24, LD_even_chr25,  LD_even_chr26,  LD_even_chr27,  LD_even_chr28,  LD_even_chr29,  LD_even_chr30,  LD_even_chr31,  LD_even_chr32,  LD_even_chr33, LD_even_chr34, LD_even_chr0)

LDdata_odd <- rbind( LD_odd_chr1,  LD_odd_chr2, LD_odd_chr3, LD_odd_chr4,LD_odd_chr5,LD_odd_chr6,  LD_odd_chr7,  LD_odd_chr8,  LD_odd_chr9,
                     LD_odd_chr10, LD_odd_chr11,  LD_odd_chr12,  LD_odd_chr13, LD_odd_chr14,  LD_odd_chr15,  LD_odd_chr16,  LD_odd_chr17,  LD_odd_chr18,  LD_odd_chr19,
                     LD_odd_chr20,  LD_odd_chr21,  LD_odd_chr22, LD_odd_chr23,  LD_odd_chr24, LD_odd_chr25,  LD_odd_chr26,  LD_odd_chr27,  LD_odd_chr28,  LD_odd_chr29,  LD_odd_chr30,  LD_odd_chr31,  LD_odd_chr32,  LD_odd_chr33, LD_odd_chr34, LD_odd_chr0)

LD_data_even_A <- rbind( LD_even_A_chr1,  LD_even_A_chr2, LD_even_A_chr3, LD_even_A_chr4, LD_even_A_chr5, LD_even_A_chr6,  LD_even_A_chr7, LD_even_A_chr8, LD_even_A_chr9,
                         LD_even_A_chr10, LD_even_A_chr11,  LD_even_A_chr12,  LD_even_A_chr13, LD_even_A_chr14,  LD_even_A_chr15,  LD_even_A_chr16,  LD_even_A_chr17,  LD_even_A_chr18,  LD_even_A_chr19,
                         LD_even_A_chr20,  LD_even_A_chr21,  LD_even_A_chr22, LD_even_A_chr23,  LD_even_A_chr24, LD_even_A_chr25,  LD_even_A_chr26,  LD_even_A_chr27,  LD_even_A_chr28,  LD_even_A_chr29,  LD_even_A_chr30,  LD_even_A_chr31,  LD_even_A_chr32,  LD_even_A_chr33, LD_even_A_chr34, LD_even_A_chr0)

LD_data_even_NA <- rbind( LD_even_NA_chr1,  LD_even_NA_chr2, LD_even_NA_chr3, LD_even_NA_chr4, LD_even_NA_chr5, LD_even_NA_chr6,  LD_even_NA_chr7, LD_even_NA_chr8, LD_even_NA_chr9,
                          LD_even_NA_chr10, LD_even_NA_chr11,  LD_even_NA_chr12,  LD_even_NA_chr13, LD_even_NA_chr14,  LD_even_NA_chr15,  LD_even_NA_chr16,  LD_even_NA_chr17,  LD_even_NA_chr18,  LD_even_NA_chr19,
                          LD_even_NA_chr20,  LD_even_NA_chr21,  LD_even_NA_chr22, LD_even_NA_chr23,  LD_even_NA_chr24, LD_even_NA_chr25,  LD_even_NA_chr26,  LD_even_NA_chr27,  LD_even_NA_chr28,  LD_even_NA_chr29,  LD_even_NA_chr30,  LD_even_NA_chr31,  LD_even_NA_chr32,  LD_even_NA_chr33, LD_even_NA_chr34, LD_even_NA_chr0)

LD_data_odd_a <- rbind( LD_odd_a_chr1,  LD_odd_a_chr2, LD_odd_a_chr3, LD_odd_a_chr4, LD_odd_a_chr5, LD_odd_a_chr6,  LD_odd_a_chr7, LD_odd_a_chr8, LD_odd_a_chr9,
                        LD_odd_a_chr10, LD_odd_a_chr11,  LD_odd_a_chr12,  LD_odd_a_chr13, LD_odd_a_chr14,  LD_odd_a_chr15,  LD_odd_a_chr16,  LD_odd_a_chr17,  LD_odd_a_chr18,  LD_odd_a_chr19,
                        LD_odd_a_chr20,  LD_odd_a_chr21,  LD_odd_a_chr22, LD_odd_a_chr23,  LD_odd_a_chr24, LD_odd_a_chr25,  LD_odd_a_chr26,  LD_odd_a_chr27,  LD_odd_a_chr28,  LD_odd_a_chr29,  LD_odd_a_chr30,  LD_odd_a_chr31,  LD_odd_a_chr32,  LD_odd_a_chr33, LD_odd_a_chr34, LD_odd_a_chr0)

LD_data_odd_na <- rbind( LD_odd_na_chr1,  LD_odd_na_chr2, LD_odd_na_chr3, LD_odd_na_chr4, LD_odd_na_chr5, LD_odd_na_chr6,  LD_odd_na_chr7, LD_odd_na_chr8, LD_odd_na_chr9,
                       LD_odd_na_chr10, LD_odd_na_chr11,  LD_odd_na_chr12,  LD_odd_na_chr13, LD_odd_na_chr14,  LD_odd_na_chr15,  LD_odd_na_chr16,  LD_odd_na_chr17,  LD_odd_na_chr18,  LD_odd_na_chr19,
                       LD_odd_na_chr20,  LD_odd_na_chr21,  LD_odd_na_chr22, LD_odd_na_chr23,  LD_odd_na_chr24, LD_odd_na_chr25,  LD_odd_na_chr26,  LD_odd_na_chr27,  LD_odd_na_chr28,  LD_odd_na_chr29,  LD_odd_na_chr30,  LD_odd_na_chr31,  LD_odd_na_chr32,  LD_odd_na_chr33, LD_odd_na_chr34, LD_odd_na_chr0)


#### 
###
#
#Calculate the rank of each of the LDs per chromosome for all 7 groups. 

#Add a column to each dataframe that is for the rank of the R2 per chromosome
LDdata_all$Rank <- NA
LDdata_even$Rank <-NA
LDdata_odd$Rank <- NA
LD_data_even_A$Rank <- NA
LD_data_even_NA$Rank <-NA
LD_data_odd_a$Rank <- NA
LD_data_odd_na$Rank <- NA

#this uses dyplr to rank the r2 per each chromosome, then orders the rank within each chromosome
LDdata_all <- LDdata_all  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LDdata_all <- as.data.frame(LDdata_all)
LDdata_all_sort <- LDdata_all[order(LDdata_all[,1],LDdata_all[,8]), ]
head(LDdata_all_sort)

LDdata_even 
LDdata_even <- LDdata_even  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LDdata_even <- as.data.frame(LDdata_even)
LDdata_even_sort <- LDdata_even[order(LDdata_even[,1],LDdata_even[,8]), ]
head(LDdata_even_sort)

LDdata_odd 
LDdata_odd <- LDdata_odd  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LDdata_odd <- as.data.frame(LDdata_odd)
LDdata_odd_sort <- LDdata_odd[order(LDdata_odd[,1],LDdata_odd[,8]), ]
head(LDdata_odd_sort)

LD_data_even_A 
LD_data_even_A <- LD_data_even_A  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LD_data_even_A <- as.data.frame(LD_data_even_A)
LD_data_even_A_sort <- LD_data_even_A[order(LD_data_even_A[,1],LD_data_even_A[,8]), ]
head(LD_data_even_A_sort)

LD_data_even_NA
LD_data_even_NA <- LD_data_even_NA  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LD_data_even_NA <- as.data.frame(LD_data_even_NA)
LD_data_even_NA_sort <- LD_data_even_NA[order(LD_data_even_NA[,1],LD_data_even_NA[,8]), ]
head(LD_data_even_NA_sort)

LD_data_odd_a 
LD_data_odd_a <- LD_data_odd_a  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LD_data_odd_a <- as.data.frame(LD_data_odd_a)
LD_data_odd_a_sort <- LD_data_odd_a[order(LD_data_odd_a[,1],LD_data_odd_a[,8]), ]
head(LD_data_odd_a_sort)

LD_data_odd_na
LD_data_odd_na <- LD_data_odd_na  %>% group_by(CHR_A) %>% mutate(Rank = order(order(R2, decreasing=FALSE)))
LD_data_odd_na <- as.data.frame(LD_data_odd_na)
LD_data_odd_na_sort <- LD_data_odd_na[order(LD_data_odd_na[,1],LD_data_odd_na[,8]), ]
head(LD_data_odd_na_sort)


#Facet wrap ggplot the rank of r2 per each chromosome for each of the 6 groups
#Color scheme to use
Blue_pink_7 <- c("All"= "#8800C3", "Even"="#000080", "Odd"="#FF00FF", "Even_NA"="#AAAAD4", "Even_A"= "#5555AA", "Odd_NA"="#FFAAFF", "Odd_A" = "#FF55FF")
plot(rep(1,7),col=Blue_pink_7,pch=19,cex=7)

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/Ranked_LD_all_groups.pdf", width = 9, height = 7)

# plot Individual Genotype Rate Post_Filtered_data, Ignore the Unassigned LD
ggplot(LDdata_all_sort[LDdata_all_sort$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LDdata_all_sort[LDdata_all_sort$CHR_A != "35",], colour = "#8800C3") + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for ALL samples")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(LDdata_even_sort[LDdata_even_sort$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LDdata_even_sort[LDdata_even_sort$CHR_A != "35",], colour = "#000080") + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for Even lineage ")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(LDdata_odd_sort[LDdata_odd_sort$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LDdata_odd_sort[LDdata_odd_sort$CHR_A != "35",], colour = "#FF00FF") + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for Odd lineage")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(LD_data_even_A_sort[LD_data_even_A_sort$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_even_A_sort[LD_data_even_A_sort$CHR_A != "35",], colour = "#000080") + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for Even Asian samples")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(LD_data_even_NA_sort[LD_data_even_NA_sort$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_even_NA_sort[LD_data_even_NA_sort$CHR_A != "35",], colour = "#000080") + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for Even North American samples")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(LD_data_odd_a_sort[LD_data_odd_a_sort$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_odd_a_sort[LD_data_odd_a_sort$CHR_A != "35",], colour = "#FF00FF") + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for Odd Asian samples")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(LD_data_odd_na_sort[LD_data_odd_na_sort$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_odd_na_sort[LD_data_odd_na_sort$CHR_A != "35",], colour = "#FF00FF") + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for Odd North American samples")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()

##### SMALLER POINT SIZE EASIER TO SEE IN PDF

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/Ranked_LD_all_groups_size.pdf", width = 9, height = 7)

# plot Individual Genotype Rate Post_Filtered_data, Ignore the Unassigned LD
ggplot(LDdata_all_sort[LDdata_all_sort$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LDdata_all_sort[LDdata_all_sort$CHR_A != "35",], colour = "#8800C3", size =.5 ) + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for ALL samples")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(LDdata_even_sort[LDdata_even_sort$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LDdata_even_sort[LDdata_even_sort$CHR_A != "35",], colour = "#000080", size =.5 ) + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for Even lineage ")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(LDdata_odd_sort[LDdata_odd_sort$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LDdata_odd_sort[LDdata_odd_sort$CHR_A != "35",], colour = "#FF00FF", size =.5 ) + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for Odd lineage")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(LD_data_even_A_sort[LD_data_even_A_sort$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_even_A_sort[LD_data_even_A_sort$CHR_A != "35",], colour = "#000080", size =.5 ) + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for Even Asian samples")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(LD_data_even_NA_sort[LD_data_even_NA_sort$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_even_NA_sort[LD_data_even_NA_sort$CHR_A != "35",], colour = "#000080", size =.5 ) + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for Even North American samples")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(LD_data_odd_a_sort[LD_data_odd_a_sort$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_odd_a_sort[LD_data_odd_a_sort$CHR_A != "35",], colour = "#FF00FF", size =.5 ) + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for Odd Asian samples")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(LD_data_odd_na_sort[LD_data_odd_na_sort$CHR_A != "35",], aes(x=Rank, y= R2 )) + 
  geom_point(data= LD_data_odd_na_sort[LD_data_odd_na_sort$CHR_A != "35",], colour = "#FF00FF", size =.5 ) + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("Ranked LD R2 value by chromosome for Odd North American samples")+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

dev.off()


#### 
###
#
#Calculate the physical distance between each of the paired loci
#Add a column to each dataframe that is the physical distance between each of the loci in the pairs
LDdata_all$Distance <- abs(LDdata_all$BP_A - LDdata_all$BP_B)
LDdata_even$Distance <- abs(LDdata_even$BP_A - LDdata_even$BP_B)
LDdata_odd$Distance <- abs(LDdata_odd$BP_A - LDdata_odd$BP_B)
LD_data_even_A$Distance <-  abs(LD_data_even_A$BP_A - LD_data_even_A$BP_B)
LD_data_even_NA$Distance <- abs(LD_data_even_NA$BP_A - LD_data_even_NA$BP_B)
LD_data_odd_a$Distance <-  abs(LD_data_odd_a$BP_A - LD_data_odd_a$BP_B)
LD_data_odd_na$Distance <-  abs(LD_data_odd_na$BP_A - LD_data_odd_na$BP_B)


## Plot LD decay plots for each of the loci. 

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/LD_Decay_by_CHR.pdf", width = 9, height = 7)

# plot Individual Genotype Rate Post_Filtered_data, Ignore the Unassigned LD
ggplot(LDdata_all[LDdata_all$CHR_A != "35",], aes(x=Distance, y= R2 )) + 
  geom_point(data= LDdata_all[LDdata_all$CHR_A != "35",], colour = "#8800C3", size =.5 ) + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 valaue by physical distance per chromosome for ALL samples")+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_continuous(labels = comma)

ggplot(LDdata_even[LDdata_even$CHR_A != "35",], aes(x=Distance, y= R2 )) + 
  geom_point(data= LDdata_even[LDdata_even$CHR_A != "35",], colour = "#000080", size =.5 ) + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 valaue by physical distance per chromosome for Even lineage ")+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_continuous(labels = comma)

ggplot(LDdata_odd[LDdata_odd$CHR_A != "35",], aes(x=Distance, y= R2 )) + 
  geom_point(data= LDdata_odd[LDdata_odd$CHR_A != "35",], colour = "#FF00FF", size =.5 ) + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 valaue by physical distance per chromosome for Odd lineage")+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_continuous(labels = comma)

ggplot(LD_data_even_A[LD_data_even_A$CHR_A != "35",], aes(x=Distance, y= R2 )) + 
  geom_point(data= LD_data_even_A[LD_data_even_A$CHR_A != "35",], colour = "#000080", size =.5 ) + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 valaue by physical distance per chromosome for Even Asian samples")+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_continuous(labels = comma)

ggplot(LD_data_even_NA[LD_data_even_NA$CHR_A != "35",], aes(x=Distance, y= R2 )) + 
  geom_point(data= LD_data_even_NA[LD_data_even_NA$CHR_A != "35",], colour = "#000080", size =.5 ) + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 valaue by physical distance per chromosome for Even North American samples")+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_continuous(labels = comma)

ggplot(LD_data_odd_a[LD_data_odd_a$CHR_A != "35",], aes(x=Distance, y= R2 )) + 
  geom_point(data= LD_data_odd_a[LD_data_odd_a$CHR_A != "35",], colour = "#FF00FF", size =.5 ) + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 valaue by physical distance per chromosome for Odd Asian samples")+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_continuous(labels = comma)

ggplot(LD_data_odd_na[LD_data_odd_na$CHR_A != "35",], aes(x=Distance, y= R2 )) + 
  geom_point(data= LD_data_odd_na[LD_data_odd_na$CHR_A != "35",], colour = "#FF00FF", size =.5 ) + facet_wrap(~CHR_A) + theme_bw() + guides(fill= FALSE) +
  ggtitle("LD R2 valaue by physical distance per chromosome for Odd North American samples")+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_x_continuous(labels = comma)

dev.off()







