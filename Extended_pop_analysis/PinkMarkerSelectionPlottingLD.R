### Plotting the LD r2 from PLINK 
###    for the pink panel of 2312 markers
###    
### Garrett McKinney and Carolyn Tarpey | May 2018
### ---------------------------------------

#This R code was originally written by Garrett McKinney to plot the output of PLINK LD r2 for each chromosome. 
#It requires the LD output from PLINK. We used the alignment of the pinks to Chinook to get the Chromosome 
#and position assignments here, so there are 34 Chromosomes, and 


#Pink data filtering
library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)

#load genepop files
LDdata<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/Panel_2312_PLINK_out.ld",sep="")
LDdata[1:5,1:5]
dim(LDdata)

LDdata_melt <- melt(LDdata)


# #This is the example code that garrett Sent me:
# ggplot()+geom_point(data=LDdata,aes(x=BP_A,y=BP_B,color=R2),shape=15)+scale_color_gradient(low="white",high="dark red")+theme_bw()
# ggplot()+geom_point(data=LDdata[LDdata$R2>=0.3,],aes(x=BP_A,y=BP_B,color=R2),shape=15)+scale_color_gradient(low="black",high="red")+theme_bw()


#Parse the data by chromosome- though we ran all chromosomes together in PLINK, we used the command inter-chr to 
# ensure that the LD was calculated only within a chromosome. 

