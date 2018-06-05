### Plotting the LD r2 from PLINK
###     Main Hierarchy
###    for the full dataset of One SNP per Tag- 23759 loci. 
###    THIS IS THE SECOND SHEET FOR THE ACTUAL PLOTS OF THE CHROMOSOMES
### Garrett McKinney and Carolyn Tarpey | May 2018
## ---------------------------------------

## THIS IS THE SECOND SHEET- YOU NEED THE MAIN SHEET TO IMPORT THE DATA FOR EACH OF THESE PLOTS

#This R code was originally written by Garrett McKinney to plot the output of PLINK LD r2 for each chromosome. 
#It requires the LD output from PLINK. We used the alignment of the pinks to Chinook to get the Chromosome 
#and position assignments here, so there are 34 Chromosomes, and chromosome 35 are those that aligned to the unassigned 

# #This is the example code that garrett Sent me:
#q +geom_point(data=LDdata,aes(x=BP_A,y=BP_B,color=R2),shape=15)+scale_color_gradient(low="black",high="dark red")+theme_bw()
#q +geom_point(data=LDdata[LDdata$R2>=0.3,],aes(x=BP_A,y=BP_B,color=R2),shape=15)+scale_color_gradient(low="black",high="red")+theme_bw()


library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)
library(reshape2)
library(plotly)
library(gridExtra)
library(scales) 
library(grid)
library(RColorBrewer)


#assign each of the Chromosomes a plot that will be called later in the grid arrange 
#this is the backbone of the ggplot that we are building
q <- ggplot() + labs(x="", y="") +scale_color_gradient(low="lightgoldenrodyellow",high="orangered3")+ theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

#The dataframe is too big to plot all the points, so we are going to import the filtered version that we made earlier

LDdata_all_filter <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_all_.2_BPA_BPB_R2.txt",  header=TRUE, 
                                stringsAsFactors = FALSE, na.strings = "-" )

#### ALL 
#facet wrap of all the chromosomes at once. 
cc <- ggplot() + labs(x="", y="") +scale_color_gradient(low="lightgoldenrodyellow",high="orangered3")+ theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+
  geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A != "35"),], aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .5, size= .7) + 
  ggtitle("LD r2 value >.2 ALL Populations") + facet_wrap(~CHR_A, scales= "free")

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_ALL_eachCHR_thumb.jpg", cc,  width = 9,  height = 7,  dpi = 1200)

c0<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "35") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Unassigned Chr") 
c1<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "1") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 1") 
c2<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "2") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 2") 
c3<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "3") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 3") 
c4<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "4") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 4") 
c5<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "5") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 5") 
c6<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "6") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 6") 
c7<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "7") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 7") 
c8<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "8") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 8") 
c9<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "9") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 9") 
c10<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "10") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 10") 
c11<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "11") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 11") 
c12<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "12") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 12") 
c13<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "13") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 13") 
c14<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "14") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 14") 
c15<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "15") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 15") 
c16<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "16") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 16") 
c17<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "17") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 17") 
c18<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "18") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 18") 
c19<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "19") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 19") 
c20<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "20") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 20") 
c21<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "21") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 21") 
c22<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "22") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 22") 
c23<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "23") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 23") 
c24<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "24") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 24") 
c25<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "25") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 25") 
c26<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "26") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 26") 
c27<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "27") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 27") 
c28<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "28") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 28") 
c29<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "29") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 29") 
c30<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "30") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 30") 
c31<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "31") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 31") 
c32<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "32") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 32") 
c33<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "33") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 33") 
c34<-q +geom_point(data=LDdata_all_filter[(LDdata_all_filter$CHR_A == "34") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 34") 

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_ALL_eachCHR.pdf", width = 9, height = 7)

cc #the thumb
c0
c1
c2
c3
c4
c5
c6
c7
c8
c9
c10
c11
c12
c13
c14
c15
c16
c17
c18
c19
c20
c21
c22
c23
c24
c25
c26
c27
c28
c29
c30
c31
c32
c33
c34
dev.off()

rm( LDdata_all_filter, cc, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27, 
    c28, c29, c30, c31, c32, c33, c34 )
gc()


#### Even
LDdata_even_filter <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_even_.2_BPA_BPB_R2.txt",  header=TRUE, 
                                 stringsAsFactors = FALSE, na.strings = "-" )

#facet wrap of all the chromosomes at once. 
ee <- ggplot() + labs(x="", y="") +scale_color_gradient(low="lightgoldenrodyellow",high="orangered3")+ theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+
  geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A != "35"),], aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .5, size= .7) + 
  ggtitle("LD r2 value >.2 Even Populations") + facet_wrap(~CHR_A, scales= "free")

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_Even_eachCHR_thumb.jpg", ee,  width = 9,  height = 7,  dpi = 1200)

#assign each of the Chromosomes a plot that will be called later in the grid arrange 
e0<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "0") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Unassigned Chr")
e1<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "1") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 1")
e2<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "2") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 2")
e3<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "3") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 3")
e4<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "4") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 4")
e5<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "5") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 5")
e6<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "6") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 6")
e7<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "7") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 7")
e8<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "8") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 8")
e9<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "9") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 9")
e10<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "10") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 10")
e11<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "11") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 11")
e12<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "12") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 12")
e13<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "13") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 13")
e14<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "14") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 14")
e15<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "15") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 15")
e16<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "16") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 16")
e17<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "17") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 17")
e18<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "18") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 18")
e19<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "19") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 19")
e20<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "20") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 20")
e21<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "21") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 21")
e22<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "22") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 22")
e23<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "23") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 23")
e24<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "24") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 24")
e25<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "25") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 25")
e26<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "26") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 26")
e27<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "27") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 27")
e28<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "28") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 28")
e29<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "29") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 29")
e30<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "30") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 30")
e31<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "31") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 31")
e32<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "32") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 32")
e33<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "33") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 33")
e34<-q +geom_point(data=LDdata_even_filter[(LDdata_even_filter$CHR_A == "34") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 34")

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_Even_eachCHR.pdf", width = 9, height = 7)

ee #the thumb
e0
e1
e2
e3
e4
e5
e6
e7
e8
e9
e10
e11
e12
e13
e14
e15
e16
e17
e18
e19
e20
e21
e22
e23
e24
e25
e26
e27
e28
e29
e30
e31
e32
e33
e34
dev.off()

rm( LDdata_even_filter, ee, e0, e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13, e14, e15, e16, e17, e18, e19, e20, e21, e22, e23, e24, e25, e26, e27, 
    e28, e29, e30, e31, e32, e33, e34 )
gc()


#### Odd
LDdata_odd_filter <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_odd_.2_BPA_BPB_R2.txt",  header=TRUE, 
                                stringsAsFactors = FALSE, na.strings = "-" )

#facet wrap of all the chromosomes at once. 
oo <- ggplot() + labs(x="", y="") +scale_color_gradient(low="lightgoldenrodyellow",high="orangered3")+ theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+
  geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A != "35"),], aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .5, size= .7) + 
  ggtitle("LD r2 value >.2 Odd Populations") + facet_wrap(~CHR_A, scales= "free")

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_Odd_eachCHR_thumb.jpg", oo,  width = 9,  height = 7,  dpi = 1200)


#assign each of the Chromosomes a plot that will be called later in the grid arrange 
o0<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "35") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Unassigned Chr")
o1<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "1") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 1")
o2<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "2") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 2")
o3<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "3") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 3")
o4<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "4") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 4")
o5<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "5") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 5")
o6<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "6") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 6")
o7<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "7") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 7")
o8<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "8") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 8")
o9<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "9") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 9")
o10<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "10") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 10")
o11<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "11") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 11")
o12<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "12") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 12")
o13<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "13") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 13")
o14<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "14") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 14")
o15<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "15") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 15")
o16<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "16") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 16")
o17<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "17") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 17")
o18<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "18") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 18")
o19<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "19") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 19")
o20<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "20") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 20")
o21<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "21") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 21")
o22<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "22") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 22")
o23<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "23") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 23")
o24<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "24") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 24")
o25<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "25") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 25")
o26<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "26") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 26")
o27<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "27") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 27")
o28<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "28") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 28")
o29<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "29") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 29")
o30<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "30") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 30")
o31<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "31") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 31")
o32<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "32") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 32")
o33<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "33") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 33")
o34<-q +geom_point(data=LDdata_odd_filter[(LDdata_odd_filter$CHR_A == "34") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 34")

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_Odd_eachCHR.pdf", width = 9, height = 7)

oo
o0
o1
o2
o3
o4
o5
o6
o7
o8
o9
o10
o11
o12
o13
o14
o15
o16
o17
o18
o19
o20
o21
o22
o23
o24
o25
o26
o27
o28
o29
o30
o31
o32
o33
o34
dev.off()

rm(LDdata_odd_filter,oo, o0,o1, o2, o3, o4, o5, o6, o7, o8, o9, o10, o11, o12, o13, o14, o15, o16, o17, o18, o19, o20, o21, o22, o23, o24, o25, o26, o27, 
             o28, o29, o30, o31, o32, o33, o34)

#### Even_No Susitna
LDdata_even_ns_filter <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_even_ns_.2_BPA_BPB_R2.txt",  header=TRUE, 
                                stringsAsFactors = FALSE, na.strings = "-" )

#facet wrap of all the chromosomes at once. 
ens <- ggplot() + labs(x="", y="") +scale_color_gradient(low="lightgoldenrodyellow",high="orangered3")+ theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+
  geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A != "35"),], aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .5, size= .7) + 
  ggtitle("LD r2 value >.2 Even No Susitna Populations") + facet_wrap(~CHR_A, scales= "free")

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_even_ns_eachCHR_thumb.jpg", ens,  width = 9,  height = 7,  dpi = 1200)

#assign ensch of the Chromosomes a plot that will be called later in the grid arrange 
ens0<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "35") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Unassigned Chr")
ens1<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "1") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 1")
ens2<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "2") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 2")
ens3<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "3") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 3")
ens4<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "4") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 4")
ens5<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "5") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 5")
ens6<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "6") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 6")
ens7<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "7") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 7")
ens8<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "8") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 8")
ens9<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "9") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 9")
ens10<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "10") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 10")
ens11<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "11") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 11")
ens12<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "12") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 12")
ens13<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "13") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 13")
ens14<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "14") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 14")
ens15<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "15") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 15")
ens16<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "16") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 16")
ens17<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "17") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 17")
ens18<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "18") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 18")
ens19<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "19") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 19")
ens20<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "20") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 20")
ens21<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "21") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 21")
ens22<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "22") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 22")
ens23<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "23") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 23")
ens24<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "24") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 24")
ens25<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "25") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 25")
ens26<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "26") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 26")
ens27<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "27") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 27")
ens28<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "28") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 28")
ens29<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "29") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 29")
ens30<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "30") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 30")
ens31<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "31") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 31")
ens32<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "32") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 32")
ens33<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "33") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 33")
ens34<-q +geom_point(data=LDdata_even_ns_filter[(LDdata_even_ns_filter$CHR_A == "34") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 34")

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_Even_NS_eachCHR.pdf", width = 9, height = 7)

ens
ens0
ens1
ens2
ens3
ens4
ens5
ens6
ens7
ens8
ens9
ens10
ens11
ens12
ens13
ens14
ens15
ens16
ens17
ens18
ens19
ens20
ens21
ens22
ens23
ens24
ens25
ens26
ens27
ens28
ens29
ens30
ens31
ens32
ens33
ens34
dev.off()


rm(LDdata_even_ns_filter, ens, ens0, ens1, ens2, ens3, ens4, ens5, ens6, ens7, ens8, ens9, ens10, ens11, ens12, ens13, ens14, ens15, ens16, ens17, ens18, ens19, ens20, ens21, ens22, ens23, ens24, ens25, ens26, ens27, 
             ens28, ens29, ens30, ens31, ens32, ens33, ens34)
gc()


####  ODD No Susitna

LDdata_odd_ns_filter <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_odd_ns_.2_BPA_BPB_R2.txt",  header=TRUE, 
                                   stringsAsFactors = FALSE, na.strings = "-" )

#facet wrap of all the chromosomes at once. 
ons <- ggplot() + labs(x="", y="") +scale_color_gradient(low="lightgoldenrodyellow",high="orangered3")+ theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+
  geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A != "35"),], aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .5, size= .7) + 
  ggtitle("LD r2 value >.2 Odd No Susitna Populations") + facet_wrap(~CHR_A, scales= "free")

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_Odd_NS_eachCHR_thumb.jpg", ons,  width = 9,  height = 7,  dpi = 1200)


#assign each of the Chromosomes a plot that will be called later in the grid arrange 
ons0<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "35") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Unassigned Chr")
ons1<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "1") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 1")
ons2<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "2") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 2")
ons3<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "3") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 3")
ons4<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "4") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 4")
ons5<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "5") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 5")
ons6<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "6") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 6")
ons7<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "7") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 7")
ons8<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "8") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 8")
ons9<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "9") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 9")
ons10<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "10") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 10")
ons11<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "11") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 11")
ons12<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "12") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 12")
ons13<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "13") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 13")
ons14<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "14") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 14")
ons15<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "15") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 15")
ons16<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "16") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 16")
ons17<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "17") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 17")
ons18<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "18") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 18")
ons19<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "19") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 19")
ons20<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "20") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 20")
ons21<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "21") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 21")
ons22<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "22") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 22")
ons23<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "23") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 23")
ons24<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "24") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 24")
ons25<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "25") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 25")
ons26<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "26") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 26")
ons27<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "27") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 27")
ons28<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "28") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 28")
ons29<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "29") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 29")
ons30<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "30") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 30")
ons31<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "31") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 31")
ons32<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "32") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 32")
ons33<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "33") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 33")
ons34<-q +geom_point(data=LDdata_odd_ns_filter[(LDdata_odd_ns_filter$CHR_A == "34") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 34")

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_Odd_NS_eachCHR.pdf", width = 9, height = 7)
ons
ons0
ons1
ons2
ons3
ons4
ons5
ons6
ons7
ons8
ons9
ons10
ons11
ons12
ons13
ons14
ons15
ons16
ons17
ons18
ons19
ons20
ons21
ons22
ons23
ons24
ons25
ons26
ons27
ons28
ons29
ons30
ons31
ons32
ons33
ons34
dev.off()

rm(LDdata_odd_ns_filter, ons, ons0, ons1, ons2, ons3, ons4, ons5, ons6, ons7, ons8, ons9, ons10, ons11, ons12, ons13, ons14, ons15, ons16, ons17, ons18, ons19, ons20, ons21, ons22, ons23, ons24, ons25, ons26, ons27, 
             ons28, ons29, ons30, ons31, ons32, ons33, ons34)
gc()

#NORTH AMERICA ALL

LD_data_na_all_filter <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_na_all_.2_BPA_BPB_R2.txt",  header=TRUE, 
                                    stringsAsFactors = FALSE, na.strings = "-" )

#facet wrap of all the chromosomes at once. 
naa <- ggplot() + labs(x="", y="") +scale_color_gradient(low="lightgoldenrodyellow",high="orangered3")+ theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+
  geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A != "35"),], aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .5, size= .7) + 
  ggtitle("LD r2 value >.2 ALL North American Populations") + facet_wrap(~CHR_A, scales= "free")

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_NA_ALL_eachCHR_thumb.jpg", naa,  width = 9,  height = 7,  dpi = 1200)


#assign each of the Chromosomes a plot that will be called later in the grid arrange
na0<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "35") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Unassigned Chr")
na1<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "1") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 1")
na2<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "2") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 2")
na3<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "3") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 3")
na4<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "4") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 4")
na5<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "5") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 5")
na6<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "6") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 6")
na7<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "7") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 7")
na8<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "8") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 8")
na9<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "9") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 9")
na10<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "10") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 10")
na11<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "11") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 11")
na12<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "12") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 12")
na13<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "13") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 13")
na14<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "14") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 14")
na15<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "15") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 15")
na16<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "16") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 16")
na17<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "17") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 17")
na18<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "18") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 18")
na19<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "19") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 19")
na20<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "20") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 20")
na21<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "21") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 21")
na22<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "22") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 22")
na23<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "23") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 23")
na24<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "24") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 24")
na25<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "25") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 25")
na26<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "26") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 26")
na27<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "27") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 27")
na28<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "28") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 28")
na29<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "29") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 29")
na30<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "30") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 30")
na31<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "31") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 31")
na32<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "32") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 32")
na33<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "33") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 33")
na34<-q +geom_point(data=LD_data_na_all_filter[(LD_data_na_all_filter$CHR_A == "34") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 34")

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_NA_ALL_eachCHR_thumb.pdf", width = 9, height = 7)

naa
na0
na1
na2
na3
na4
na5
na6
na7
na8
na9
na10
na11
na12
na13
na14
na15
na16
na17
na18
na19
na20
na21
na22
na23
na24
na25
na26
na27
na28
na29
na30
na31
na32
na33
na34
dev.off()

#list the plots and call them in their layout
rm(LD_data_na_all_filter, naa, na0,na1, na2, na3, na4, na5, na6, na7, na8, na9, na10, na11, na12, na13, na14, na15, na16, na17, na18, na19, na20, na21, na22, na23, na24, na25, na26, na27, 
             na28, na29, na30, na31, na32, na33, na34)
gc()


#All Asian Populations 
LD_data_a_all_filter <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_a_all_.2_BPA_BPB_R2.txt",  header=TRUE, 
                                   stringsAsFactors = FALSE, na.strings = "-" )

#facet wrap of all the chromosomes at once. 
aa <- ggplot() + labs(x="", y="") +scale_color_gradient(low="lightgoldenrodyellow",high="orangered3")+ theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+
  geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A != "35"),], aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .5, size= .7) + 
  ggtitle("LD r2 value >.2 ALL Asian Populations") + facet_wrap(~CHR_A, scales= "free")

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_ASIA_ALL_eachCHR_thumb.jpg", aa,  width = 9,  height = 7,  dpi = 1200)

#assign each of the Chromosomes a plot that will be called later in the grid arrange
aa0<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "35") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Unassigned Chr")
aa1<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "1") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 1")
aa2<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "2") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 2")
aa3<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "3") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 3")
aa4<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "4") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 4")
aa5<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "5") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 5")
aa6<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "6") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 6")
aa7<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "7") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 7")
aa8<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "8") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 8")
aa9<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "9") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 9")
aa10<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "10") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 10")
aa11<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "11") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 11")
aa12<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "12") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 12")
aa13<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "13") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 13")
aa14<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "14") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 14")
aa15<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "15") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 15")
aa16<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "16") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 16")
aa17<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "17") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 17")
aa18<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "18") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 18")
aa19<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "19") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 19")
aa20<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "20") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 20")
aa21<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "21") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 21")
aa22<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "22") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 22")
aa23<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "23") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 23")
aa24<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "24") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 24")
aa25<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "25") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 25")
aa26<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "26") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 26")
aa27<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "27") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 27")
aa28<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "28") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 28")
aa29<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "29") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 29")
aa30<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "30") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 30")
aa31<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "31") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 31")
aa32<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "32") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 32")
aa33<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "33") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 33")
aa34<-q +geom_point(data=LD_data_a_all_filter[(LD_data_a_all_filter$CHR_A == "34") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 34")

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_ASIA_ALL_eachCHR_thumb.pdf", width = 9, height = 7)
aa
aa0
aa1
aa2
aa3
aa4
aa5
aa6
aa7
aa8
aa9
aa10
aa11
aa12
aa13
aa14
aa15
aa16
aa17
aa18
aa19
aa20
aa21
aa22
aa23
aa24
aa25
aa26
aa27
aa28
aa29
aa30
aa31
aa32
aa33
aa34
dev.off()

#list the plots and call them in their layout
rm(LD_data_a_all_filter, aa, aa0, aa1, aa2, aa3, aa4, aa5, aa6, aa7, aa8, aa9, aa10, aa11, aa12, aa13, aa14, aa15, aa16, aa17, aa18, aa19, aa20, aa21, aa22, aa23, aa24, aa25, aa26, aa27, 
             aa28, aa29, aa30, aa31, aa32, aa33, aa34)

