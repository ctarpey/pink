### Plotting the LD r2 from PLINK
###     REGIONAL GROUPINGS 
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
LDdata_even_A_filter <-  read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_even_A_.2_BPA_BPB_R2.txt",  header=TRUE, 
                             stringsAsFactors = FALSE, na.strings = "-" )

#### ALL 
#facet wrap of all the chromosomes at once. 
ea <- ggplot() + labs(x="", y="") +scale_color_gradient(low="lightgoldenrodyellow",high="orangered3")+ theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+
  geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A != "35"),], aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .5, size= .7) + 
  ggtitle("LD r2 value >.2 Even Asian Populations") + facet_wrap(~CHR_A, scales= "free")

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_even_A_eachCHR_thumb.jpg", ea,  width = 9,  height = 7,  dpi = 1200)

ea0<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "35") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Unassigned Chr") 
ea1<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "1") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 1") 
ea2<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "2") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 2") 
ea3<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "3") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 3") 
ea4<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "4") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 4") 
ea5<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "5") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 5") 
ea6<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "6") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 6") 
ea7<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "7") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 7") 
ea8<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "8") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 8") 
ea9<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "9") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 9") 
ea10<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "10") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 10") 
ea11<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "11") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 11") 
ea12<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "12") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 12") 
ea13<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "13") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 13") 
ea14<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "14") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 14") 
ea15<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "15") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 15") 
ea16<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "16") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 16") 
ea17<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "17") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 17") 
ea18<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "18") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 18") 
ea19<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "19") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 19") 
ea20<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "20") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 20") 
ea21<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "21") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 21") 
ea22<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "22") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 22") 
ea23<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "23") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 23") 
ea24<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "24") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 24") 
ea25<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "25") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 25") 
ea26<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "26") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 26") 
ea27<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "27") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 27") 
ea28<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "28") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 28") 
ea29<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "29") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 29") 
ea30<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "30") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 30") 
ea31<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "31") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 31") 
ea32<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "32") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 32") 
ea33<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "33") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 33") 
ea34<-q +geom_point(data=LDdata_even_A_filter[(LDdata_even_A_filter$CHR_A == "34") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Chr. 34") 

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LDdata_even_A_eachCHR.pdf", width = 9, height = 7)
ea #the thumb
ea0
ea1
ea2
ea3
ea4
ea5
ea6
ea7
ea8
ea9
ea10
ea11
ea12
ea13
ea14
ea15
ea16
ea17
ea18
ea19
ea20
ea21
ea22
ea23
ea24
ea25
ea26
ea27
ea28
ea29
ea30
ea31
ea32
ea33
ea34
dev.off()

rm( LDdata_even_A_filter, ea, ea0, ea1, ea2, ea3, ea4, ea5, ea6, ea7, ea8, ea9, ea10, ea11, ea12, ea13, ea14, ea15, ea16, ea17, ea18, ea19, ea20, ea21, ea22, ea23, ea24, ea25, ea26, ea27, 
    ea28, ea29, ea30, ea31, ea32, ea33, ea34 )
gc()


#### Even North America
LDdata_even_NA_filter <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_even_NA_.2_BPA_BPB_R2.txt",  header=TRUE, 
                                    stringsAsFactors = FALSE, na.strings = "-" )


#facet wrap of all the chromosomes at once. 
ena <- ggplot() + labs(x="", y="") +scale_color_gradient(low="lightgoldenrodyellow",high="orangered3")+ theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+
  geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A != "35"),], aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .5, size= .7) + 
  ggtitle("LD r2 value >.2 Even North American Populations") + facet_wrap(~CHR_A, scales= "free")

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_Even_NA_eachCHR_thumb.jpg", ena,  width = 9,  height = 7,  dpi = 1200)

#assign each of the Chromosomes a plot that will be called later in the grid arrange 
ena0<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "35") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65) + ggtitle("Unassigned Chr") 
ena1<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "1") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 1")
ena2<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "2") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 2")
ena3<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "3") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 3")
ena4<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "4") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 4")
ena5<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "5") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 5")
ena6<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "6") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 6")
ena7<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "7") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 7")
ena8<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "8") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 8")
ena9<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "9") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 9")
ena10<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "10") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 10")
ena11<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "11") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 11")
ena12<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "12") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 12")
ena13<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "13") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 13")
ena14<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "14") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 14")
ena15<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "15") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 15")
ena16<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "16") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 16")
ena17<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "17") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 17")
ena18<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "18") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 18")
ena19<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "19") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 19")
ena20<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "20") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 20")
ena21<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "21") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 21")
ena22<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "22") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 22")
ena23<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "23") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 23")
ena24<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "24") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 24")
ena25<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "25") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 25")
ena26<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "26") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 26")
ena27<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "27") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 27")
ena28<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "28") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 28")
ena29<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "29") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 29")
ena30<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "30") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 30")
ena31<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "31") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 31")
ena32<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "32") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 32")
ena33<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "33") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 33")
ena34<-q +geom_point(data=LDdata_even_NA_filter[(LDdata_even_NA_filter$CHR_A == "34") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 34")

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_Even_NA_eachCHR.pdf", width = 9, height = 7)

ena #the thumb
ena0
ena1
ena2
ena3
ena4
ena5
ena6
ena7
ena8
ena9
ena10
ena11
ena12
ena13
ena14
ena15
ena16
ena17
ena18
ena19
ena20
ena21
ena22
ena23
ena24
ena25
ena26
ena27
ena28
ena29
ena30
ena31
ena32
ena33
ena34
dev.off()

rm( LDdata_even_NA_filter, ena, ena0, ena1, ena2, ena3, ena4, ena5, ena6, ena7, ena8, ena9, ena10, ena11, ena12, ena13, ena14, ena15, ena16, ena17, ena18, ena19, ena20, ena21, ena22, ena23, ena24, ena25, ena26, ena27, 
    ena28, ena29, ena30, ena31, ena32, ena33, ena34 )
gc()


#### Odd Asian 
LDdata_odd_a_filter<-  read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_odd_a_.2_BPA_BPB_R2.txt",  header=TRUE, 
                                  stringsAsFactors = FALSE, na.strings = "-" )

#facet wrap of all the chromosomes at once. 
oa <- ggplot() + labs(x="", y="") +scale_color_gradient(low="lightgoldenrodyellow",high="orangered3")+ theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+
  geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A != "35"),], aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .5, size= .7) + 
  ggtitle("LD r2 value >.2 Odd Asian Populations") + facet_wrap(~CHR_A, scales= "free")

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_Odd_A_eachCHR_thumb.jpg", oa,  width = 9,  height = 7,  dpi = 1200)


#assign each of the Chromosomes a plot that will be called later in the grid arrange 
oa0<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "35") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Unassigned Chr")
oa1<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "1") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 1")
oa2<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "2") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 2")
oa3<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "3") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 3")
oa4<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "4") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 4")
oa5<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "5") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 5")
oa6<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "6") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 6")
oa7<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "7") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 7")
oa8<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "8") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 8")
oa9<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "9") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 9")
oa10<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "10") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 10")
oa11<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "11") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 11")
oa12<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "12") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 12")
oa13<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "13") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 13")
oa14<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "14") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 14")
oa15<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "15") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 15")
oa16<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "16") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 16")
oa17<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "17") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 17")
oa18<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "18") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 18")
oa19<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "19") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 19")
oa20<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "20") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 20")
oa21<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "21") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 21")
oa22<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "22") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 22")
oa23<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "23") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 23")
oa24<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "24") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 24")
oa25<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "25") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 25")
oa26<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "26") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 26")
oa27<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "27") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 27")
oa28<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "28") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 28")
oa29<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "29") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 29")
oa30<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "30") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 30")
oa31<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "31") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 31")
oa32<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "32") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 32")
oa33<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "33") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 33")
oa34<-q +geom_point(data=LDdata_odd_a_filter[(LDdata_odd_a_filter$CHR_A == "34") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 34")

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_Odd_A_eachCHR.pdf", width = 9, height = 7)
oa
oa0
oa1
oa2
oa3
oa4
oa5
oa6
oa7
oa8
oa9
oa10
oa11
oa12
oa13
oa14
oa15
oa16
oa17
oa18
oa19
oa20
oa21
oa22
oa23
oa24
oa25
oa26
oa27
oa28
oa29
oa30
oa31
oa32
oa33
oa34
dev.off()

rm(LDdata_odd_a_filter,oa, oa0, oa1, oa2, oa3, oa4, oa5, oa6, oa7, oa8, oa9, oa10, oa11, oa12, oa13, oa14, oa15, oa16, oa17, oa18, oa19, oa20, oa21, oa22, oa23, oa24, oa25, oa26, oa27, 
   oa28, oa29, oa30, oa31, oa32, oa33, oa34)

#### Odd North America 
LDdata_odd_na_filter <-  read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_odd_na_.2_BPA_BPB_R2.txt",  header=TRUE, 
                                    stringsAsFactors = FALSE, na.strings = "-" )

#facet wrap of all the chromosomes at once. 
ona <- ggplot() + labs(x="", y="") +scale_color_gradient(low="lightgoldenrodyellow",high="orangered3")+ theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+
  geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A != "35"),], aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .5, size= .7) + 
  ggtitle("LD r2 value >.2 Odd North American Populations") + facet_wrap(~CHR_A, scales= "free")

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_Odd_NA_eachCHR_thumb.jpg", ona,  width = 9,  height = 7,  dpi = 1200)

#assign onach of the Chromosomes a plot that will be called later in the grid arrange 
ona0<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "35") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Unassigned Chr")
ona1<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "1") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 1")
ona2<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "2") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 2")
ona3<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "3") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 3")
ona4<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "4") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 4")
ona5<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "5") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 5")
ona6<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "6") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 6")
ona7<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "7") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 7")
ona8<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "8") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 8")
ona9<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "9") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 9")
ona10<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "10") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 10")
ona11<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "11") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 11")
ona12<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "12") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 12")
ona13<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "13") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 13")
ona14<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "14") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 14")
ona15<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "15") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 15")
ona16<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "16") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 16")
ona17<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "17") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 17")
ona18<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "18") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 18")
ona19<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "19") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 19")
ona20<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "20") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 20")
ona21<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "21") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 21")
ona22<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "22") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 22")
ona23<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "23") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 23")
ona24<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "24") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 24")
ona25<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "25") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 25")
ona26<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "26") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 26")
ona27<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "27") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 27")
ona28<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "28") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 28")
ona29<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "29") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 29")
ona30<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "30") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 30")
ona31<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "31") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 31")
ona32<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "32") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 32")
ona33<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "33") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 33")
ona34<-q +geom_point(data=LDdata_odd_na_filter[(LDdata_odd_na_filter$CHR_A == "34") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 34")

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_Odd_NA_eachCHR.pdf", width = 9, height = 7)

ona
ona0
ona1
ona2
ona3
ona4
ona5
ona6
ona7
ona8
ona9
ona10
ona11
ona12
ona13
ona14
ona15
ona16
ona17
ona18
ona19
ona20
ona21
ona22
ona23
ona24
ona25
ona26
ona27
ona28
ona29
ona30
ona31
ona32
ona33
ona34
dev.off()


rm(LDdata_odd_na_filter, ona, ona0, ona1, ona2, ona3, ona4, ona5, ona6, ona7, ona8, ona9, ona10, ona11, ona12, ona13, ona14, ona15, ona16, ona17, ona18, ona19, ona20, ona21, ona22, ona23, ona24, ona25, ona26, ona27, 
   ona28, ona29, ona30, ona31, ona32, ona33, ona34)
gc()


####  Even North America  No Susitna

LDdata_even_na_ns_filter <-  read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_even_na_ns_.2_BPA_BPB_R2.txt",  header=TRUE, 
                                        stringsAsFactors = FALSE, na.strings = "-" )

#facet wrap of all the chromosomes at once. 
enans <- ggplot() + labs(x="", y="") +scale_color_gradient(low="lightgoldenrodyellow",high="orangered3")+ theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+
  geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A != "35"),], aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .5, size= .7) + 
  ggtitle("LD r2 value >.2 Even NA No Susitna Populations") + facet_wrap(~CHR_A, scales= "free")

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_Even_NA_NS_eachCHR_thumb.jpg", enans,  width = 9,  height = 7,  dpi = 1200)


#assign each of the Chromosomes a plot that will be called later in the grid arrange 
enans0<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "35") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Unassigned Chr")
enans1<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "1") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 1")
enans2<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "2") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 2")
enans3<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "3") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 3")
enans4<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "4") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 4")
enans5<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "5") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 5")
enans6<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "6") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 6")
enans7<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "7") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 7")
enans8<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "8") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 8")
enans9<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "9") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 9")
enans10<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "10") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 10")
enans11<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "11") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 11")
enans12<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "12") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 12")
enans13<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "13") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 13")
enans14<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "14") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 14")
enans15<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "15") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 15")
enans16<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "16") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 16")
enans17<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "17") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 17")
enans18<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "18") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 18")
enans19<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "19") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 19")
enans20<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "20") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 20")
enans21<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "21") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 21")
enans22<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "22") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 22")
enans23<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "23") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 23")
enans24<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "24") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 24")
enans25<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "25") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 25")
enans26<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "26") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 26")
enans27<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "27") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 27")
enans28<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "28") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 28")
enans29<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "29") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 29")
enans30<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "30") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 30")
enans31<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "31") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 31")
enans32<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "32") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 32")
enans33<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "33") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 33")
enans34<-q +geom_point(data=LDdata_even_na_ns_filter[(LDdata_even_na_ns_filter$CHR_A == "34") ,],aes (x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 34")

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_EVEN_NA_NS_eachCHR.pdf", width = 9, height = 7)
enans
enans0
enans1
enans2
enans3
enans4
enans5
enans6
enans7
enans8
enans9
enans10
enans11
enans12
enans13
enans14
enans15
enans16
enans17
enans18
enans19
enans20
enans21
enans22
enans23
enans24
enans25
enans26
enans27
enans28
enans29
enans30
enans31
enans32
enans33
enans34
dev.off()

rm(LDdata_even_na_ns_filter, enans, enans0, enans1, enans2, enans3, enans4, enans5, enans6, enans7, enans8, enans9, enans10, enans11, enans12, enans13, enans14, enans15, enans16, enans17, enans18, enans19, enans20, enans21, enans22, enans23, enans24, enans25, enans26, enans27, 
   enans28, enans29, enans30, enans31, enans32, enans33, enans34)
gc()

## Odd North America No Susitna
LDdata_odd_na_ns_filter <-  read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_data_odd_na_ns_.2_BPA_BPB_R2.txt",  header=TRUE, 
                                       stringsAsFactors = FALSE, na.strings = "-" )

#facet wrap of all the chromosomes at once. 
onans <- ggplot() + labs(x="", y="") +scale_color_gradient(low="lightgoldenrodyellow",high="orangered3")+ theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+
  geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A != "35"),], aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .5, size= .7) + 
  ggtitle("LD r2 value >.2 Odd NA No Susitna Populations") + facet_wrap(~CHR_A, scales= "free")

ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_ONA_NS_eachCHR_thumb.jpg", onans,  width = 9,  height = 7,  dpi = 1200)


#assign each of the Chromosomes a plot that will be called later in the grid arrange
onans0<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "35") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Unasssigned Chr")
onans1<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "1") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 1")
onans2<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "2") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 2")
onans3<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "3") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 3")
onans4<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "4") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 4")
onans5<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "5") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 5")
onans6<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "6") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 6")
onans7<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "7") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 7")
onans8<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "8") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 8")
onans9<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "9") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 9")
onans10<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "10") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 10")
onans11<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "11") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 11")
onans12<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "12") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 12")
onans13<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "13") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 13")
onans14<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "14") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 14")
onans15<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "15") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 15")
onans16<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "16") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 16")
onans17<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "17") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 17")
onans18<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "18") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 18")
onans19<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "19") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 19")
onans20<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "20") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 20")
onans21<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "21") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 21")
onans22<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "22") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 22")
onans23<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "23") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 23")
onans24<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "24") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 24")
onans25<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "25") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 25")
onans26<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "26") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 26")
onans27<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "27") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 27")
onans28<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "28") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 28")
onans29<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "29") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 29")
onans30<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "30") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 30")
onans31<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "31") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 31")
onans32<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "32") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 32")
onans33<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "33") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 33")
onans34<-q +geom_point(data=LDdata_odd_na_ns_filter[(LDdata_odd_na_ns_filter$CHR_A == "34") ,],aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .65)+ ggtitle("Chr. 34")

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/LD_r2_ONA_NS_eachCHR_thumb.pdf", width = 9, height = 7)

onans
onans0
onans1
onans2
onans3
onans4
onans5
onans6
onans7
onans8
onans9
onans10
onans11
onans12
onans13
onans14
onans15
onans16
onans17
onans18
onans19
onans20
onans21
onans22
onans23
onans24
onans25
onans26
onans27
onans28
onans29
onans30
onans31
onans32
onans33
onans34

dev.off()

#list the plots and call them in their layout
rm(LDdata_odd_na_ns_filter, onans, onans0,onans1, onans2, onans3, onans4, onans5, onans6, onans7, onans8, onans9, onans10, onans11, onans12, onans13, onans14, onans15, onans16, onans17, onans18, onans19, onans20, onans21, onans22, onans23, onans24, onans25, onans26, onans27, 
   onans28, onans29, onans30, onans31, onans32, onans33, onans34)
gc()

