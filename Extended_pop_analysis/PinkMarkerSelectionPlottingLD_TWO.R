### Plotting the LD r2 from PLINK
###    for the pink panel of 2312 markers
###    THIS IS THE SECOND SHEET FOR THE ACTUAL PLOTS OF THE CHROMOSOMES
### Garrett McKinney and Carolyn Tarpey | May 2018
## ---------------------------------------

## THIS IS THE SECOND SHEET- YOU NEED THE MAIN SHEET TO IMPORT THE DATA FOR EACH OF THESE PLOTS

#This R code was originally written by Garrett McKinney to plot the output of PLINK LD r2 for each chromosome. 
#It requires the LD output from PLINK. We used the alignment of the pinks to Chinook to get the Chromosome 
#and position assignments here, so there are 34 Chromosomes, and chromosome 35 are those that aligned to the unassigned 

# #This is the example code that garrett Sent me:
# ggplot()+geom_point(data=LDdata,aes(x=BP_A,y=BP_B,color=R2),shape=15)+scale_color_gradient(low="black",high="dark red")+theme_bw()
# ggplot()+geom_point(data=LDdata[LDdata$R2>=0.3,],aes(x=BP_A,y=BP_B,color=R2),shape=15)+scale_color_gradient(low="black",high="red")+theme_bw()


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

#load PLINK LD output files
#### 2312 ALL 


#assign each of the Chromosomes a plot that will be called later in the grid arrange 
c0<- ggplot()+geom_point(data=LDdata_chr0,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Unassigned Chr") 
c1<- ggplot()+geom_point(data=LDdata_chr1,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 1") 
c2<- ggplot()+geom_point(data=LDdata_chr2,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 2") 
c3<- ggplot()+geom_point(data=LDdata_chr3,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 3") 
c4<- ggplot()+geom_point(data=LDdata_chr4,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 4") 
c5<- ggplot()+geom_point(data=LDdata_chr5,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 5") 
c6<- ggplot()+geom_point(data=LDdata_chr6,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 6") 
c7<- ggplot()+geom_point(data=LDdata_chr7,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 7") 
c8<- ggplot()+geom_point(data=LDdata_chr8,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 8") 
c9<- ggplot()+geom_point(data=LDdata_chr9,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 9") 
c10<- ggplot()+geom_point(data=LDdata_chr10,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 10") 
c11<- ggplot()+geom_point(data=LDdata_chr11,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 11") 
c12<- ggplot()+geom_point(data=LDdata_chr12,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 12") 
c13<- ggplot()+geom_point(data=LDdata_chr13,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 13") 
c14<- ggplot()+geom_point(data=LDdata_chr14,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 14") 
c15<- ggplot()+geom_point(data=LDdata_chr15,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 15") 
c16<- ggplot()+geom_point(data=LDdata_chr16,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 16") 
c17<- ggplot()+geom_point(data=LDdata_chr17,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 17") 
c18<- ggplot()+geom_point(data=LDdata_chr18,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 18") 
c19<- ggplot()+geom_point(data=LDdata_chr19,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 19") 
c20<- ggplot()+geom_point(data=LDdata_chr20,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 20") 
c21<- ggplot()+geom_point(data=LDdata_chr21,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 21") 
c22<- ggplot()+geom_point(data=LDdata_chr22,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 22") 
c23<- ggplot()+geom_point(data=LDdata_chr23,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 23") 
c24<- ggplot()+geom_point(data=LDdata_chr24,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 24") 
c25<- ggplot()+geom_point(data=LDdata_chr25,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 25") 
c26<- ggplot()+geom_point(data=LDdata_chr26,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 26") 
c27<- ggplot()+geom_point(data=LDdata_chr27,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 27") 
c28<- ggplot()+geom_point(data=LDdata_chr28,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 28") 
c29<- ggplot()+geom_point(data=LDdata_chr29,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 29") 
c30<- ggplot()+geom_point(data=LDdata_chr30,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 30") 
c31<- ggplot()+geom_point(data=LDdata_chr31,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 31") 
c32<- ggplot()+geom_point(data=LDdata_chr32,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 32") 
c33<- ggplot()+geom_point(data=LDdata_chr33,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 33") 
c34<- ggplot()+geom_point(data=LDdata_chr34,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("All individuals, Chr. 34") 

#list the plots and call them in their layout
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/LD_r2_ALL_pink_panel_2312.pdf", width = 9, height = 7)

grid.arrange(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27, 
             c28, c29, c30, c31, c32, c33, c34, top="LD r2 value between ALL pink panel markers; positions based on Chinook chromosomes")

#plot the LD of the markers that did not align to any of the Chinook chromosomes sepearate: 
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

#### 2312 Even

#assign each of the Chromosomes a plot that will be called later in the grid arrange 
e0<- ggplot()+geom_point(data=LD_even_chr0 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Unassigned Chr")
e1<- ggplot()+geom_point(data=LD_even_chr1 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 1")
e2<- ggplot()+geom_point(data=LD_even_chr2 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 2")
e3<- ggplot()+geom_point(data=LD_even_chr3 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 3")
e4<- ggplot()+geom_point(data=LD_even_chr4 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 4")
e5<- ggplot()+geom_point(data=LD_even_chr5 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 5")
e6<- ggplot()+geom_point(data=LD_even_chr6 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 6")
e7<- ggplot()+geom_point(data=LD_even_chr7 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 7")
e8<- ggplot()+geom_point(data=LD_even_chr8 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 8")
e9<- ggplot()+geom_point(data=LD_even_chr9 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 9")
e10<- ggplot()+geom_point(data=LD_even_chr10,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 10")
e11<- ggplot()+geom_point(data=LD_even_chr11,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 11")
e12<- ggplot()+geom_point(data=LD_even_chr12,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 12")
e13<- ggplot()+geom_point(data=LD_even_chr13,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 13")
e14<- ggplot()+geom_point(data=LD_even_chr14,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 14")
e15<- ggplot()+geom_point(data=LD_even_chr15,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 15")
e16<- ggplot()+geom_point(data=LD_even_chr16,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 16")
e17<- ggplot()+geom_point(data=LD_even_chr17,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 17")
e18<- ggplot()+geom_point(data=LD_even_chr18,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 18")
e19<- ggplot()+geom_point(data=LD_even_chr19,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 19")
e20<- ggplot()+geom_point(data=LD_even_chr20,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 20")
e21<- ggplot()+geom_point(data=LD_even_chr21,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 21")
e22<- ggplot()+geom_point(data=LD_even_chr22,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 22")
e23<- ggplot()+geom_point(data=LD_even_chr23,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 23")
e24<- ggplot()+geom_point(data=LD_even_chr24,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 24")
e25<- ggplot()+geom_point(data=LD_even_chr25,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 25")
e26<- ggplot()+geom_point(data=LD_even_chr26,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 26")
e27<- ggplot()+geom_point(data=LD_even_chr27,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 27")
e28<- ggplot()+geom_point(data=LD_even_chr28,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 28")
e29<- ggplot()+geom_point(data=LD_even_chr29,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 29")
e30<- ggplot()+geom_point(data=LD_even_chr30,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 30")
e31<- ggplot()+geom_point(data=LD_even_chr31,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 31")
e32<- ggplot()+geom_point(data=LD_even_chr32,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 32")
e33<- ggplot()+geom_point(data=LD_even_chr33,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 33")
e34<- ggplot()+geom_point(data=LD_even_chr34,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even individuals, Chr. 34")

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/LD_r2_aLL_pink_panel_EVEN_2312.pdf", width = 9, height = 7)

grid.arrange(e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12, e13, e14, e15, e16, e17, e18, e19, e20, e21, e22, e23, e24, e25, e26, e27, 
             e28, e29, e30, e31, e32, e33, e34,  top="LD r2 value between all pink panel markers in the even lineage; positions based on Chinook chromosomes")

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

#### 2312 Odd

#assign each of the Chromosomes a plot that will be called later in the grid arrange 
o0<- ggplot()+geom_point(data=LD_odd_chr0 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Unassigned Chr")
o1<- ggplot()+geom_point(data=LD_odd_chr1 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 1")
o2<- ggplot()+geom_point(data=LD_odd_chr2 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 2")
o3<- ggplot()+geom_point(data=LD_odd_chr3 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 3")
o4<- ggplot()+geom_point(data=LD_odd_chr4 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 4")
o5<- ggplot()+geom_point(data=LD_odd_chr5 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 5")
o6<- ggplot()+geom_point(data=LD_odd_chr6 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 6")
o7<- ggplot()+geom_point(data=LD_odd_chr7 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 7")
o8<- ggplot()+geom_point(data=LD_odd_chr8 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 8")
o9<- ggplot()+geom_point(data=LD_odd_chr9 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 9")
o10<- ggplot()+geom_point(data=LD_odd_chr10,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 10")
o11<- ggplot()+geom_point(data=LD_odd_chr11,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 11")
o12<- ggplot()+geom_point(data=LD_odd_chr12,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 12")
o13<- ggplot()+geom_point(data=LD_odd_chr13,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 13")
o14<- ggplot()+geom_point(data=LD_odd_chr14,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 14")
o15<- ggplot()+geom_point(data=LD_odd_chr15,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 15")
o16<- ggplot()+geom_point(data=LD_odd_chr16,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 16")
o17<- ggplot()+geom_point(data=LD_odd_chr17,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 17")
o18<- ggplot()+geom_point(data=LD_odd_chr18,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 18")
o19<- ggplot()+geom_point(data=LD_odd_chr19,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 19")
o20<- ggplot()+geom_point(data=LD_odd_chr20,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 20")
o21<- ggplot()+geom_point(data=LD_odd_chr21,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 21")
o22<- ggplot()+geom_point(data=LD_odd_chr22,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 22")
o23<- ggplot()+geom_point(data=LD_odd_chr23,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 23")
o24<- ggplot()+geom_point(data=LD_odd_chr24,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 24")
o25<- ggplot()+geom_point(data=LD_odd_chr25,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 25")
o26<- ggplot()+geom_point(data=LD_odd_chr26,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 26")
o27<- ggplot()+geom_point(data=LD_odd_chr27,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 27")
o28<- ggplot()+geom_point(data=LD_odd_chr28,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 28")
o29<- ggplot()+geom_point(data=LD_odd_chr29,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 29")
o30<- ggplot()+geom_point(data=LD_odd_chr30,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 30")
o31<- ggplot()+geom_point(data=LD_odd_chr31,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 31")
o32<- ggplot()+geom_point(data=LD_odd_chr32,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 32")
o33<- ggplot()+geom_point(data=LD_odd_chr33,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 33")
o34<- ggplot()+geom_point(data=LD_odd_chr34,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd individuals, Chr. 34")

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/LD_r2_aLL_pink_panel_Odd_2312.pdf", width = 9, height = 7)
#list the plots and call them in their layout
grid.arrange(o1, o2, o3, o4, o5, o6, o7, o8, o9, o10, o11, o12, o13, o14, o15, o16, o17, o18, o19, o20, o21, o22, o23, o24, o25, o26, o27, 
             o28, o29, o30, o31, o32, o33, o34, top="LD r2 value between all pink panel markers in the odd lineage; positions based on Chinook chromosomes")

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

#### 2312 Even Asia

#assign each of the Chromosomes a plot that will be called later in the grid arrange 
ea0<- ggplot()+geom_point(data=LD_even_A_chr0 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Unassigned Chr")
ea1<- ggplot()+geom_point(data=LD_even_A_chr1 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 1")
ea2<- ggplot()+geom_point(data=LD_even_A_chr2 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 2")
ea3<- ggplot()+geom_point(data=LD_even_A_chr3 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 3")
ea4<- ggplot()+geom_point(data=LD_even_A_chr4 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 4")
ea5<- ggplot()+geom_point(data=LD_even_A_chr5 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 5")
ea6<- ggplot()+geom_point(data=LD_even_A_chr6 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 6")
ea7<- ggplot()+geom_point(data=LD_even_A_chr7 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 7")
ea8<- ggplot()+geom_point(data=LD_even_A_chr8 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 8")
ea9<- ggplot()+geom_point(data=LD_even_A_chr9 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 9")
ea10<- ggplot()+geom_point(data=LD_even_A_chr10,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 10")
ea11<- ggplot()+geom_point(data=LD_even_A_chr11,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 11")
ea12<- ggplot()+geom_point(data=LD_even_A_chr12,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 12")
ea13<- ggplot()+geom_point(data=LD_even_A_chr13,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 13")
ea14<- ggplot()+geom_point(data=LD_even_A_chr14,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 14")
ea15<- ggplot()+geom_point(data=LD_even_A_chr15,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 15")
ea16<- ggplot()+geom_point(data=LD_even_A_chr16,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 16")
ea17<- ggplot()+geom_point(data=LD_even_A_chr17,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 17")
ea18<- ggplot()+geom_point(data=LD_even_A_chr18,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 18")
ea19<- ggplot()+geom_point(data=LD_even_A_chr19,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 19")
ea20<- ggplot()+geom_point(data=LD_even_A_chr20,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 20")
ea21<- ggplot()+geom_point(data=LD_even_A_chr21,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 21")
ea22<- ggplot()+geom_point(data=LD_even_A_chr22,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 22")
ea23<- ggplot()+geom_point(data=LD_even_A_chr23,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 23")
ea24<- ggplot()+geom_point(data=LD_even_A_chr24,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 24")
ea25<- ggplot()+geom_point(data=LD_even_A_chr25,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 25")
ea26<- ggplot()+geom_point(data=LD_even_A_chr26,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 26")
ea27<- ggplot()+geom_point(data=LD_even_A_chr27,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 27")
ea28<- ggplot()+geom_point(data=LD_even_A_chr28,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 28")
ea29<- ggplot()+geom_point(data=LD_even_A_chr29,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 29")
ea30<- ggplot()+geom_point(data=LD_even_A_chr30,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 30")
ea31<- ggplot()+geom_point(data=LD_even_A_chr31,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 31")
ea32<- ggplot()+geom_point(data=LD_even_A_chr32,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 32")
ea33<- ggplot()+geom_point(data=LD_even_A_chr33,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 33")
ea34<- ggplot()+geom_point(data=LD_even_A_chr34,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even Asia individuals, Chr. 34")

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/LD_r2_aLL_pink_panel_EVEN_Asia_2312.pdf", width = 9, height = 7)
#list the plots and call them in their layout
grid.arrange(ea1, ea2, ea3, ea4, ea5, ea6, ea7, ea8, ea9, ea10, ea11, ea12, ea13, ea14, ea15, ea16, ea17, ea18, ea19, ea20, ea21, ea22, ea23, ea24, ea25, ea26, ea27, 
             ea28, ea29, ea30, ea31, ea32, ea33, ea34, top="LD r2 value between all pink panel markers in the even lineage, Asian populations; positions based on Chinook chromosomes")

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

#### 2312 Even North America

#assign each of the Chromosomes a plot that will be called later in the grid arrange 
ena0<- ggplot()+geom_point(data=LD_even_NA_chr0 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Unassigned Chr")
ena1<- ggplot()+geom_point(data=LD_even_NA_chr1 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 1")
ena2<- ggplot()+geom_point(data=LD_even_NA_chr2 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 2")
ena3<- ggplot()+geom_point(data=LD_even_NA_chr3 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 3")
ena4<- ggplot()+geom_point(data=LD_even_NA_chr4 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 4")
ena5<- ggplot()+geom_point(data=LD_even_NA_chr5 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 5")
ena6<- ggplot()+geom_point(data=LD_even_NA_chr6 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 6")
ena7<- ggplot()+geom_point(data=LD_even_NA_chr7 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 7")
ena8<- ggplot()+geom_point(data=LD_even_NA_chr8 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 8")
ena9<- ggplot()+geom_point(data=LD_even_NA_chr9 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 9")
ena10<- ggplot()+geom_point(data=LD_even_NA_chr10,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 10")
ena11<- ggplot()+geom_point(data=LD_even_NA_chr11,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 11")
ena12<- ggplot()+geom_point(data=LD_even_NA_chr12,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 12")
ena13<- ggplot()+geom_point(data=LD_even_NA_chr13,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 13")
ena14<- ggplot()+geom_point(data=LD_even_NA_chr14,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 14")
ena15<- ggplot()+geom_point(data=LD_even_NA_chr15,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 15")
ena16<- ggplot()+geom_point(data=LD_even_NA_chr16,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 16")
ena17<- ggplot()+geom_point(data=LD_even_NA_chr17,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 17")
ena18<- ggplot()+geom_point(data=LD_even_NA_chr18,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 18")
ena19<- ggplot()+geom_point(data=LD_even_NA_chr19,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 19")
ena20<- ggplot()+geom_point(data=LD_even_NA_chr20,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 20")
ena21<- ggplot()+geom_point(data=LD_even_NA_chr21,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 21")
ena22<- ggplot()+geom_point(data=LD_even_NA_chr22,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 22")
ena23<- ggplot()+geom_point(data=LD_even_NA_chr23,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 23")
ena24<- ggplot()+geom_point(data=LD_even_NA_chr24,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 24")
ena25<- ggplot()+geom_point(data=LD_even_NA_chr25,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 25")
ena26<- ggplot()+geom_point(data=LD_even_NA_chr26,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 26")
ena27<- ggplot()+geom_point(data=LD_even_NA_chr27,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 27")
ena28<- ggplot()+geom_point(data=LD_even_NA_chr28,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 28")
ena29<- ggplot()+geom_point(data=LD_even_NA_chr29,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 29")
ena30<- ggplot()+geom_point(data=LD_even_NA_chr30,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 30")
ena31<- ggplot()+geom_point(data=LD_even_NA_chr31,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 31")
ena32<- ggplot()+geom_point(data=LD_even_NA_chr32,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 32")
ena33<- ggplot()+geom_point(data=LD_even_NA_chr33,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 33")
ena34<- ggplot()+geom_point(data=LD_even_NA_chr34,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Even NA individuals, Chr. 34")

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/LD_r2_aLL_pink_panel_EVEN_NA_2312.pdf", width = 9, height = 7)
#list the plots and call them in their layout
grid.arrange(ena1, ena2, ena3, ena4, ena5, ena6, ena7, ena8, ena9, ena10, ena11, ena12, ena13, ena14, ena15, ena16, ena17, ena18, ena19, ena20, ena21, ena22, ena23, ena24, ena25, ena26, ena27, 
             ena28, ena29, ena30, ena31, ena32, ena33, ena34, top="LD r2 value between all pink panel markers in the even lineage, North American populations; positions based on Chinook chromosomes")

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

#assign each of the Chromosomes a plot that will be called later in the grid arrange
oa0<- ggplot()+geom_point(data=LD_odd_a_chr0 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Unassigned Chr")
oa1<- ggplot()+geom_point(data=LD_odd_a_chr1 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 1")
oa2<- ggplot()+geom_point(data=LD_odd_a_chr2 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 2")
oa3<- ggplot()+geom_point(data=LD_odd_a_chr3 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 3")
oa4<- ggplot()+geom_point(data=LD_odd_a_chr4 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 4")
oa5<- ggplot()+geom_point(data=LD_odd_a_chr5 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 5")
oa6<- ggplot()+geom_point(data=LD_odd_a_chr6 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 6")
oa7<- ggplot()+geom_point(data=LD_odd_a_chr7 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 7")
oa8<- ggplot()+geom_point(data=LD_odd_a_chr8 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 8")
oa9<- ggplot()+geom_point(data=LD_odd_a_chr9 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 9")
oa10<- ggplot()+geom_point(data=LD_odd_a_chr10,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 10")
oa11<- ggplot()+geom_point(data=LD_odd_a_chr11,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 11")
oa12<- ggplot()+geom_point(data=LD_odd_a_chr12,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 12")
oa13<- ggplot()+geom_point(data=LD_odd_a_chr13,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 13")
oa14<- ggplot()+geom_point(data=LD_odd_a_chr14,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 14")
oa15<- ggplot()+geom_point(data=LD_odd_a_chr15,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 15")
oa16<- ggplot()+geom_point(data=LD_odd_a_chr16,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 16")
oa17<- ggplot()+geom_point(data=LD_odd_a_chr17,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 17")
oa18<- ggplot()+geom_point(data=LD_odd_a_chr18,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 18")
oa19<- ggplot()+geom_point(data=LD_odd_a_chr19,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 19")
oa20<- ggplot()+geom_point(data=LD_odd_a_chr20,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 20")
oa21<- ggplot()+geom_point(data=LD_odd_a_chr21,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 21")
oa22<- ggplot()+geom_point(data=LD_odd_a_chr22,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 22")
oa23<- ggplot()+geom_point(data=LD_odd_a_chr23,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 23")
oa24<- ggplot()+geom_point(data=LD_odd_a_chr24,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 24")
oa25<- ggplot()+geom_point(data=LD_odd_a_chr25,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 25")
oa26<- ggplot()+geom_point(data=LD_odd_a_chr26,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 26")
oa27<- ggplot()+geom_point(data=LD_odd_a_chr27,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 27")
oa28<- ggplot()+geom_point(data=LD_odd_a_chr28,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 28")
oa29<- ggplot()+geom_point(data=LD_odd_a_chr29,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 29")
oa30<- ggplot()+geom_point(data=LD_odd_a_chr30,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 30")
oa31<- ggplot()+geom_point(data=LD_odd_a_chr31,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 31")
oa32<- ggplot()+geom_point(data=LD_odd_a_chr32,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 32")
oa33<- ggplot()+geom_point(data=LD_odd_a_chr33,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 33")
oa34<- ggplot()+geom_point(data=LD_odd_a_chr34,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd Asia individuals, Chr. 34")

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/LD_r2_aLL_pink_panel_Odd_Asia_2312.pdf", width = 9, height = 7)

#list the plots and call them in their layout
grid.arrange(oa1, oa2, oa3, oa4, oa5, oa6, oa7, oa8, oa9, oa10, oa11, oa12, oa13, oa14, oa15, oa16, oa17, oa18, oa19, oa20, oa21, oa22, oa23, oa24, oa25, oa26, oa27, 
             oa28, oa29, oa30, oa31, oa32, oa33, oa34, top="LD r2 value between all pink panel markers in the odd lineage, Asian populations; positions based on Chinook chromosomes")

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
#assign each of the Chromosomes a plot that will be called later in the grid arrange
ona0<- ggplot()+geom_point(data=LD_odd_na_chr0 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Unassigned Chr")
ona1<- ggplot()+geom_point(data=LD_odd_na_chr1 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 1")
ona2<- ggplot()+geom_point(data=LD_odd_na_chr2 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 2")
ona3<- ggplot()+geom_point(data=LD_odd_na_chr3 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 3")
ona4<- ggplot()+geom_point(data=LD_odd_na_chr4 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 4")
ona5<- ggplot()+geom_point(data=LD_odd_na_chr5 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 5")
ona6<- ggplot()+geom_point(data=LD_odd_na_chr6 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 6")
ona7<- ggplot()+geom_point(data=LD_odd_na_chr7 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 7")
ona8<- ggplot()+geom_point(data=LD_odd_na_chr8 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 8")
ona9<- ggplot()+geom_point(data=LD_odd_na_chr9 ,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 9")
ona10<- ggplot()+geom_point(data=LD_odd_na_chr10,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 10")
ona11<- ggplot()+geom_point(data=LD_odd_na_chr11,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 11")
ona12<- ggplot()+geom_point(data=LD_odd_na_chr12,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 12")
ona13<- ggplot()+geom_point(data=LD_odd_na_chr13,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 13")
ona14<- ggplot()+geom_point(data=LD_odd_na_chr14,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 14")
ona15<- ggplot()+geom_point(data=LD_odd_na_chr15,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 15")
ona16<- ggplot()+geom_point(data=LD_odd_na_chr16,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 16")
ona17<- ggplot()+geom_point(data=LD_odd_na_chr17,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 17")
ona18<- ggplot()+geom_point(data=LD_odd_na_chr18,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 18")
ona19<- ggplot()+geom_point(data=LD_odd_na_chr19,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 19")
ona20<- ggplot()+geom_point(data=LD_odd_na_chr20,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 20")
ona21<- ggplot()+geom_point(data=LD_odd_na_chr21,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 21")
ona22<- ggplot()+geom_point(data=LD_odd_na_chr22,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 22")
ona23<- ggplot()+geom_point(data=LD_odd_na_chr23,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 23")
ona24<- ggplot()+geom_point(data=LD_odd_na_chr24,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 24")
ona25<- ggplot()+geom_point(data=LD_odd_na_chr25,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 25")
ona26<- ggplot()+geom_point(data=LD_odd_na_chr26,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 26")
ona27<- ggplot()+geom_point(data=LD_odd_na_chr27,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 27")
ona28<- ggplot()+geom_point(data=LD_odd_na_chr28,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 28")
ona29<- ggplot()+geom_point(data=LD_odd_na_chr29,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 29")
ona30<- ggplot()+geom_point(data=LD_odd_na_chr30,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 30")
ona31<- ggplot()+geom_point(data=LD_odd_na_chr31,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 31")
ona32<- ggplot()+geom_point(data=LD_odd_na_chr32,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 32")
ona33<- ggplot()+geom_point(data=LD_odd_na_chr33,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 33")
ona34<- ggplot()+geom_point(data=LD_odd_na_chr34,aes(x=BP_A,y=BP_B,color=R2),shape=15, alpha= .35)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(),  axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Odd NA individuals, Chr. 34")

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/LD_r2_aLL_pink_panel_Odd_NA_2312.pdf", width = 9, height = 7)
#list the plots and call them in their layout
grid.arrange(ona1, ona2, ona3, ona4, ona5, ona6, ona7, ona8, ona9, ona10, ona11, ona12, ona13, ona14, ona15, ona16, ona17, ona18, ona19, ona20, ona21, ona22, ona23, ona24, ona25, ona26, ona27, 
             ona28, ona29, ona30, ona31, ona32, ona33, ona34, top="LD r2 value between all pink panel markers in the odd lineage, North American populations; positions based on Chinook chromosomes")

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