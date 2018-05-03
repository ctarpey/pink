### Plotting the LD r2 from PLINK 
###    for the pink panel of 2312 markers
###    
### Garrett McKinney and Carolyn Tarpey | May 2018
### ---------------------------------------

#This R code was originally written by Garrett McKinney to plot the output of PLINK LD r2 for each chromosome. 
#It requires the LD output from PLINK. We used the alignment of the pinks to Chinook to get the Chromosome 
#and position assignments here, so there are 34 Chromosomes, and chromosome 35 are those that aligned to the unassigned 


# #This is the example code that garrett Sent me:
# ggplot()+geom_point(data=LDdata,aes(x=BP_A,y=BP_B,color=R2),shape=15)+scale_color_gradient(low="#E2E2E2",high="dark red")+theme_bw()
# ggplot()+geom_point(data=LDdata[LDdata$R2>=0.3,],aes(x=BP_A,y=BP_B,color=R2),shape=15)+scale_color_gradient(low="black",high="red")+theme_bw()

#Pink data filtering
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

#load genepop files
LDdata<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/Panel_2312_PLINK_out.ld",sep="")
head(LDdata)
dim(LDdata)

#Parse the data by chromosome- we want each chromosome by itself. 
LDdata_chr0 <- LDdata[which(LDdata$CHR_A == 0),]
LDdata_chr1 <- LDdata[which(LDdata$CHR_A == 1),]
LDdata_chr2 <- LDdata[which(LDdata$CHR_A == 2),]
LDdata_chr3 <- LDdata[which(LDdata$CHR_A == 3),]
LDdata_chr4 <- LDdata[which(LDdata$CHR_A == 4),]
LDdata_chr5 <- LDdata[which(LDdata$CHR_A == 5),]
LDdata_chr6 <- LDdata[which(LDdata$CHR_A == 6),]
LDdata_chr7 <- LDdata[which(LDdata$CHR_A == 7),]
LDdata_chr8 <- LDdata[which(LDdata$CHR_A == 8),]
LDdata_chr9 <- LDdata[which(LDdata$CHR_A == 9),]
LDdata_chr10 <- LDdata[which(LDdata$CHR_A == 10),]
LDdata_chr11 <- LDdata[which(LDdata$CHR_A == 11),]
LDdata_chr12 <- LDdata[which(LDdata$CHR_A == 12),]
LDdata_chr13 <- LDdata[which(LDdata$CHR_A == 13),]
LDdata_chr14 <- LDdata[which(LDdata$CHR_A == 14),]
LDdata_chr15 <- LDdata[which(LDdata$CHR_A == 15),]
LDdata_chr16 <- LDdata[which(LDdata$CHR_A == 16),]
LDdata_chr17 <- LDdata[which(LDdata$CHR_A == 17),]
LDdata_chr18 <- LDdata[which(LDdata$CHR_A == 18),]
LDdata_chr19 <- LDdata[which(LDdata$CHR_A == 19),]
LDdata_chr20 <- LDdata[which(LDdata$CHR_A == 20),]
LDdata_chr21 <- LDdata[which(LDdata$CHR_A == 21),]
LDdata_chr22 <- LDdata[which(LDdata$CHR_A == 22),]
LDdata_chr23 <- LDdata[which(LDdata$CHR_A == 23),]
LDdata_chr24 <- LDdata[which(LDdata$CHR_A == 24),]
LDdata_chr25 <- LDdata[which(LDdata$CHR_A == 25),]
LDdata_chr26 <- LDdata[which(LDdata$CHR_A == 26),]
LDdata_chr27 <- LDdata[which(LDdata$CHR_A == 27),]
LDdata_chr28 <- LDdata[which(LDdata$CHR_A == 28),]
LDdata_chr29 <- LDdata[which(LDdata$CHR_A == 29),]
LDdata_chr30 <- LDdata[which(LDdata$CHR_A == 30),]
LDdata_chr31 <- LDdata[which(LDdata$CHR_A == 31),]
LDdata_chr32 <- LDdata[which(LDdata$CHR_A == 32),]
LDdata_chr33 <- LDdata[which(LDdata$CHR_A == 33),]
LDdata_chr34 <- LDdata[which(LDdata$CHR_A == 34),]
LDdata_chr35 <- LDdata[which(LDdata$CHR_A == 35),]


#assign each of the Chromosomes a plot that will be called later in the grid arrange 
c0<- ggplot()+geom_point(data=LDdata_chr0,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("No matching Alignments") 
c1<- ggplot()+geom_point(data=LDdata_chr1,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 1") 
c2<- ggplot()+geom_point(data=LDdata_chr2,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 2") 
c3<- ggplot()+geom_point(data=LDdata_chr3,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 3") 
c4<- ggplot()+geom_point(data=LDdata_chr4,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 4") 
c5<- ggplot()+geom_point(data=LDdata_chr5,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 5") 
c6<- ggplot()+geom_point(data=LDdata_chr6,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 6") 
c7<- ggplot()+geom_point(data=LDdata_chr7,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 7") 
c8<- ggplot()+geom_point(data=LDdata_chr8,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 8") 
c9<- ggplot()+geom_point(data=LDdata_chr9,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 9") 
c10<- ggplot()+geom_point(data=LDdata_chr10,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 10") 
c11<- ggplot()+geom_point(data=LDdata_chr11,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 11") 
c12<- ggplot()+geom_point(data=LDdata_chr12,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 12") 
c13<- ggplot()+geom_point(data=LDdata_chr13,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 13") 
c14<- ggplot()+geom_point(data=LDdata_chr14,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 14") 
c15<- ggplot()+geom_point(data=LDdata_chr15,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 15") 
c16<- ggplot()+geom_point(data=LDdata_chr16,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 16") 
c17<- ggplot()+geom_point(data=LDdata_chr17,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 17") 
c18<- ggplot()+geom_point(data=LDdata_chr18,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 18") 
c19<- ggplot()+geom_point(data=LDdata_chr19,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 19") 
c20<- ggplot()+geom_point(data=LDdata_chr20,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 20") 
c21<- ggplot()+geom_point(data=LDdata_chr21,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 21") 
c22<- ggplot()+geom_point(data=LDdata_chr22,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 22") 
c23<- ggplot()+geom_point(data=LDdata_chr23,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 23") 
c24<- ggplot()+geom_point(data=LDdata_chr24,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 24") 
c25<- ggplot()+geom_point(data=LDdata_chr25,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 25") 
c26<- ggplot()+geom_point(data=LDdata_chr26,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 26") 
c27<- ggplot()+geom_point(data=LDdata_chr27,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 27") 
c28<- ggplot()+geom_point(data=LDdata_chr28,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 28") 
c29<- ggplot()+geom_point(data=LDdata_chr29,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 29") 
c30<- ggplot()+geom_point(data=LDdata_chr30,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 30") 
c31<- ggplot()+geom_point(data=LDdata_chr31,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 31") 
c32<- ggplot()+geom_point(data=LDdata_chr32,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 32") 
c33<- ggplot()+geom_point(data=LDdata_chr33,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 33") 
c34<- ggplot()+geom_point(data=LDdata_chr34,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Chr. 34") 
c35<- ggplot()+geom_point(data=LDdata_chr35,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")+ ggtitle("Unassigned") 


#list the plots and call them in their layout
grid.arrange(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27, 
             c28, c29, c30, c31, c32, c33, c34, c35, nrow=5, top="LD r2 value between pink panel markers assigned to same Chinook chromosome")

#this is supposed to give them all one legend, and it does, but IDK how to make it be in the correct spot. 
grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + theme(legend.position="right"))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  grid.arrange(
    do.call(arrangeGrob, lapply(plots, function(x)
      x + theme(legend.position="none",
                strip.background = element_blank(),
                strip.text.x = element_blank()))),
    legend,
    ncol = 1,
    heights = unit.c(unit(1, "npc") - lheight, lheight))
}

grid_arrange_shared_legend(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27, 
             c28, c29, c30, c31, c32, c33, c34, c35, nrow=5, top="LD r2 value between pink panel markers assigned to same Chinook chromosome")


## plotting each of the chromosomes individually so that you can flip through and understand what is happening in each. 


#assign each of the Chromosomes a plot that will be called later in the grid arrange 
ggplot()+geom_point(data=LDdata_chr0,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("No matching Alignments") 
ggplot()+geom_point(data=LDdata_chr1,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 1") 
ggplot()+geom_point(data=LDdata_chr2,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 2") 
ggplot()+geom_point(data=LDdata_chr3,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 3") 
ggplot()+geom_point(data=LDdata_chr4,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 4") 
ggplot()+geom_point(data=LDdata_chr5,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 5") 
ggplot()+geom_point(data=LDdata_chr6,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 6") 
ggplot()+geom_point(data=LDdata_chr7,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 7") 
ggplot()+geom_point(data=LDdata_chr8,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 8") 
ggplot()+geom_point(data=LDdata_chr9,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 9") 
ggplot()+geom_point(data=LDdata_chr10,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 10") 
ggplot()+geom_point(data=LDdata_chr11,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 11") 
ggplot()+geom_point(data=LDdata_chr12,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 12") 
ggplot()+geom_point(data=LDdata_chr13,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 13") 
ggplot()+geom_point(data=LDdata_chr14,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 14") 
ggplot()+geom_point(data=LDdata_chr15,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 15") 
ggplot()+geom_point(data=LDdata_chr16,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 16") 
ggplot()+geom_point(data=LDdata_chr17,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 17") 
ggplot()+geom_point(data=LDdata_chr18,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 18") 
ggplot()+geom_point(data=LDdata_chr19,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 19") 
ggplot()+geom_point(data=LDdata_chr20,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 20") 
ggplot()+geom_point(data=LDdata_chr21,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 21") 
ggplot()+geom_point(data=LDdata_chr22,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 22") 
ggplot()+geom_point(data=LDdata_chr23,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 23") 
ggplot()+geom_point(data=LDdata_chr24,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 24") 
ggplot()+geom_point(data=LDdata_chr25,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 25") 
ggplot()+geom_point(data=LDdata_chr26,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 26") 
ggplot()+geom_point(data=LDdata_chr27,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 27") 
ggplot()+geom_point(data=LDdata_chr28,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 28") 
ggplot()+geom_point(data=LDdata_chr29,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 29") 
ggplot()+geom_point(data=LDdata_chr30,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 30") 
ggplot()+geom_point(data=LDdata_chr31,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 31") 
ggplot()+geom_point(data=LDdata_chr32,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 32") 
ggplot()+geom_point(data=LDdata_chr33,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 33") 
ggplot()+geom_point(data=LDdata_chr34,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Chr. 34") 
ggplot()+geom_point(data=LDdata_chr35,aes(x=BP_A,y=BP_B,color=R2),shape=15)+ labs(x="", y="") +scale_color_gradient(low="#E2E2E2",high="dark red")+ theme_bw()+  ggtitle("Unassigned") 

