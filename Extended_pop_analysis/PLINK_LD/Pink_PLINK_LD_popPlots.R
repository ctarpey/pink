### Create 6 lists of which samples to keep for PLINK
###    Splits by lineage then by region within lineage
###    Plots the differences in the numbers between groups
### Carolyn Tarpey | May 2018 #updated may 2018 to make sure Lakel07 is even and lakel06 is odd!
### ---------------------------------------

#It requires the list of the populations and idividuals that are listed in the PED file of PLINK. I wrote a little python code that takes the first two columns from that file
# Subset_PED_for_KEEP.py 

#install.packages("gridExtra")
#library(colorRampPalette)
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

#load the txt file that has the pops and the individuals in each: 
PLINK_PED <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/PLINK_23759_465_POP_INDS.txt",sep="", header = FALSE)
colnames(PLINK_PED)<- c("Pop", "Inds")
head(PLINK_PED)

POP_INFO <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/PlottingLD_pop_info.txt",sep="", header =TRUE)
head(POP_INFO)



##### Lists of the populations in each of the 6 groupings: 
e_list <- POP_INFO[POP_INFO$LINEAGE=="Even",]
o_list<- POP_INFO[POP_INFO$LINEAGE=="Odd",]
e_na_list<- POP_INFO[(POP_INFO$LINEAGE=="Even") & (POP_INFO$CONTINENT=="North_America"),]
e_a_list <- POP_INFO[(POP_INFO$LINEAGE=="Even") & (POP_INFO$CONTINENT=="Asia"),]
o_na_list<- POP_INFO[(POP_INFO$LINEAGE=="Odd") & (POP_INFO$CONTINENT=="North_America"),]
o_a_list <- POP_INFO[(POP_INFO$LINEAGE=="Odd") & (POP_INFO$CONTINENT=="Asia"),]

e_NS_list <- POP_INFO[(POP_INFO$LINEAGE=="Even") & (POP_INFO$POPNAME !="SUSIT14"),]
o_NS_list <- POP_INFO[(POP_INFO$LINEAGE=="Odd") & (POP_INFO$POPNAME !="SUSIT13"),]
e_na_NS_list <- POP_INFO[(POP_INFO$LINEAGE=="Even") & (POP_INFO$CONTINENT=="North_America") & (POP_INFO$POPNAME !="SUSIT14"),]
o_na_NS_list<- POP_INFO[(POP_INFO$LINEAGE=="Odd") & (POP_INFO$CONTINENT=="North_America") & (POP_INFO$POPNAME !="SUSIT13"),]

na_all_list <- POP_INFO[POP_INFO$CONTINENT =="North_America",]
na_NS_list <- POP_INFO[(POP_INFO$CONTINENT =="North_America") & (POP_INFO$POPNAME !="SUSIT13" & POP_INFO$POPNAME !="SUSIT14" ),]
a_all_list <- POP_INFO[POP_INFO$CONTINENT =="Asia",]


#Subset the list of all the populations and individuals by the 6 groups based on population
Even_inds <- PLINK_PED[PLINK_PED$Pop%in%e_list$POP,]
Odd_inds <- PLINK_PED[PLINK_PED$Pop%in%o_list$POP,]
Even_NS_inds <- PLINK_PED[PLINK_PED$Pop%in%e_NS_list$POP,]
Odd_NS_inds <- PLINK_PED[PLINK_PED$Pop%in%o_NS_list$POP,]

e_na_inds <- PLINK_PED[PLINK_PED$Pop%in%e_na_list$POP,]
e_na_NS_inds <- PLINK_PED[PLINK_PED$Pop%in%e_na_NS_list$POP,]
e_a_inds <- PLINK_PED[PLINK_PED$Pop%in%e_a_list$POP,]

o_na_inds <- PLINK_PED[PLINK_PED$Pop%in%o_na_list$POP,]
o_na_NS_inds <- PLINK_PED[PLINK_PED$Pop%in%o_na_NS_list$POP,]
o_a_inds <- PLINK_PED[PLINK_PED$Pop%in%o_a_list$POP,]

na_all_inds<- PLINK_PED[PLINK_PED$Pop%in%na_all_list$POP,]
na_NS_inds <- PLINK_PED[PLINK_PED$Pop%in%na_NS_list$POP,]
a_all_inds <- PLINK_PED[PLINK_PED$Pop%in%a_all_list$POP,]

####Write The lists of the individuals in each population group to a file
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Even_inds_all.txt", "wb")
write.table(Even_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Odd_inds_all.txt", "wb")
write.table(Odd_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Even_inds_NS.txt", "wb")
write.table(Even_NS_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Odd_inds_NS.txt", "wb")
write.table(Odd_NS_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/e_na_inds_all.txt", "wb")
write.table(e_na_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/e_a_inds.txt", "wb")
write.table(e_a_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/o_na_inds_all.txt", "wb")
write.table(o_na_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/o_a_inds.txt", "wb")
write.table(o_a_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/e_na_inds_NS.txt", "wb")
write.table(e_na_NS_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/o_na_inds_NS.txt", "wb")
write.table(o_na_NS_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/na_all_inds.txt", "wb")
write.table(na_all_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/na_NS_inds.txt", "wb")
write.table(na_NS_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/a_all_inds.txt", "wb")
write.table(a_all_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)



#Plot the difference in the population group numbers 
#summary of all the individuals by population
Inds_by_pop<- table(PLINK_PED[,1])
Inds_by_pop <- as.data.frame(Inds_by_pop)
colnames(Inds_by_pop) <- c("POP","Ind_counts")
POP_counts <- merge(Inds_by_pop, POP_INFO, by = "POP")
head(POP_counts)
dim(POP_counts)

#summary of all the individuals by the 10 groupings
Ind_counts <- c(dim(Even_inds)[1], dim(Odd_inds)[1], dim(Even_NS_inds)[1], dim(Odd_NS_inds )[1], dim(e_na_inds)[1], dim(e_a_inds)[1], 
                dim(o_na_inds)[1], dim(o_a_inds)[1], dim(e_na_NS_inds)[1], dim(o_na_NS_inds)[1])
Pop_groups <- c("Even_inds", "Odd_inds", "Even_NS_inds", "Odd_NS_inds" , "e_na_inds", "e_a_inds", 
                "o_na_inds", "o_a_inds", "e_na_NS_inds", "o_na_NS_inds")
Lineage <- c("Even", "Odd", "Even", "Odd" , "Even", "Even", "Odd", "Odd", "Even", "Odd")
Region <-c("All", "All", "All", "All" , "America", "Asia", "America", "Asia", "America", "America")
Plotting <- c(1 ,6 ,2 ,7, 3 ,5 ,8 ,10 ,4 ,9 )
Ind_counts <- as.numeric(Ind_counts)
Summary_Ped <- data.frame(Ind_counts, Pop_groups, Lineage, Region,Plotting)
head(Summary_Ped)
Summary_Ped <-Summary_Ped[order(Lineage,Region),]
 

#colors
Histo_14_colors<- c("#6e016b","#6e016b","#88419d","#88419d","#8c6bb1","#8c6bb1","#8c96c6","#8c96c6","#7bccc4","#7bccc4",
                    '#4eb3d3','#4eb3d3','#2b8cbe','#2b8cbe')
plot(rep(1,14),col=Histo_14_colors,pch=19,cex=7)

Histo_16_colors<- c("#6e016b","#6e016b","#88419d","#88419d","#8c6bb1","#8c6bb1","#8c96c6","#8c96c6","#83B1C5", "#83B1C5",
                    "#7bccc4","#7bccc4", '#4eb3d3','#4eb3d3','#2b8cbe','#2b8cbe')
plot(rep(1,16),col=Histo_16_colors,pch=19,cex=7)

Histo_18_colors<- c("#6e016b","#6e016b","#88419d","#88419d","#8c6bb1","#8c6bb1","#8c96c6","#8c96c6","#83B1C5", "#83B1C5",
                    "#7bccc4","#7bccc4", '#4eb3d3','#4eb3d3','#2b8cbe','#2b8cbe','#1b5b7c','#1b5b7c')
plot(rep(1,18),col=Histo_18_colors,pch=19,cex=7)

Blue_pink_6 <- c("Even"="#000080", "Odd"="#FF00FF", "Even_NA"="#AAAAD4", "Even_A"= "#5555AA", "Odd_NA"="#FFAAFF", "Odd_A" = "#FF55FF")
plot(rep(1,6),col=Blue_pink_6,pch=19,cex=7)

Blue_pink_7 <-c("ALL" = "#8800C3", "Even"="#000080", "Odd"="#FF00FF", "Even_NA"="#AAAAD4", "Even_A"= "#5555AA", "Odd_NA"="#FFAAFF", 
                "Odd_A" = "#FF55FF")
plot(rep(1,7),col=Blue_pink_7,pch=19,cex=7)

Blue_pink_8 <-c("ALL" = "#57007c", "ALL_NS" = "#8800C3", "Even"="#000080", "Odd"="#FF00FF", "Even_NA"="#AAAAD4", "Even_A"= "#5555AA", 
                "Odd_NA"="#FFAAFF", "Odd_A" = "#FF55FF")
plot(rep(1,8),col=Blue_pink_8,pch=19,cex=7)

Blue_pink_10 <-c("Even_inds"="#000080", "Even_NS_inds"="#3434af", "Odd_inds"="#FF00FF", "Odd_NS_inds"="#e067e0",
                 "e_na_inds"="#8b8baf", "e_na_NS_inds"="#AAAAD4", 
                 "e_a_inds"= "#5555AA", "o_na_inds"="#FF55FF", "o_na_NS_inds"= "#FFAAFF", "o_a_inds" = "#cc88cc")
plot(rep(1,10),col=Blue_pink_10,pch=19,cex=7)

Blue_pink_12 <-c("ALL" = "#57007c", "ALL_NS" = "#8800C3", "Even"="#000080", "Even_NS"="#3434af", "Odd"="#FF00FF", "Odd_NS"="#e067e0",
                 "Even_NA"="#8b8baf", "Even_NA_NS"="#AAAAD4", 
                 "Even_A"= "#5555AA", "Odd_NA"="#FF55FF", "Odd_NA_NS"= "#FFAAFF", "Odd_A" = "#cc88cc")
plot(rep(1,12),col=Blue_pink_12,pch=19,cex=7)

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/Population_group_counts.pdf", width = 9, height = 7)

#plot of all the individuals by the  groupings
ggplot(data=Summary_Ped, aes(Summary_Ped$Plotting, Summary_Ped$Ind_counts)) + theme_bw() +
  geom_bar(aes(fill= Summary_Ped$Pop_groups),  stat="identity", alpha= .75) + theme(legend.position="none") +
  xlab("Population Group") + ylab("Count of Individuals") + ggtitle("Number of Individuals per Group") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +  scale_fill_manual(values = Blue_pink_10) +
  scale_x_discrete(limits = c("Even", "Even_NS","Even_NA","Even_NA_NS", "Even_A","Odd", "Odd_NS", "Odd_NA", "Odd_NA_NS", "Odd_A" ))

#plot of all the individuals by their populations
ggplot(data=POP_counts, aes(POP_counts$LINEAGE, POP_counts$Ind_counts)) + theme_bw() +
  geom_bar(aes(fill= POP_counts$POP), stat="identity", position="dodge", alpha= .75) + theme(legend.position="none") +
  xlab("Population") + ylab("Count of Individuals") + ggtitle("Number of Individuals per Population") +
  scale_fill_manual(values = Histo_18_colors) + geom_text(aes(x = POP_counts$LINEAGE, y = POP_counts$Ind_counts, 
          label = POP_counts$POPNAME, group = POP_counts$POP),position = position_dodge(width = .9), vjust = -0.5, size = 2)

dev.off()
