### Create 6 lists of which samples to keep for PLINK
###    Splits by lineage then by region within lineage
###    
### Carolyn Tarpey | May 2018
### ---------------------------------------

#It requires the list of the populations and idividuals that are listed in the PED file of PLINK. I wrote a little python code that takes the first two columns from that file
# Subset_PED_for_KEEP.py 

install.packages("RColorBrewer")
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
library(colorRampPalette)
library(RColorBrewer)

#load the txt file that has the pops and the individuals in each: 

PLINK_PED <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/makingPlinkFiles/Panel_2312_popInds.txt",sep="", header = TRUE)
head(PLINK_PED)

##### Lists of the populations in each of the 6 groupings: 
e_list <- c("pop_1", "pop_4", "pop_6", "pop_7", "pop_16", "pop_9", "pop_12", "pop_14")
o_list <- c("pop_2", "pop_3", "pop_5", "pop_8", "pop_15", "pop_10", "pop_11", "pop_13")
e_na_list <- c("pop_9", "pop_12", "pop_14")
e_a_list <- c("pop_1", "pop_4", "pop_6", "pop_7", "pop_16")
o_na_list <- c("pop_10", "pop_11", "pop_13")
o_a_list <- c("pop_2", "pop_3", "pop_5", "pop_8", "pop_15")

#Subset the list of all the populations and individuals by the 6 groups based on population
Even_inds <- PLINK_PED[PLINK_PED$Pop%in%e_list,]
Odd_inds <- PLINK_PED[PLINK_PED$Pop%in%o_list,]
e_na_inds <- PLINK_PED[PLINK_PED$Pop%in%e_na_list,]
e_a_inds <- PLINK_PED[PLINK_PED$Pop%in%e_a_list,]
o_na_inds <- PLINK_PED[PLINK_PED$Pop%in%o_na_list,]
o_a_inds <- PLINK_PED[PLINK_PED$Pop%in%o_a_list,]


####Write The lists of the individuals in each population group to a file
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/Even_inds.txt", "wb")
write.table(Even_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/Odd_inds.txt", "wb")
write.table(Odd_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/e_na_inds.txt", "wb")
write.table(e_na_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/e_a_inds.txt", "wb")
write.table(e_a_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/o_na_inds.txt", "wb")
write.table(o_na_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/o_a_inds.txt", "wb")
write.table(o_a_inds,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)



#Plot the difference in the population group numbers 
#summary of all the idividuals by population
Inds_by_pop<- table(PLINK_PED[,1])
Inds_by_pop <- as.data.frame(Inds_by_pop)
colnames(Inds_by_pop) <- c("PopNum","Ind_counts")
Inds_by_pop$Pop <- c("Amur10", "Lakel07", "Nome91","Nome94","Snoh03","Snoh96","Tauy09","Tauy12","Amur11", "Hayly09", "Hayly10", "Koppe91", "Koppe96", "Kushi06", "Kushi07", "Lakel06")
Inds_by_pop$Lineage <- c("Even", "Odd", "Odd", "Even", "Odd", "Even", "Odd", "Even", "Odd", "Odd", "Even", "Odd", "Even", "Even","Odd", "Even")
Inds_by_pop <- as.data.frame(Inds_by_pop)

#summary of all the individuals by the 6 groupings
Ind_counts <- c(dim(Even_inds)[1], dim(Odd_inds)[1], dim(e_na_inds)[1], dim(e_a_inds)[1], dim(o_na_inds)[1], dim(o_a_inds)[1])
Pop_groups <- c("Even", "Odd", "Even_NA", "Even_A", "Odd_NA", "Odd_A")
Ind_counts <- as.numeric(Ind_counts)
Summary_Ped <- data.frame(Ind_counts, Pop_groups)

#colors
Histo_14_colors<- c("#6e016b","#6e016b","#88419d","#88419d","#8c6bb1","#8c6bb1","#8c96c6","#8c96c6","#7bccc4","#7bccc4",
                    '#4eb3d3','#4eb3d3','#2b8cbe','#2b8cbe')
plot(rep(1,14),col=Histo_14_colors,pch=19,cex=7)

Histo_16_colors<- c("#6e016b","#6e016b","#88419d","#88419d","#8c6bb1","#8c6bb1","#8c96c6","#8c96c6","#83B1C5", "#83B1C5",
                    "#7bccc4","#7bccc4", '#4eb3d3','#4eb3d3','#2b8cbe','#2b8cbe')
plot(rep(1,16),col=Histo_16_colors,pch=19,cex=7)


Blue_pink_6 <- c("Even"="#000080", "Odd"="#FF00FF", "Even_NA"="#AAAAD4", "Even_A"= "#5555AA", "Odd_NA"="#FFAAFF", "Odd_A" = "#FF55FF")
plot(rep(1,6),col=Blue_pink_6,pch=19,cex=7)

#plot of all the individuals by the 6 groupings
ggplot(data=Summary_Ped, aes(Summary_Ped$Pop_groups, Summary_Ped$Ind_counts)) + theme_bw() +
  geom_bar(aes(fill= Summary_Ped$Pop_groups),  stat="identity",alpha= .75) + theme(legend.position="none") +
  xlab("Population Group") + ylab("Count of Individuals") + ggtitle("Number of Individuals per Group") +
  scale_fill_manual(values = Blue_pink_6, labels = c("Even", "Odd", "Even_NA", "Even_A", "Odd_NA", "Odd_A")) 


#plot of all the individuals by their populations
ggplot(data=Inds_by_pop, aes(Inds_by_pop$Lineage, Inds_by_pop$Ind_counts)) + theme_bw() +
  geom_bar(aes(fill= Inds_by_pop$Pop), stat="identity", position="dodge", alpha= .75) + theme(legend.position="none") +
  xlab("Population") + ylab("Count of Individuals") + ggtitle("Number of Individuals per Population") +
  scale_fill_manual(values = Histo_16_colors) + geom_text(aes(x = Inds_by_pop$Lineage, y = Inds_by_pop$Ind_counts, 
          label = Inds_by_pop$Pop, group = Inds_by_pop$Pop),position = position_dodge(width = .9), vjust = -0.5, size = 2)



