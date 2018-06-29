### Plotting code, especially for colors
###    
###   
###   Carolyn Tarpey | may 2018
### ---------------------------------------



library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)
library(plotly)
library(gridExtra)
library(scales) 

############ Importing 

#import the original dataset that has only been filtered for 80% genotype rate
OG_80PCfiltered <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/FirstFiltering/filteredGenos_just80PCgenorate_76762.txt",colClasses="factor", header = TRUE)
OG_80PCfiltered[1:5,1:5]
dim(OG_80PCfiltered)

All_col<- c("Hokkaido"="#61005e","Amur"="#cd5490","Magadan"="#e0957e","Kamchatka"="#c9cf4a","Norton Sound"="#677f3e","Cook Inlet" = "#563e7f", 
            "Prince William Sound"="#7297ee","Skeena" = "#2171b5", "Puget Sound"="#00158a")

##Show the colors used for a default palette of 18
show_col(hue_pal()(18))

#colors
Histo_14_colors<- c("#6e016b","#6e016b","#88419d","#88419d","#8c6bb1","#8c6bb1","#8c96c6","#8c96c6","#7bccc4","#7bccc4",'#4eb3d3','#4eb3d3','#2b8cbe','#2b8cbe')
plot(rep(1,14),col=Histo_14_colors,pch=19,cex=7)

Pop_groups_6<- c("#007adf",  "#6e84bf", "#498dbd","#a940c6", "#da0085", "#a26dd5")
plot(rep(1,6),col=Pop_groups_6,pch=19,cex=7)

colfunc <- colorRampPalette(c("#c9cf4a", "#2171b5"))
susitna_3 <- colfunc(4)
plot(rep(1,4),col=susitna_3,pch=19,cex=7)

colfunc <- colorRampPalette(c("navy", "magenta"))
colfunc(16)
x <- colfunc(16)
plot(rep(1,16),col=x,pch=19,cex=7)

colfunc <- colorRampPalette(c("white", "magenta"))
colfunc(4)
pink_3 <- colfunc(4)
plot(rep(1,4),col=pink_3,pch=19,cex=7)

colfunc <- colorRampPalette(c("white", "navy"))
colfunc(4)
blue_3<- colfunc(4)
plot(rep(1,4),col=blue_3,pch=19,cex=7)

Spectral_6<-brewer.pal(6, "Spectral")
plot(rep(1,6),col=Spectral_6,pch=19,cex=7)

Spectral_7<-brewer.pal(7, "Spectral")
plot(rep(1,7),col=Spectral_7,pch=19,cex=7)

#colors
Histo_14_colors<- c("#6e016b","#6e016b","#88419d","#88419d","#8c6bb1","#8c6bb1","#8c96c6","#8c96c6","#7bccc4","#7bccc4",
                    '#4eb3d3','#4eb3d3','#2b8cbe','#2b8cbe')
plot(rep(1,14),col=Histo_14_colors,pch=19,cex=7)

Histo_16_colors<- c("#6e016b","#6e016b","#88419d","#88419d","#8c6bb1","#8c6bb1","#8c96c6","#8c96c6","#83B1C5", "#83B1C5",
                    "#7bccc4","#7bccc4", '#4eb3d3','#4eb3d3','#2b8cbe','#2b8cbe')
plot(rep(1,16),col=Histo_16_colors,pch=19,cex=7)


Blue_pink_6 <- c("Even"="#000080", "Odd"="#FF00FF", "Even_NA"="#AAAAD4", "Even_A"= "#5555AA", "Odd_NA"="#FFAAFF", "Odd_A" = "#FF55FF")
plot(rep(1,6),col=Blue_pink_6,pch=19,cex=7)

######################################## Little bits of code
#this is to change the x axis so that the numbers are angled and not scientific method. 
require(grid) ?
theme(axis.text.x = element_text(angle = 45, hjust = 1))+ scale_x_continuous(labels = comma)

#to not include any of the x or y axis information- super sleek plots good for facet wrapping where the scales dont matter
+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position="none")

#for facet wrapping and having each plot be on its own scale
facet_wrap(~CHR_A, scale="free") 

#for subsetting a certain set of the data without having it explicitly pulled out already

ggplot(LDdata_all[(LDdata_all$CHR_A != "35") & (LDdata_all$R2 >=0.05),], aes(x=Rank, y= R2 )) + 
  geom_point(data= LDdata_all[(LDdata_all$CHR_A != "35") & (LDdata_all$R2 >=0.05),], colour = "#D53E4F", size =.5 ) + 

# for saving as a jpeg because the plot is too big to be rendered in pdf each time
  ggsave("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/LD/Ranked_LD_main_hierarchy_LDdata_all.jpg", 
     ggplot(LDdata_all[LDdata_all$CHR_A != "35",], aes(x=Rank, y= R2 )) + geom_point(data= LDdata_all[LDdata_all$CHR_A != "35",], colour = "black", size =.5 ) + 
     facet_wrap(~CHR_A, scale="free") + theme_bw(), width = 9,  height = 7,  dpi = 1200)

scale_x_continuous(limits = c(-5000, 5000))
+ xlim(-5000, 5000)

######################################## Full on plotting 
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

#plot of all the individuals by the 6 groupings
ggplot(data=Summary_Ped, aes(Summary_Ped$Pop_groups, Summary_Ped$Ind_counts)) + theme_bw() +
  geom_bar(aes(fill= Summary_Ped$Pop_groups),  stat="identity") + theme(legend.position="none") +
  xlab("Population Group") + ylab("Count of Individuals") + ggtitle("Number of Individuals per Group") +
  scale_fill_manual(values = Pop_groups_6, labels = c("Even", "Odd", "Even_NA", "Even_A", "Odd_NA", "Odd_A")) 


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


