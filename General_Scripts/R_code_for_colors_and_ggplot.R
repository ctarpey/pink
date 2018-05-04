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





#colors
Histo_14_colors<- c("#6e016b","#6e016b","#88419d","#88419d","#8c6bb1","#8c6bb1","#8c96c6","#8c96c6","#7bccc4","#7bccc4",'#4eb3d3','#4eb3d3','#2b8cbe','#2b8cbe')
plot(rep(1,6),col=Pop_groups_6,pch=19,cex=7)

Pop_groups_6<- c("#007adf",  "#6e84bf", "#498dbd","#a940c6", "#da0085", "#a26dd5")
plot(rep(1,6),col=Pop_groups_6,pch=19,cex=7)

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


