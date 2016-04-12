### Plotting the outliers using a Circos plot
###  
### Carolyn Tarpey | April 2016 
### ---------------------------------------

library("RCircos")

source("https://bioconductor.org/biocLite.R")
biocLite("OmicCircos")
library("OmicCircos")

source("https://bioconductor.org/biocLite.R")
biocLite("ggbio")

library("RColorBrewer")
library(ggplot2)


setwd('G:/Analysis/Mapping/Outliers_Manhattan_Plot/Arlequin_LFMM')

par()  ##view default par settings
opar<-par()   ##save a copy of default par settings so that I can revert back to them while plotting

# Import Fst values and save values for each pop (example of full code)
Fst_mapped <-read.table("Outliers_Master.txt", header=TRUE) ##This file has cumulative map positions
head(Fst_mapped)