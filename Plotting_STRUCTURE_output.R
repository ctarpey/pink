### Plotting STRUCTURE output with the program STRUCTURE PLOT
###    runs of STRUCTURE with k 1: 18 
### Carolyn Tarpey | December  2015 
### ---------------------------------------

library(shiny)
library(markdown)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

#the input file for the plot is the meanQ with a header line and two additional
#columns at the beginning of the file, the first one is the population name
#the second is the name of the sample- the whole thing is a txt file


setwd("G:/Analysis/Pop_analysis/Populations_b3_may/fastSTRUCTURE/out" )

runApp(appDir = getwd(), port = getOption("shiny.port"))


