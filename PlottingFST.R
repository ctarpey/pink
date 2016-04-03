### plotting Fst heatmap
###    Fst from genepop 16681 is all 15996 are nondups
### Carolyn Tarpey | March 2016
### ---------------------------------------

#install.packages("ggplot2")
#install.packages("reshape2")
#install.packages("stringi")
install.packages("clusterGenomics")

library(adegenet)
library(ggplot2)
library(ape)
library(reshape2)
library(clusterGenomics)


setwd('G:/Analysis/Pop_analysis/Populations_b3_may/PlottingFST')


#making dataframes of 16681 Fst

odd16881 <- read.table ("16681_fst_odd_ordered.txt", sep = "\t", header = TRUE)


head(odd16881)

row.names(odd16881) <- odd16881$X16681_odd
odd16881.n <-  odd16881[,2:8]

odd16881.matrix <- data.matrix(odd16881.n)

odd_heatmap <- heatmap(odd16881.matrix, Rowv=NA, Colv=NA, col = cm.colors(256), scale="column", margins=c(5,10))






odd16881 <- as.matrix(odd16881)

heatmap(odd16881, Colv= FALSE, scale = 'none')

melted_16681_odd <- melt(odd16881)
head(melted_16681_odd)

ggplot(data = melted_16681_odd, aes(x=X16681_odd, y=variable, fill=value)) +  geom_tile() +
  scale_fill_gradient2(low = "white", high = "steelblue")
                       
                       
#set the row names to row names
rownames(odd16881)<- odd16881[,1]
odd16881.n <-  odd16881[,2:8]

# Get lower triangle of the correlation matrix
get_lower_tri<-function(x){
  x[upper.tri(x)] <- NA
  return(x)
}

# Get upper triangle of the correlation matrix
get_upper_tri <- function(x){
  x[lower.tri(x)]<- NA
  return(x)
}

#return the upper and lower triangles

upper_tri <- get_upper_tri(odd16881.n)
upper_tri

lower_tri <-get_lower_tri(odd16881.n)
lower_tri

melted_16681_odd <- melt(lower_tri, na.rm = TRUE)
head(melted_16681_odd)

ggplot(data = melted_16681_odd, aes(x=melted_16681_odd, y=variable, fill=value)) +  geom_tile() +
  scale_fill_gradient2(low = "white", high = "steelblue")

plotHeatmap(odd16881.n, colorscale = "blue-white-red")





even16881 <- read.table ("16681_fst_MATRIX_even.txt", sep = "\t")
even16881 <- as.data.frame(even16881)  



#all 

all16881 <- read.table ("16681_fst_MATRIX.txt", sep = "\t", header= TRUE)
#all16881 <- as.data.frame(all16881)  

melted_16681 <- melt(all16881)
head(melted_16681)
ggplot(data = melted_16681) + geom_tile(aes(x=X16681, y=variable, fill=value))



ggplot(data = melted_16681, aes(x=X16681, y=variable, fill=value)) +  geom_tile() +
  scale_fill_gradient2(low = "white", high = "steelblue")










###########no duplicates

# 
# 
# #making dataframes of 15996 Fst
# 
# odd15996 <- read.table ("15996_fst_MATRIX_odd.txt", sep = "\t")
# odd15996 <- as.data.frame(odd15996)  
# 
# even15996 <- read.table ("15996_fst_MATRIX_even.txt", sep = "\t")
# even15996 <- as.data.frame(even15996)  
# 
# all15996 <- read.table ("15996_fst_MATRIX.txt", sep = "\t")
# all15996 <- as.data.frame(all15996)  
