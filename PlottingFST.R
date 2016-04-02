### plotting Fst heatmap
###    Fst from genepop 16681 is all 15996 are nondups
### Carolyn Tarpey | March 2016
### ---------------------------------------


library(adegenet)
library(ggplot2)
library(ape)
library(reshape2)


setwd('G:/Analysis/Pop_analysis/Populations_b3_may/PlottingFST')


#making dataframes of 16681 Fst

odd16881 <- read.table ("16681_fst_MATRIX_odd.txt", sep = "\t", header = TRUE)

melted_16681odd <- melt(odd16881)
head(melted_16681odd)
ggplot(data = melted_16681odd) + geom_tile(aes(x=X16681odd, y=variable, fill=value))


ggplot(data = melted_16681odd, aes(x=X16681odd, y=variable, fill=value)) +  geom_tile() +
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

melted_16681_odd <- melt(upper_tri, na.rm = TRUE)


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







#making dataframes of 15996 Fst

odd15996 <- read.table ("15996_fst_MATRIX_odd.txt", sep = "\t")
odd15996 <- as.data.frame(odd15996)  

even15996 <- read.table ("15996_fst_MATRIX_even.txt", sep = "\t")
even15996 <- as.data.frame(even15996)  

all15996 <- read.table ("15996_fst_MATRIX.txt", sep = "\t")
all15996 <- as.data.frame(all15996)  
