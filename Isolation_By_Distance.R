### Isolation by Distance calculations for populations
###    Distances in Kilometers, from ArcGIS
### Carolyn Tarpey | March 2016
### ---------------------------------------

#this is the genepop data set, if you want to put it in adegenet
#pops_16681_edit <- read.genepop("G:/Analysis/Pop_analysis/Populations_b3_may/Adegenet/16681_80_edit.gen")


library(adegenet)
library(ggplot2)
library(ape)

setwd('G:/Analysis/Pop_analysis/Populations_b3_may/IBD')

#making matrices of 16681 Fst

odd16881 <- read.table ("16681_fst_MATRIX_odd.txt", sep = "\t")
odd16881 <- as.matrix(odd16881)  

even16881 <- read.table ("16681_fst_MATRIX_even.txt", sep = "\t")
even16881 <- as.matrix(even16881)  

all16881 <- read.table ("16681_fst_MATRIX.txt", sep = "\t")
all16881 <- as.matrix(all16881)  

#making matrices of 15996 Fst

odd15996 <- read.table ("15996_fst_MATRIX_odd.txt", sep = "\t")
odd15996 <- as.matrix(odd15996)  

even15996 <- read.table ("15996_fst_MATRIX_even.txt", sep = "\t")
even15996 <- as.matrix(even15996)  

all15996 <- read.table ("15996_fst_MATRIX.txt", sep = "\t")
all15996 <- as.matrix(all15996)  

#making geography matrices

allgeo<- read.table ("geography_MATRIX.txt", sep = "\t")
allgeo <- as.matrix(allgeo)  

oddgeo <- read.table ("geography_MATRIX_odd.txt", sep = "\t")
oddgeo <- as.matrix(oddgeo)  

evengeo <- read.table ("geography_MATRIX_even.txt", sep = "\t")
evengeo <- as.matrix(evengeo)  


####mantel tests of the fst matrices against the geographic distance


odd16881Mantel <- mantel.test(odd16881, oddgeo, nperm = 999 )

even16881Mantel <- mantel.test(even16881, evengeo, nperm = 999)

all16881Mantel <- mantel.test(all16881, allgeo, nperm = 999)

odd15996Mantel <- mantel.test(odd15996, oddgeo, nperm = 999 )

even15996Mantel <- mantel.test(even15996, evengeo, nperm = 999)

all15996Mantel <- mantel.test(all15996, allgeo, nperm = 999)

# Results

odd16881Mantel 
even16881Mantel
all16881Mantel

odd15996Mantel
even15996Mantel
all15996Mantel
