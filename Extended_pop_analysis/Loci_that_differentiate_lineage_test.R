### Test for loci that differentiate Lineage 
###   To investigate the true lineage designation of Lakel samples
###   
###  Carolyn Tarpey | Dec 2017
### ---------------------------------------

#This is R code requires a genepop file for each of the genotype sets and should be stripped of the header line that Stacks puts in.  
#The second line should start with a tab, the loci names should be tab deliminated, and there should be no "pop" between the samples. 
#This formatting ensures R imports the genepop file as a table. 

#Pink data filtering
library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)
library(tidyr)


OLD_E_O <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/Lakel_lineage_test/OLD_E_O_allele_freq.txt", header= TRUE, na.strings = "NA")
NEW_E_O <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/Lakel_lineage_test/NEW_E_O_allele_freq.txt", header= TRUE, na.strings = "NA" )




#####################Upload the genotypes !!!!This takes up all the memory!!!!
#using the 16662 loci that were found between both old and new sets
# 
# ALL16662 <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/Lakel_lineage_test/batch_4_16662_all.genepop.txt",colClasses="factor")
# NEW16662 <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/Lakel_lineage_test/batch_4_16662_NEW.genepop.txt",colClasses="factor")
# OLD16662 <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/Lakel_lineage_test/batch_4_16662_old.genepop.txt",colClasses="factor")

popmap_all <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/Lakel_lineage_test/PopMap.txt")
popmap_old <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/Lakel_lineage_test/PopMap_old.txt")
popmap_new <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/Lakel_lineage_test/PopMap_new.txt")

ALL16662_pops<-cbind(popmap_all,ALL16662)
ALL16662_pops$V1<- NULL
colnames(ALL16662_pops)[1]<-"Pop"
ALL16662_pops[1:5,1:5]
length(ALL16662_pops)
dim(ALL16662_pops)

NEW16662_pops<-cbind(popmap_new,NEW16662)
NEW16662_pops$V1<- NULL
colnames(NEW16662_pops)[1]<-"Pop"
NEW16662_pops[1:5,1:5]
length(NEW16662_pops)
dim(NEW16662_pops)

OLD16662_pops<-cbind(popmap_old,OLD16662)
OLD16662_pops$V1<- NULL
colnames(OLD16662_pops)[1]<-"Pop"
OLD16662_pops[1:5,1:5]
length(OLD16662_pops)
dim(OLD16662_pops)
OLD16662_pops[1:5,1:5]

###Split each into Even and Odd Lineages

Lakel06 <- "PLAKEL06"
Lakel07 <- "PLAKEL07"

#OLD_16662_even_pops <- c("PAMUR10","PHAYLY10","PKOPE96","PKUSHI06","PNOME94","PSNOH96","PTAUY12")
#OLD_16662_odd_pops <-c("PAMUR11","PHAYLY09","PKOPE91","PKUSHI07","PNOME91","PSNOH03","PTAUY09")

NEW_16662_Lakel_even <-subset(NEW16662_pops, Pop == Lakel06)
NEW_16662_Lakel_even[1:5,1:5]
dim(NEW_16662_Lakel_even)

NEW_16662_Lakel_odd <-subset(NEW16662_pops, Pop == Lakel07)
NEW_16662_Lakel_odd[1:5,1:5]
dim(NEW_16662_Lakel_odd)

OLD_16662_even <-subset(OLD16662_pops,Pop == "PAMUR10" | Pop == "PHAYLY10" | Pop == "PKOPE96" | Pop == "PKUSHI06" | Pop == "PNOME94" | Pop == "PSNOH96" | Pop == "PTAUY12")
OLD_16662_even[1:5,1:5]
dim(OLD_16662_even)

OLD_16662_odd <-subset(OLD16662_pops, Pop == "PAMUR11" | Pop == "PHAYLY09" | Pop == "PKOPE91" | Pop == "PKUSHI07" | Pop == "PNOME91" | Pop == "PSNOH03" | Pop == "PTAUY09")
OLD_16662_odd[1:5,1:5]
dim(OLD_16662_odd)

NEW_16662_Lakel_even_genos <- NEW_16662_Lakel_even[,-1]
NEW_16662_Lakel_odd_genos <- NEW_16662_Lakel_odd[,-1]
OLD_16662_even_genos <- OLD_16662_even[,-1]
OLD_16662_odd_genos <- OLD_16662_odd[,-1]


### Calculate allele frequencies with Adegenet

NEW_16662_genepop <- read.genepop("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/Lakel_lineage_test/batch_4_NEW_Even_Odd.gen")
OLD_16662_genepop <- read.genepop("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/Lakel_lineage_test/batch_4_OLD_Even_Odd.gen")


NEW_freq <- makefreq(NEW_16662_genepop, missing = NA)
OLD_freq <- makefreq(OLD_16662_genepop, missing = NA)


# Stop Here 
# 
# ###################################
# 
# ###Script to calculate, output and graph allele frequencies
# ###Created by Mark Christie on 3/1/2012
# 
# #gdata<-read.table(file="BeechData_east.txt", header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)  #load data ; here genotypes of Beech data from dryad
# #head(gdata) #look at the data
# 
# #population=gdata[,-c(1:3)]  #Here we are accessing only the genotypes (i.e., removing columns 1-2)
# 
# NEW_16662_Lakel_even_genos <- NEW_16662_Lakel_even[,-1]
# NEW_16662_Lakel_odd_genos <- NEW_16662_Lakel_odd[,-1]
# OLD_16662_even_genos <- OLD_16662_even[,-1]
# OLD_16662_odd_genos <- OLD_16662_odd[,-1]
# 
# NEW_16662_Lakel_even_genos[1:6,1:6]
# dim(NEW_16662_Lakel_even_genos)
# 
# NEW_16662_Lakel_odd_genos[1:6,1:6]
# dim(NEW_16662_Lakel_odd_genos)
# 
# OLD_16662_even_genos[1:6,1:6]
# dim(OLD_16662_even_genos)
# 
# OLD_16662_odd_genos[1:6,1:6]
# dim(OLD_16662_odd_genos)
# 
# #L=ncol(population)   #how many columns are there?
# a <- ncol(NEW_16662_Lakel_even_genos)
# b <- ncol(NEW_16662_Lakel_odd_genos)
# c <- ncol(OLD_16662_even_genos)
# d <- ncol(OLD_16662_odd_genos)
# 
# 
# #locus_positions=(2*(unique(round((1:(L-2))/2)))+1)   #find the starting column number for each locus
# 
# a_locus_positions <- (2*(unique(round((1:(a-2))/2)))+1)
# b_locus_positions <- (2*(unique(round((1:(b-2))/2)))+1)
# c_locus_positions <- (2*(unique(round((1:(c-2))/2)))+1)
# d_locus_positions <- (2*(unique(round((1:(d-2))/2)))+1)
# 
# #lnames=colnames(population)                          #locus names, from the header
# anames <- colnames(NEW_16662_Lakel_even_genos)                          
# bnames <- colnames(NEW_16662_Lakel_odd_genos) 
# cnames <- colnames(OLD_16662_even_genos) 
# dnames <- colnames(OLD_16662_odd_genos)
# 
# # OUT=NULL        #create a null dataset to append allele freqs to
# # 
# # for (x in locus_positions) {                       #begin for loop, to calculate frequencies for each locus
# #   alleles=c(population[,x],population[,x+1])        #For example, combine columns 1 and 2 for locus 1 (two columns because they are diploid)
# #   alleles2=as.data.frame(table(alleles))             #count each allele at locus x
# #   missing=alleles2[which(alleles2[,1]==0),2]          #count missing data at locus x, entered as '0' in this dataset (not used further for simplicity)
# #   alleles3=alleles2[-which(alleles2[,1]==0),]          #remove missing data (otherwise 0 would be counted in total number of alleles)
# #   alleles4=cbind(alleles3,alleles3[,2]/sum(alleles3[,2])) #calculate frequencies
# #   output=cbind(x,lnames[x],alleles4)                        #combine x, locusname, and frequencies
# #   
# #   
# #   OUT <<- rbind(OUT,output)
# #   
# # }
# # 
# # 
# # colnames(OUT) <- c("Number","Locus","allele","count","frequency") #add column headers
# # Allelefreqs <- OUT[,-1]
# # write.table(Allelefreqs,file="Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/Lakel_lineage_test/Allelefrequencies.txt",row.names=FALSE,col.names=TRUE,sep="\t",append=FALSE)

#________________________________________________________________________
locipassedMAF_temp<-locipassedMAF_temp[,1]
head(locipassedMAF_temp)
length(locipassedMAF_temp)
#locipassedMAF_temp<-locipassedMAF #to re-set
intersect_filteredWhitelist_16681 <- intersect(locipassedMAF_temp, pink16681)
diff_filteredWhitelist_16681 <-setdiff(pink16681,locipassedMAF_temp)
diff_16681_filteredWhitelist <-setdiff(locipassedMAF_temp, pink16681)

length(intersect_filteredWhitelist_16681)
head(locipassedMAF_temp)

length(diff_filteredWhitelist_16681)
head(diff_filteredWhitelist_16681)

length(diff_16681_filteredWhitelist)
head(diff_16681_filteredWhitelist)

#########################################