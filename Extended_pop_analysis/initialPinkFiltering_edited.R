### Initial Filtering of Pink Stacks Population Output
###    combines whitelists, and filters on MAF and locus genotype rate to create a final whitelist
###    Uses a per population MAF instead of global filter and a per locus genotype rate
### Garrett McKinney and Carolyn Tarpey | October 2017
### ---------------------------------------

#This is R code originally written by Garrett McKinney to combine and filter the output of Stacks populations. 
#It requires a genepop file for each of the whitelists used to run the populations command in Stacks. 
#The genepop files should be stripped of the header line that Stacks puts in and the second line should start with a tab so that 
#R will load it in as a table and the loci names should be tab deliminated, and there should be no "pop" between the samples.

#install.packages("vcfR")
#install.packages("stringr")
#install.packages("ggplot2")
#install.packages("lazyeval")
#install.packages("adegenet")
#install.packages("dplyr")
#install.packages("DBI")

#Pink data filtering
library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)

#load genepop files
genos1<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/OriginalWhitelists/Whitelist_Genepops/whitelist_1.genepop",colClasses="factor")
genos1[1:5,1:5]
dim(genos1)
genos2<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/OriginalWhitelists/Whitelist_Genepops/whitelist_2.genepop",colClasses="factor")
genos2[1:5,1:5]
dim(genos2)
genos3<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/OriginalWhitelists/Whitelist_Genepops/whitelist_3.genepop",colClasses="factor")
genos3[1:5,1:5]
dim(genos3)
genos4<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/OriginalWhitelists/Whitelist_Genepops/whitelist_4.genepop",colClasses="factor")
genos4[1:5,1:5]
dim(genos4)
genos5<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/OriginalWhitelists/Whitelist_Genepops/whitelist_5.genepop",colClasses="factor")
genos5[1:5,1:5]
dim(genos5)
genos6<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/OriginalWhitelists/Whitelist_Genepops/whitelist_6.genepop",colClasses="factor")
genos6[1:5,1:5]
dim(genos6)
genos7<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/OriginalWhitelists/Whitelist_Genepops/whitelist_7.genepop",colClasses="factor")
genos7[1:5,1:5]
dim(genos7)
genos8<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/OriginalWhitelists/Whitelist_Genepops/whitelist_8.genepop",colClasses="factor")
genos8[1:5,1:5]
dim(genos8)
genos9<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/OriginalWhitelists/Whitelist_Genepops/whitelist_9.genepop",colClasses="factor")
genos9[1:5,1:5]
dim(genos9)
genos10<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/OriginalWhitelists/Whitelist_Genepops/whitelist_10.genepop",colClasses="factor")
genos10[1:5,1:5]
dim(genos10)
genos11<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/OriginalWhitelists/Whitelist_Genepops/whitelist_11.genepop",colClasses="factor")
genos11[1:5,1:5]
dim(genos11)
genos12<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/OriginalWhitelists/Whitelist_Genepops/whitelist_12.genepop",colClasses="factor")
genos12[1:5,1:5]

#use cbind to combine data into single dataset
allGenos<-cbind(genos1,genos2,genos3,genos4,genos5,genos6,genos7,genos8,genos9,genos10,genos11,genos12)
allGenos[1:5,1:5]
dim(allGenos)

#get unique tags from full dataset
allLoci<-data.frame(str_split_fixed(colnames(allGenos),"_",2))
colnames(allLoci)<-c("Tag","SNP")
allLoci$Locus<-colnames(allGenos)
length(allLoci$Locus)
length(unique(allLoci$Tag))
#write.table(allLoci,"Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/RawLoci_unfiltered.txt",quote=FALSE,row.names=FALSE)
#write.table(allGenos,"Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/RawGenos_unfiltered.txt",quote=FALSE,row.names=FALSE)

#get genotype rate per locus
locusGenoRate<-apply(allGenos,2,function(x) 1-(sum(x=="0000")/dim(allGenos)[1]))
locusGenoRate<-data.frame(keyName=names(locusGenoRate), value=locusGenoRate, row.names=NULL)
colnames(locusGenoRate)<-c("Locus","GenoRate")
dim(locusGenoRate)
head(locusGenoRate)
#write.table(locusGenoRate,"Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/LocusGenoRate.txt",quote=FALSE,row.names=FALSE)

#filterLoci with >= 80% genotype rate
filteredLoci<-locusGenoRate[locusGenoRate$GenoRate>=0.80,]
dim(filteredLoci)
head(filteredLoci)

###########################################
#if we want to see which of our 16681 snps is included at this point
loci_16681 <- readLines("Z:/WORK/TARPEY/Pink_Populations/listof16681LOCI.txt")
head(loci_16681)

#list of all the loci names in allGenos, no "X" in the name 
allGenos_loci_test <- gsub("X","",colnames(allGenos))

#loci that are in 16681 and 110433 
in_16681_and_110433 <- intersect (loci_16681, allGenos_loci_test)
length (in_16681_and_110433)
setdiff(loci_16681, allGenos_loci_test)


#if we want to see what was eliminated here
excludedLoci_rate<-locusGenoRate[locusGenoRate$GenoRate<=0.80,]
dim(excludedLoci_rate)
head(excludedLoci_rate)

excludedLoci<- as.vector(excludedLoci_rate[,1])
head(excludedLoci)
class(excludedLoci)

#if we want to see which of our 16681 snps is included at this point
#loci that are in 16681 and 76761 (80% locus genotype rate)
#list of all the loci names in allGenos, no "X" in the name 
filteredLoci_test <- gsub("X","", filteredLoci$Locus)

in_16681_and_76761 <- intersect (loci_16681, filteredLoci_test)
length(in_16681_and_110433)
setdiff(loci_16681, filteredLoci_test)

#We want to know the number of loci at the 31485 unique tags
allLoci<-data.frame(str_split_fixed(colnames(allGenos),"_",2))
colnames(allLoci)<-c("Tag","SNP")
allLoci$Locus<-colnames(allGenos)
length(allLoci$Locus)
head(allLoci)


#length(unique(allLoci$Tag))

############################################


#get unique tags from 80% genoRate filtered dataset
retainedLoci_80perc<-data.frame(str_split_fixed(filteredLoci$Locus,"_",2))
colnames(retainedLoci_80perc)<-c("Tag","SNP")
retainedLoci_80perc$Locus<-filteredLoci$Locus
head(retainedLoci_80perc)
length(unique(retainedLoci_80perc$Tag))

#format dataset to remove X from locus names
filteredLociIDs<-filteredLoci$Locus
filteredLociIDs<-gsub("X","",filteredLociIDs)

#filter dataset to include only filtered loci
filteredGenos<-allGenos[,colnames(allGenos)%in%filteredLoci$Locus]
dim(filteredGenos)

allGenos[1:5,1:5]
filteredGenos[1:5,1:5]
dim(filteredGenos)
write.table(filteredGenos,"Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/FirstFiltering/filteredGenos_just80PCgenorate_76761.txt",quote=FALSE,row.names=TRUE)



#get minor allele frequency for loci
calculateMAF<-function(genotypes){
  genotypeList<-sort(unique(genotypes))
  genotypeList<-genotypeList[genotypeList != "0000"]
  allelesList1<-substr(genotypeList,1,2)
  allelesList2<-substr(genotypeList,3,4)
  allelesList<-unique(c(allelesList1,allelesList2))
  allele1Counts<-sum(str_count(genotypes,allelesList[1]))
  allele2Counts<-sum(str_count(genotypes,allelesList[2]))
  if(length(allelesList)==1){
    MAF=0
  }else if(allele1Counts>=allele2Counts){
    MAF<-allele2Counts/(allele1Counts+allele2Counts)
  }else{
    MAF<-allele1Counts/(allele1Counts+allele2Counts)
  }
  return(MAF)
}

#assign each sample to a population, delete the first column (duplicate sample name) and rename the new column Pop
popmap<- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/PopMap.txt")
filteredGenos<-cbind(popmap,filteredGenos)
filteredGenos$V1<- NULL
colnames(filteredGenos)[1]<-"Pop"
filteredGenos[1:5,1:5]
#create a dataframe for our MAF results
npops<- length(unique(filteredGenos$Pop))
nloci<- length(filteredGenos[,-1])

filteredGenos_POP_MAF<-matrix(nrow=npops, ncol=nloci) 
rownames(filteredGenos_POP_MAF)<-as.vector(unique(filteredGenos$Pop)) #name the rows by the population names
colnames(filteredGenos_POP_MAF)<-colnames(filteredGenos[,-1])

for (i in 1:length(unique(filteredGenos$Pop))){
  tempPop<-(unique(filteredGenos$Pop)[i])
  tempGeno<-subset(filteredGenos, Pop == tempPop)
  #print(tempPop)
  #print(tempGeno[1:5,1:5])
  tempMAF<-apply(tempGeno[,-1],2,calculateMAF)
  #print(head(tempMAF))
  filteredGenos_POP_MAF[i,]<-c(tempMAF)
  #print(filteredGenos_POP_MAF[1:i,1:7])
}

filteredGenos_POP_MAF[1:5,1:5]
dim(filteredGenos_POP_MAF)
#colnames(filteredGenos_POP_MAF[,1])

##### Filter the loci by the MAF results by population
##retain loci that were at least 0.05 in any of the 18 populations
lociMAF<-vector()
lociMAF<-apply(filteredGenos_POP_MAF,2,function(x) sum(x >=0.05, na.rm=TRUE))
lociMAF<-data.frame(keynames=names(lociMAF), value=lociMAF, row.names = NULL)
colnames(lociMAF)<-c("Locus", "PopsMAF")

#if the loci does not pass the test, delete it
locipassedMAF<-vector()
locipassedMAF<-subset(lociMAF, lociMAF[,2]!=0)
dim(locipassedMAF)
head(locipassedMAF)

# a second copy to use later to see what matches with our 16681
locipassedMAF_temp<-vector()
locipassedMAF_temp<-subset(lociMAF, lociMAF[,2]!=0)
head(locipassedMAF_temp)
dim(locipassedMAF_temp)
#locipassedMAF <-locipassedMAF_temp #reset

#get unique tags from initial filtered dataset
locipassedMAF<-data.frame(str_split_fixed(locipassedMAF$Locus,"_",2))
colnames(locipassedMAF)<-c("Tag","SNP")

## to figure out how many loci from the original data set will be genotyped at these tags
UNIQUEtags31485<-unique(locipassedMAF$Tag)
length(UNIQUEtags31485)

#locipassedMAF_tags$Locus<-locipassedMAF$Locus
head(locipassedMAF)
# to get the unique tags for one SNP per tag
locipassedMAF<-unique(locipassedMAF$Tag)
length(locipassedMAF)






###########################################
#####################Test to see that what we got in our whitelist has the 16681 loci that we want

#if we want to see which of our 16681 snps is included at this point
loci_16681 <- readLines("Z:/WORK/TARPEY/Pink_Populations/listof16681LOCI.txt")
head(loci_16681)

#list of all the loci names in locipassedMAF_temp, no "X" in the name 
#locipassedMAF_temp$Locus <- gsub("X","",locipassedMAF_temp$Locus)
MAFlocilist <- gsub("X","",locipassedMAF_temp$Locus)
length(MAFlocilist)

#loci that are in 16662 and 48687 
length(in_16681_and_110433)
in_16662_and_48687 <- intersect(in_16681_and_110433, MAFlocilist)
length(in_16662_and_48687)

setdiff(in_16681_and_110433, MAFlocilist)


### to figure out which loci will be included in the STACKS output files for the unique 31485 tags
head(UNIQUEtags31485)
allLoci[1:2,1:2]
#UNIQUEtags31485 <- as.list(UNIQUEtags31485)
loci4Tags <-allLoci[allLoci$Tag%in%UNIQUEtags31485,1:2]
#filteredGenos<-allGenos[,colnames(allGenos)%in%filteredLoci$Locus]
head(loci4Tags)
length(loci4Tags$Tag)

##Compare the 63758 loci that are in this set that is pulled from the 31485 loci with the 16681
loci4Tags$Tag<- gsub("X","",loci4Tags$Tag) #remove the X in the tag name
head(loci4Tags)
loci4Tags$Locus<- paste(loci4Tags$Tag,loci4Tags$SNP, sep="_") #this concatenates the two columns with the _ in the middle
LociIN31485tags<- loci4Tags$Locus
head(LociIN31485tags)
in_16662_and_63758 <- intersect(in_16681_and_110433, LociIN31485tags)
length(in_16662_and_63758)

###########################################
#########################################

#output list of filtered loci, need to open as binary file (using write binary, "wb") to give unix line endings
#format dataset to remove X from locus names
firstfilteredLociIDs<-gsub("X","",locipassedMAF)
head(firstfilteredLociIDs)
outputFile<-file("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/firstfilteredLociIDs", "wb")
write.table(firstfilteredLociIDs,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)


