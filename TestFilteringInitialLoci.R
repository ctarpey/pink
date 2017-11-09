### Identifies the # of loci or samples lost at each filtering step of the 
### Initial Filtering of Pink Stacks Population Output
### 
###  Garrett McKinney and Carolyn Tarpey | October 2017
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
genos1<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_1.genepop",colClasses="factor")
genos1[1:5,1:5]
dim(genos1)
genos2<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_2.genepop",colClasses="factor")
genos2[1:5,1:5]
dim(genos2)
genos3<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_3.genepop",colClasses="factor")
genos3[1:5,1:5]
dim(genos3)
genos4<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_4.genepop",colClasses="factor")
genos4[1:5,1:5]
dim(genos4)
genos5<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_5.genepop",colClasses="factor")
genos5[1:5,1:5]
dim(genos5)
genos6<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_6.genepop",colClasses="factor")
genos6[1:5,1:5]
dim(genos6)
genos7<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_7.genepop",colClasses="factor")
genos7[1:5,1:5]
dim(genos7)
genos8<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_8.genepop",colClasses="factor")
genos8[1:5,1:5]
dim(genos8)
genos9<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_9.genepop",colClasses="factor")
genos9[1:5,1:5]
dim(genos9)
genos10<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_10.genepop",colClasses="factor")
genos10[1:5,1:5]
dim(genos10)
genos11<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_11.genepop",colClasses="factor")
genos11[1:5,1:5]
dim(genos11)
genos12<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_12.genepop",colClasses="factor")
genos12[1:5,1:5]


#use cbind to combine data into single dataset
allGenos<-cbind(genos1,genos2,genos3,genos4,genos5,genos6,genos7,genos8,genos9,genos10,genos11,genos12)
allGenos[1:5,1:5]
dim(allGenos)

#get unique tags from full dataset
allLoci<-data.frame(str_split_fixed(colnames(allGenos),"_",2))
colnames(allLoci)<-c("Tag","SNP")
allLoci$Locus<-colnames(allGenos)
head(allLoci)
length(allLoci$Locus)
length(unique(allLoci$Tag))

#get genotype rate per locus
locusGenoRate<-apply(allGenos,2,function(x) 1-(sum(x=="0000")/dim(allGenos)[1]))
locusGenoRate<-data.frame(keyName=names(locusGenoRate), value=locusGenoRate, row.names=NULL)
colnames(locusGenoRate)<-c("Locus","GenoRate")
dim(locusGenoRate)

#get genotype rate per sample
sampleGenoRate<-apply(allGenos,1,function(x) 1-(sum(x=="0000")/dim(allGenos)[2]))
sampleGenoRate<-data.frame(keyName=names(sampleGenoRate), value=sampleGenoRate, row.names=NULL)
colnames(sampleGenoRate)<-c("Sample","GenoRate")
dim(sampleGenoRate)
head(sampleGenoRate)

#plot ranked genotype rate for loci
locusGenoRate_ranked<-locusGenoRate[order(locusGenoRate$GenoRate),]
locusGenoRate_ranked$rank<-seq(1,dim(locusGenoRate_ranked)[1],by=1)
head(locusGenoRate_ranked)

##################### Testing where the loci that we expected to be retained were filtered out
#input the list of the 7,615 loci that were filtered out through this process that we want to retain.

missing7615 <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/missing7615.txt")
missing7615 <- unlist(missing7615, use.names = FALSE)
length(missing7615)

#filterLoci with >= 70% genotype rate
filteredLoci<-locusGenoRate[locusGenoRate$GenoRate>=0.70,]
dim(filteredLoci)

#flip that to keep a list of the loci that were filtered out here
excludedLoci_rate<-locusGenoRate[locusGenoRate$GenoRate<=0.70,]
dim(excludedLoci_rate)
head(excludedLoci_rate)

excludedLoci<- as.vector(excludedLoci_rate[,1])
head(excludedLoci)
class(excludedLoci)

#Compare the two vectors to see what there is in common
Ex_locusGenRate70 <- intersect(excludedLoci, missing7615)
length(Ex_locusGenRate70)

#Ex_locusGenRate70 <- union(excludedLoci, missing7615) #test to make sure its working like we think


#get unique tags from 70% genoRate filtered dataset
retainedLoci_70perc<-data.frame(str_split_fixed(filteredLoci$Locus,"_",2))
colnames(retainedLoci_70perc)<-c("Tag","SNP")
retainedLoci_70perc$Locus<-filteredLoci$Locus
head(retainedLoci_70perc)
length(unique(retainedLoci_70perc$Tag))

#format dataset to remove X from locus names
filteredLociIDs<-filteredLoci$Locus
filteredLociIDs<-gsub("X","",filteredLociIDs)

#filter dataset to include only filtered loci
filteredGenos<-allGenos[,colnames(allGenos)%in%filteredLoci$Locus]
dim(filteredGenos)

#get genotype rate per sample for filtered loci
sampleGenoRate_filteredLoci<-apply(filteredGenos,1,function(x) 1-(sum(x=="0000")/dim(filteredGenos)[2]))
sampleGenoRate_filteredLoci<-data.frame(keyName=names(sampleGenoRate_filteredLoci), value=sampleGenoRate_filteredLoci, row.names=NULL)
colnames(sampleGenoRate_filteredLoci)<-c("Sample","GenoRate")
dim(sampleGenoRate_filteredLoci)
head(sampleGenoRate_filteredLoci)


#plot ranked genotype rate for samples with filtered loci
sampleGenoRate_filteredLoci_ranked<-sampleGenoRate_filteredLoci[order(sampleGenoRate_filteredLoci$GenoRate),]
sampleGenoRate_filteredLoci_ranked$rank<-seq(1,dim(sampleGenoRate_filteredLoci_ranked)[1],by=1)

#filter samples with >=70% genotype rate
filteredSamples<-sampleGenoRate_filteredLoci[sampleGenoRate_filteredLoci$GenoRate>=0.70,]
filteredGenos_filteredSamples<-filteredGenos[rownames(filteredGenos)%in%filteredSamples$Sample,]
dim(filteredGenos_filteredSamples)

##### See how many loci we lose at the filter for genotype rate per sample

#flip that to keep a list of the samples that were filtered out here
excludedSamples<-sampleGenoRate_filteredLoci[sampleGenoRate_filteredLoci$GenoRate<=0.70,]
excludedSamples

#get genotype rate per locus for filtered sample dataset
locusGenoRate_filteredSamples<-apply(filteredGenos_filteredSamples,2,function(x) 1-(sum(x=="0000")/dim(allGenos)[1]))
locusGenoRate_filteredSamples<-data.frame(keyName=names(locusGenoRate_filteredSamples), value=locusGenoRate_filteredSamples, row.names=NULL)
colnames(locusGenoRate_filteredSamples)<-c("Locus","GenoRate")
dim(locusGenoRate_filteredSamples)
head(locusGenoRate_filteredSamples)

#plot ranked genotype rate for loci
locusGenoRate_filteredSamples_ranked<-locusGenoRate_filteredSamples[order(locusGenoRate_filteredSamples$GenoRate),]
locusGenoRate_filteredSamples_ranked$rank<-seq(1,dim(locusGenoRate_filteredSamples_ranked)[1],by=1)
head(locusGenoRate_filteredSamples_ranked)

##### See how many loci we lose at the filter for Global MAF


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

filteredGenos_filteredSamples_MAF<-apply(filteredGenos_filteredSamples,2,calculateMAF)
filteredGenos_filteredSamples_MAF<-data.frame(keyName=names(filteredGenos_filteredSamples_MAF), value=filteredGenos_filteredSamples_MAF, row.names=NULL)
colnames(filteredGenos_filteredSamples_MAF)<-c("Locus","MAF")
head(filteredGenos_filteredSamples_MAF)

#combine MAF and genotype rate summaries
filteredLoci_summary<-merge(locusGenoRate_filteredSamples,filteredGenos_filteredSamples_MAF,by="Locus")
head(filteredLoci_summary)

#filter loci to retain only those with averall MAF>=0.05
MAF05_loci<-filteredGenos_filteredSamples_MAF[filteredGenos_filteredSamples_MAF$MAF>=0.05,]
filteredGenos_filteredSamples_MAF05<-filteredGenos_filteredSamples[,colnames(filteredGenos_filteredSamples)%in%MAF05_loci$Locus]
dim(filteredGenos_filteredSamples_MAF05)

###See how many loci we lose in the Global MAF filter
ex_MAF05_loci<-filteredGenos_filteredSamples_MAF[filteredGenos_filteredSamples_MAF$MAF<=0.05,]
head(ex_MAF05_loci)

#head(excludedSamples_filteredSamples)
ex_MAF05_locinames<- as.vector(ex_MAF05_loci[,1])
head(ex_MAF05_locinames)
class(ex_MAF05_locinames)

intersect_exMAF05_Missing <- intersect(ex_MAF05_locinames, missing7615)
length(intersect_exMAF05_Missing)


