
### Secondary Filtering of Pink Stacks Population Output
###   This code takes the final genepop files from the Stacks file and filters the genotypes 
###    
###  Carolyn Tarpey | November 2017
### ---------------------------------------

#This is adapted from R code written by Garrett McKinney to combine and filter the output of Stacks populations. 
#It requires a genepop file stripped of the header line, the second line should start with a tab so that 
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
all_newGenos<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/batch_4-genepop_4R.txt",colClasses="factor")
all_newGenos[1:5,1:5]
dim(all_newGenos)


#get minor allele frequency for loci- want to keep SNP with highest MAF at each tag

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


allGenos_oneSNP<-apply(all_newGenos,2,calculateMAF)
head(allGenos_oneSNP)


#put the results in a dataframe
allGenos_oneSNP_temp<-data.frame(value=allGenos_oneSNP,row.names=names(all_newGenos))
colnames(allGenos_oneSNP_temp)<-c("MAF")
head(allGenos_oneSNP_temp)

#concatenate the locus names to the results. 
allLoci<-colnames(all_newGenos)
length(allLoci)
allGenos_oneSNP_temp$Locus<-allLoci
head(allGenos_oneSNP_temp)

#split the tags and snp positions
allGenos_oneSNP_temp_tags<-data.frame(str_split_fixed(allGenos_oneSNP_temp$Locus,"_",2))
colnames(allGenos_oneSNP_temp_tags)<-c("Tag","SNP")
head(allGenos_oneSNP_temp_tags)
allGenos_oneSNP_temp$Tag<-allGenos_oneSNP_temp_tags$Tag
allGenos_oneSNP_temp$SNP<-allGenos_oneSNP_temp_tags$SNP
head(allGenos_oneSNP_temp)
#write.table(allGenos_oneSNP_temp,"Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/allGenos_oneSNP_temp.txt",quote=FALSE,row.names=FALSE)

#Retain the SNP with the highest MAF per tag
# This one gives an output with the correct number of unique tags, and spot checked in excel
oneMAF<- allGenos_oneSNP_temp %>% group_by(Tag) %>% slice(which.max(MAF))
head(oneMAF)
#write.table(oneMAF,"Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/oneMAF.txt",quote=FALSE,row.names=FALSE)
oneMAF<-as.data.frame(oneMAF)
head(oneMAF)


#####################Test to see that what we got in One SNP per tag by MAF has the 16681 loci that we want

pink16681 <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/pink16681.txt")
pink16681 <- unlist(pink16681, use.names = FALSE)
length(pink16681)

oneMAF_temp<-oneMAF$Locus
head(oneMAF_temp)
#locipassedMAF_temp<-locipassedMAF #to re-set
intersect_oneMAF_16681 <- intersect(oneMAF_temp, pink16681)
diff_oneMAF_16681 <-setdiff(pink16681,oneMAF_temp)

length(intersect_oneMAF_16681)
head(intersect_oneMAF_16681)

length(diff_oneMAF_16681)
head(diff_oneMAF_16681)

###########################################

#filter dataset to include only one SNP per tag loci
length(oneMAF$Locus)

filteredGenos_oneSNP<-all_newGenos[,oneMAF$Locus]
dim(filteredGenos_oneSNP)

all_newGenos[1:5,1:5]
filteredGenos_oneSNP[1:5,1:5]


############# Write a txt file with the list of loci- use as whitelist or subset genepop file

#subset genepop file (for use with subset_genepop_by_SNPs.pl):

#format dataset to remove X from locus names
oneMAF_loci<-gsub("X","",oneMAF$Locus)
head(oneMAF_loci)

outputFile<-file("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/oneSNP_subsetgenepop.txt", "wb")
write.table(oneMAF_loci,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

#use as whitelist in Stacks: 

oneMAF_loci_whitelist<-data.frame(str_split_fixed(oneMAF_loci,"_",2))

outputFile<-file("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/oneSNP_whitelist.txt", "wb")
write.table(oneMAF_loci_whitelist[,1],outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

##################### Run Genepop with the genepop file for HWE and Fis and import the results here for filtering based on those 














write.table(allLoci,"Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/RawLoci_unfiltered.txt",quote=FALSE,row.names=FALSE)
#write.table(allGenos,"Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/RawGenos_unfiltered.txt",quote=FALSE,row.names=FALSE)

#setdiff(x,y) : Set difference between x and y , consisting of all elements of x that are not in y.

#output list of filtered loci, need to open as binary file (using write binary, "wb") to give unix line endings
#format dataset to remove X from locus names
filteredLociIDs<-filteredLoci$Locus
filteredLociIDs<-gsub("X","",filteredLociIDs)
outputFile<-file("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/70per_genoRateLoci.txt", "wb")
write.table(filteredLociIDs,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

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
#ggplot()+geom_point(data=sampleGenoRate_filteredLoci_ranked,aes(x=rank,y=GenoRate))+geom_hline(aes(yintercept=0.70),lty="dashed")+ggtitle("Sample Genotype Rate of Loci Genotyped in 70% of Samples")
#ggsave(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/SampleGenotypeRateofLociGenotypedin70Samples.pdf")
#ggplot()+geom_point(data=sampleGenoRate_filteredLoci_ranked,aes(x=rank,y=GenoRate))+xlab("Sample Genotype Rate Rank")+ylab("Genotype Rate")+theme_bw()+ggtitle("Sample Genotype Rate of Loci Genotyped in 70% of Samples")
##plot histogram of genotype rate filter for samples with filtered loci
#ggplot()+geom_histogram(data=sampleGenoRate_filteredLoci,aes(x=GenoRate),binwidth=0.01)+ggtitle("Sample Genotype Rate of Loci Genotyped in 70% of Samples")

#filter samples with >=80% genotype rate
filteredSamples<-sampleGenoRate_filteredLoci[sampleGenoRate_filteredLoci$GenoRate>=0.70,]
filteredGenos_filteredSamples<-filteredGenos[rownames(filteredGenos)%in%filteredSamples$Sample,]
dim(filteredGenos_filteredSamples)

#get genotype rate per locus for filtered sample dataset
locusGenoRate_filteredSamples<-apply(filteredGenos_filteredSamples,2,function(x) 1-(sum(x=="0000")/dim(allGenos)[1]))
locusGenoRate_filteredSamples<-data.frame(keyName=names(locusGenoRate_filteredSamples), value=locusGenoRate_filteredSamples, row.names=NULL)
colnames(locusGenoRate_filteredSamples)<-c("Locus","GenoRate")
dim(locusGenoRate_filteredSamples)
head(locusGenoRate_filteredSamples)



#_____________________________________________________________________________________________________________________
#plot histogram of MAF 0.05 filtered loci
ggplot()+geom_histogram(data=filteredGenos_filteredSamples_MAF[filteredGenos_filteredSamples_MAF$MAF>=0.05,],aes(x=MAF),binwidth=0.01)+ggtitle("Histogram of 0.05 MAF Filtered Loci")
#plot genotype rate per locus relative to MAF for MAF 0.02 filtered loci
ggplot()+geom_point(data=filteredLoci_summary[filteredLoci_summary$MAF>=0.05,],aes(x=GenoRate,y=MAF),alpha=0.1)

#filter loci to retain only those with averall MAF>=0.05
MAF05_loci<-filteredGenos_filteredSamples_MAF[filteredGenos_filteredSamples_MAF$MAF>=0.05,]
filteredGenos_filteredSamples_MAF05<-filteredGenos_filteredSamples[,colnames(filteredGenos_filteredSamples)%in%MAF05_loci$Locus]
dim(filteredGenos_filteredSamples_MAF05)

#get unique tags from filtered dataset
retainedLoci<-data.frame(str_split_fixed(colnames(filteredGenos_filteredSamples_MAF05),"_",2))
colnames(retainedLoci)<-c("Tag","SNP")
retainedLoci$Locus<-colnames(filteredGenos_filteredSamples_MAF05)
head(retainedLoci)

#plot histogram of SNP position
#histogram won't work because it's discrete
#ggplot()+geom_histogram(data=retainedLoci,aes(x=SNP),binwidth=1)
#order SNP position for plotting
#sort((unique(retainedLoci$SNP)))
#retainedLoci$SNP <- factor(retainedLoci$SNP, levels = c(...))
#ggplot()+geom_bar(data=retainedLoci,aes(x=SNP))
#get number of retained tags and SNPs

retainedTagCount<-length(unique(retainedLoci$Tag))
retainedTagCount
retainedSNPcount<-length(retainedLoci$Locus)
#get number of SNPs per locus
SNPdistribution<-table(retainedLoci$Tag)
SNPdistribution<-as.data.frame(SNPdistribution)
colnames(SNPdistribution)<-c("Tag","SNPnum")
head(SNPdistribution)
#plot histogram of SNP number
ggplot()+geom_histogram(data=SNPdistribution,aes(x=SNPnum),binwidth=1)


#assign individuals to populations
samples<-rownames(allGenos)
samplePops<-data.frame(str_split_fixed(samples,"_",2))
colnames(samplePops)<-c("pop","sample")
samplePops$sample<-rownames(allGenos)
head(samplePops)


# #use cbind bind to the population map to the dataset
# popmap<- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/PopMap.txt")
# colnames(popmap)<-c("Sample","POP")
# allGenos<-cbind(popmap,allGenos)
# allGenos$Sample <- NULL #this deletes the first column which was a repeat of the sample name. 
# allGenos[1:5,1:5]
# dim(allGenos)

