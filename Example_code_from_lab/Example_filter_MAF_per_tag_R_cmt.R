###  Chinook panel: filtering raw RAD genotypes 
###    for new RAD loci to include in a refined panel
###    
### Carolyn Tarpey | June 2018 
### ---------------------------------------


library(vcfR)
library(stringr)
library(dplyr)
library(gdata)
library(reshape2)
library(plotly)


#load the raw genepop file, edited for R. Includes all individuals and loci:  
chin_raw_genepop <-read.delim("Z:/WORK/TARPEY/ChinookPanel/RAD_data/NoCookInlet/batch_13.allCombined_NoCookInlet_R_genepop.txt", sep="", header = TRUE, colClasses="factor")
dim(chin_raw_genepop)


######################################
#remove the commas from the sample names, and the X from the column names
chin_raw_genepop[1:15,1:15]
rownames(chin_raw_genepop)<-gsub(",","",rownames(chin_raw_genepop))
colnames(chin_raw_genepop)<-gsub("X","",colnames(chin_raw_genepop))
chin_raw_genepop[1:15,1:15]

####################### Filtering for highest MAF SNP per tag is usually one of the last things I do, so by that point 
####I have another name for the filtered genepop file. To help the rest of the code make sense, just rename  your 
### raw genepop file filtered_MAF_Genos

filtered_MAF_Genos <- chin_raw_genepop

################## Create a dataset that is one SNP per tag- Highest MAF ############
# use the genotypes from filtered_MAF_Genos

#the function to calculate MAF
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

filtered_MAF_Genos[1:5,1:15]

#this calculates the MAF overall. ###If I have added meta data to the genepop file, I have to change the 
## Column numbers here: c(1:(dim(filtered_MAF_Genos)[2])) so that it skips columns that arent genotypes
filtered_MAF_Genos_oneSNP<-apply(filtered_MAF_Genos[,c(1:(dim(filtered_MAF_Genos)[2]))],2,calculateMAF)
head(filtered_MAF_Genos_oneSNP)

#put the results in a dataframe
filtered_MAF_Genos_oneSNP_temp<-data.frame(value=filtered_MAF_Genos_oneSNP,row.names=names(filtered_MAF_Genos_oneSNP))
colnames(filtered_MAF_Genos_oneSNP_temp)<-c("MAF")
head(filtered_MAF_Genos_oneSNP_temp)

#concatenate the locus names to the results. ###If I have added meta data to the genepop file, I have to change the 
## Column numbers here: c(1:(dim(filtered_MAF_Genos)[2])) so that it doesnt rename columns that arent genotypes
Loci_temp<-colnames(filtered_MAF_Genos[,c(1:(dim(filtered_MAF_Genos)[2]))])
length(Loci_temp)
filtered_MAF_Genos_oneSNP_temp$Locus<-Loci_temp
head(filtered_MAF_Genos_oneSNP_temp)

#split the tags and snp positions
filtered_MAF_Genos_oneSNP_temptags<-data.frame(str_split_fixed(filtered_MAF_Genos_oneSNP_temp$Locus,"_",2))
colnames(filtered_MAF_Genos_oneSNP_temptags)<-c("Tag","SNP")
head(filtered_MAF_Genos_oneSNP_temptags)
filtered_MAF_Genos_oneSNP_temp$Tag<-filtered_MAF_Genos_oneSNP_temptags$Tag
filtered_MAF_Genos_oneSNP_temp$SNP<-filtered_MAF_Genos_oneSNP_temptags$SNP
head(filtered_MAF_Genos_oneSNP_temp)
dim(filtered_MAF_Genos_oneSNP_temp)

#Retain the SNP with the highest MAF per tag
# This one gives an output with the correct number of unique tags, and spot checked in excel
oneMAF<- filtered_MAF_Genos_oneSNP_temp %>% group_by(Tag) %>% slice(which.max(MAF))
head(oneMAF)
oneMAF<-as.data.frame(oneMAF)
head(oneMAF)
dim(oneMAF)

#write.table(oneMAF,"Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/OneMAF.txt",quote=FALSE,row.names=FALSE, header=TRUE)

#filter dataset to include only one SNP per tag loci
length(oneMAF$Locus)

filtered_MAF_Genos_oneSNP<-filtered_MAF_Genos[,oneMAF$Locus]
dim(filtered_MAF_Genos_oneSNP)

##I tested the difference between them using columns that I knew included tags with multiple snps
filtered_MAF_Genos[1:5,7:15]
filtered_MAF_Genos_oneSNP[1:5,7:15]

#write.table(filtered_MAF_Genos_oneSNP,"Z:/WORK/TARPEY/ChinookPanel/FilteringRawGenotypes/filtered_MAF_Genos_oneSNP.txt",quote=FALSE,row.names=TRUE, header=TRUE)
