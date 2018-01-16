
### Secondary Filtering of Pink Stacks Population Output
###   This code takes the final genepop file from Stacks and filters the genotypes: 
###    highest MAF SNP at each tag, HWE and filter for individuals
###  Carolyn Tarpey | November 2017
### ---------------------------------------

#This is adapted from R code written by Garrett McKinney to combine and filter the output of Stacks populations. 
#It requires a genepop file stripped of the header line, the second line should start with a tab so that 
#R will load it in as a table and the loci names should be tab deliminated, and there should be no "pop" between the samples. 
#The python script EditGenepopForR.py will do this formatting for a Genepop file that doesn't have one snp per line

#install.packages("vcfR")
#install.packages("stringr")
#install.packages("ggplot2")
#install.packages("lazyeval")
#install.packages("adegenet")
#install.packages("dplyr")
#install.packages("DBI")
#install.packages("tibble")

#Pink data filtering
library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)

#load genepop files
all_newGenos<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/batch_4_31485TAGS_genepop_4R.txt", colClasses="factor")
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
dim(oneMAF)


#####################Test to see that what we got in One SNP per tag by MAF has the 16681 loci that we want

#if we want to see which of our 16681 snps is included at this point
loci_16681 <- readLines("Z:/WORK/TARPEY/Pink_Populations/listof16681LOCI.txt")
head(loci_16681)

#list of all the loci names in oneMAF_temp
oneMAF_temp<- gsub("X","", oneMAF$Locus)
head(oneMAF_temp)

#loci that are in 16662 and 31845 
length(in_16681_and_110433)
in_16662_and_31485 <- intersect(in_16681_and_110433, oneMAF_temp)
length(in_16662_and_31485)

diff_16662_31485 <- setdiff(in_16681_and_110433, oneMAF_temp)
length(diff_16662_31485)

###########################################
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

outputFile<-file("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/oneSNPpertag_subsetgenepop.txt", "wb")
write.table(oneMAF_loci,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)


##################### Run Genepop with the genepop file for HWE and Fis and import the results here for filtering based on those 

#### Take the original Genepop file that was formatted for import at the beginning of this code and make it on snp per line
###  Then subset it with the list that was exported above, and run it in Genepop for HWE
###  Take the output of that and use the perl script get_HWresults_from_genepop.pl to get the results from the INF output file

HWE_table_Pre_Filtered <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/batch_4_31485LOCI_HWE_results.txt" ,
                                     stringsAsFactors = FALSE, row.names= 1)
HWE_pops <-c("AMUR_10", "AMUR_11", "SUSIT_13", "HAYLY_09", "HAYLY_10", "KOPE_91", "KOPE_96", "KUSHI_06", "KUSHI_07",
             "LAKEL_06", "LAKEL_07", "NOME_91", "NOME_94", "SNOH_03", "SNOH_96", "SUSIT_14", "TAUY_09", "TAUY_12")

colnames(HWE_table_Pre_Filtered)<-c("AMUR_10", "AMUR_11", "SUSIT_13", "HAYLY_09", "HAYLY_10", "KOPE_91", "KOPE_96", "KUSHI_06", "KUSHI_07",
                                    "LAKEL_06", "LAKEL_07", "NOME_91", "NOME_94", "SNOH_03", "SNOH_96", "SUSIT_14", "TAUY_09", "TAUY_12")
head(HWE_table_Pre_Filtered)
dim(HWE_table_Pre_Filtered)
nrow(HWE_table_Pre_Filtered)

##### Filter the loci by the HWE
##retain loci that were at least 0.05 in at least 9 of the populations

loci_HWE_blank_test<-vector()
loci_HWE_blank_test<-apply(HWE_table,1,function(x) sum(x <=0.05, na.rm=TRUE)-sum(x =="-", na.rm=TRUE))
loci_HWE_blank_test<-data.frame(keynames=names(loci_HWE_blank_test), value=loci_HWE_blank_test, row.names = NULL)
colnames(loci_HWE_blank_test)<-c("Locus", "failedHWE-blanks")
head(loci_HWE_blank_test)
dim(loci_HWE_blank_test)

#if the loci does not pass the test, delete it from this group
locipassedHWE_test<-vector()
locipassedHWE_test<-subset(loci_HWE_blank_test, loci_HWE_blank_test[,2]<=9)
dim(locipassedHWE_test)
head(locipassedHWE_test)

# #if the loci does not pass the test put it into a list for later
# lociFailedHWE_test<-vector()
# lociFailedHWE_test<-subset(loci_HWE_blank_test, loci_HWE_blank_test[,2]>=9)
# dim(lociFailedHWE_test)
# head(lociFailedHWE_test)
# outputFile<-file("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/Loci_FAILED_HWE.txt", "wb")
# write.table(lociFailedHWE_test[,1],outputFile,quote=FALSE,row.names=TRUE,col.names=FALSE,eol="\n")
# close(outputFile)

# #Retain the Loci that Passed HWE
# HWE_table_Post_Filtered <-  HWE_table_Pre_Filtered[row.names(HWE_table_Pre_Filtered)%in%locipassedHWE_test$Locus, ] #this subsets the matrix by the row names in the list
# dim(HWE_table_Post_Filtered)
# head(HWE_table_Post_Filtered)
 

#####################Test to see that what we got in One SNP per tag by MAF has the 16681 loci that we want

#if we want to see which of our 16681 snps is included at this point
loci_16681 <- readLines("Z:/WORK/TARPEY/Pink_Populations/listof16681LOCI.txt")
head(loci_16681)

#list of all the loci names in oneMAF_temp
locipassedHWE_test_comp<- locipassedHWE_test$Locus
head(locipassedHWE_test_comp)
length(locipassedHWE_test_comp)

#loci that are in 14637 and 30088, the loci that passed HWE
length(in_16662_and_31485)
in_14637_and_30088 <- intersect(in_16662_and_31485, locipassedHWE_test_comp)
length(in_14637_and_30088)

diff_14637_and_30088 <- setdiff(in_16662_and_31485, locipassedHWE_test_comp)
length(diff_14637_and_30088)

###########################################

##################subsample the genotypes based on the Loci that passed HWE
#all_newGenos <- temp_allnewGenos #reset

#format genotype dataset to remove X from locus names
#temp_allnewGenos <- all_newGenos
colnames(all_newGenos)<-gsub("X","",colnames(all_newGenos))
all_newGenos[1:5,1:5]

#filter dataset to include only loci that passed HWE
all_newGenos_HWE_filtered<-all_newGenos[,colnames(all_newGenos)%in%locipassedHWE_test$Locus]
dim(all_newGenos_HWE_filtered)
all_newGenos_HWE_filtered[1:5,1:5]


####Write list of loci that passed all filters so far:
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/Loci_passed_HWE_filtered.txt", "wb")
write.table(colnames(all_newGenos_HWE_filtered),outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)



#get genotype rate per sample for filtered loci
all_newGenos_HWE_inds<-apply(all_newGenos_HWE_filtered,1,function(x) 1-(sum(x=="0000")/dim(all_newGenos_HWE_filtered)[2]))
all_newGenos_HWE_inds<-data.frame(keyName=names(all_newGenos_HWE_inds), value=all_newGenos_HWE_inds, row.names=NULL)
colnames(all_newGenos_HWE_inds)<-c("Sample","GenoRate")
dim(all_newGenos_HWE_inds)
head(all_newGenos_HWE_inds)


#plot ranked genotype rate for samples with filtered loci
all_newGenos_HWE_inds_ranked<-all_newGenos_HWE_inds[order(all_newGenos_HWE_inds$GenoRate),]
all_newGenos_HWE_inds_ranked$rank<-seq(1,dim(all_newGenos_HWE_inds_ranked)[1],by=1)
ggplot()+geom_point(data=all_newGenos_HWE_inds_ranked,aes(x=rank,y=GenoRate))+theme_bw() + 
  geom_hline(aes(yintercept=0.80),lty="dashed")+ggtitle("Sample Genotype Rate One SNP per Tag; Line Shows 80% Genotype Rate")
#ggsave(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/SampleGenotypeRateofLociGenotypedin80PCSamples.pdf")


#Keep samples with >=80% genotype rate
filteredSamples<-all_newGenos_HWE_inds[all_newGenos_HWE_inds$GenoRate>=0.80,]
filteredGenos_filteredSamples<-all_newGenos_HWE_filtered[rownames(all_newGenos_HWE_filtered)%in%filteredSamples$Sample,]
dim(filteredGenos_filteredSamples)


#get genotype rate per locus for filtered sample dataset
locusGenoRate_filteredSamples<-apply(filteredGenos_filteredSamples,2,function(x) 1-(sum(x=="0000")/dim(all_newGenos_HWE_inds)[1]))
locusGenoRate_filteredSamples<-data.frame(keyName=names(locusGenoRate_filteredSamples), value=locusGenoRate_filteredSamples, row.names=NULL)
colnames(locusGenoRate_filteredSamples)<-c("Locus","GenoRate")
dim(locusGenoRate_filteredSamples)
head(locusGenoRate_filteredSamples)

####Write a table of the genotypes/individuals that have passed all the filters so far:
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/FilteredGenos_FilteredInds.txt", "wb")
write.table(filteredGenos_filteredSamples,outputFile,quote=FALSE,row.names=TRUE,col.names=TRUE,eol="\n")
close(outputFile)

####Write a list of the individuals that have passed all the filters so far:
filteredSample_names <- row.names(filteredGenos_filteredSamples)
filteredSample_names <- as.list(filteredSamples$Sample)
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/filteredSample_names.txt", "wb")
write.table(filteredSample_names,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

####Write a list of the individuals that have FAILED the filters so far:
#samples with >=80% genotype rate
FAILED_filteredSamples<-all_newGenos_HWE_inds[all_newGenos_HWE_inds$GenoRate<0.80,]
dim(FAILED_filteredSamples)
head(FAILED_filteredSamples)
FAILED_filteredSample_names <-as.list(FAILED_filteredSamples$Sample)
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/FAILED_filteredSamples.txt", "wb")
write.table(FAILED_filteredSample_names,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)


######## Should we re-filter the loci based on genotype rate now that there are indv missing? NO

#plot ranked genotype rate for samples with filtered loci
x_ranked<-locusGenoRate_filteredSamples[order(locusGenoRate_filteredSamples$GenoRate),]
x_ranked$rank<-seq(1,dim(x_ranked)[1],by=1)
ggplot()+geom_point(data=x_ranked,aes(x=rank,y=GenoRate))+ggtitle("Locus Genotype Rate of Filtered 465 Samples")+theme_bw()
#ggsave(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/LocusGenotypeRateofFiltered465Samples.pdf")

#Which loci are genotyped at less than 80% Genotype rate? 
x_80<-locusGenoRate_filteredSamples[locusGenoRate_filteredSamples$GenoRate<0.80,]
x_80_sort <-x_80[order(x_80$GenoRate, decreasing = FALSE),]
dim(x_80_sort)
head(x_80_sort)




