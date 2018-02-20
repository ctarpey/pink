### Initial Filtering of Pink Stacks Population Output
###    combines 29 whitelists and filters to create a final whitelist
### Garrett McKinney and Carolyn Tarpey | October 2017
### ---------------------------------------

#This is R code originally written by Garrett McKinney to combine and filter the output of Stacks populations. 
#It requires a genepop file for each of the whitelists used to run the populations command in Stacks. 
#The genepop files should be stripped of the header line that Stacks puts in and the second line should start with a tab so that 
#R will load it in as a table. 

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
dim(genos12)
genos13<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_13.genepop",colClasses="factor")
genos13[1:5,1:5]
dim(genos13)
genos14<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_14.genepop",colClasses="factor")
genos14[1:5,1:5]
dim(genos14)
genos15<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_15.genepop",colClasses="factor")
genos15[1:5,1:5]
dim(genos15)
genos16<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_16.genepop",colClasses="factor")
genos16[1:5,1:5]
dim(genos16)
genos17<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_17.genepop",colClasses="factor")
genos17[1:5,1:5]
dim(genos17)
genos18<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_18.genepop",colClasses="factor")
genos18[1:5,1:5]
dim(genos18)
genos19<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_19.genepop",colClasses="factor")
genos19[1:5,1:5]
dim(genos19)
genos20<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_20.genepop",colClasses="factor")
genos20[1:5,1:5]
dim(genos20)
genos21<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_21.genepop",colClasses="factor")
genos21[1:5,1:5]
dim(genos21)
genos22<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_22.genepop",colClasses="factor")
genos22[1:5,1:5]
dim(genos22)
genos23<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_23.genepop",colClasses="factor")
genos23[1:5,1:5]
dim(genos23)
genos24<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_24.genepop",colClasses="factor")
genos24[1:5,1:5]
dim(genos24)
genos25<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_25.genepop",colClasses="factor")
genos25[1:5,1:5]
dim(genos25)
genos26<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_26.genepop",colClasses="factor")
genos26[1:5,1:5]
dim(genos26)
genos27<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_27.genepop",colClasses="factor")
genos27[1:5,1:5]
dim(genos27)
genos28<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_28.genepop",colClasses="factor")
genos28[1:5,1:5]
dim(genos28)
genos29<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/whitelist_29.genepop",colClasses="factor")
genos29[1:5,1:5]
dim(genos29)


#use cbind to combine data into single dataset
allGenos<-cbind(genos1,genos2,genos3,genos4,genos5,genos6,genos7,genos8,genos9,genos10,genos11,genos12,genos13,genos14,genos15,genos16,genos17,genos18,genos19,genos20,genos21,genos22,genos23,genos24,genos25,genos26,genos27,genos28,genos29)
allGenos[1:5,1:5]
dim(allGenos)

#get unique tags from full dataset
allLoci<-data.frame(str_split_fixed(colnames(allGenos),"_",2))
colnames(allLoci)<-c("Tag","SNP")
allLoci$Locus<-colnames(allGenos)
head(allLoci)
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
ggplot()+geom_point(data=locusGenoRate_ranked,aes(x=rank,y=GenoRate))+xlab("Locus Genotype Rank")+ylab("Genotype Rate")+theme_bw()+ggtitle("Locus Genotype Rate by Rank")
#ggsave(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/LocusGenotypeRate_all.pdf")
#plot barchart of genotype rate for loci
ggplot()+geom_bar(data=locusGenoRate_ranked,aes(x=GenoRate))+ggtitle("Count of Loci per Genotype Rate ")
#plot histogram of genotype rate for loci
ggplot()+geom_histogram(data=locusGenoRate,aes(x=GenoRate),binwidth=0.001)+ggtitle("Count of Loci per Genotype Rate")
#ggsave(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/LocusCountPerGenotypeRate_all.pdf")


##################### Genotype rate >70%
#filterLoci with >= 70% genotype rate
filteredLoci<-locusGenoRate[locusGenoRate$GenoRate>=0.70,]
dim(filteredLoci)

#get unique tags from 70% genoRate filtered dataset
retainedLoci_70perc<-data.frame(str_split_fixed(filteredLoci$Locus,"_",2))
colnames(retainedLoci_75perc)<-c("Tag","SNP")
retainedLoci_70perc$Locus<-filteredLoci$Locus
head(retainedLoci_70perc)
length(unique(retainedLoci_70perc$Tag))

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
ggplot()+geom_point(data=sampleGenoRate_filteredLoci_ranked,aes(x=rank,y=GenoRate))+geom_hline(aes(yintercept=0.70),lty="dashed")+ggtitle("Sample Genotype Rate of Loci Genotyped in 70% of Samples")
ggsave(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/SampleGenotypeRateofLociGenotypedin70Samples.pdf")
ggplot()+geom_point(data=sampleGenoRate_filteredLoci_ranked,aes(x=rank,y=GenoRate))+xlab("Sample Genotype Rate Rank")+ylab("Genotype Rate")+theme_bw()+ggtitle("Sample Genotype Rate of Loci Genotyped in 70% of Samples")
#plot histogram of genotype rate filter for samples with filtered loci
ggplot()+geom_histogram(data=sampleGenoRate_filteredLoci,aes(x=GenoRate),binwidth=0.01)+ggtitle("Sample Genotype Rate of Loci Genotyped in 70% of Samples")

#filter samples with >=70% genotype rate
filteredSamples<-sampleGenoRate_filteredLoci[sampleGenoRate_filteredLoci$GenoRate>=0.70,]
filteredGenos_filteredSamples<-filteredGenos[rownames(filteredGenos)%in%filteredSamples$Sample,]
dim(filteredGenos_filteredSamples)


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
ggplot()+geom_point(data=locusGenoRate_filteredSamples_ranked,aes(x=rank,y=GenoRate))+ggtitle("Locus Genotype Rate",subtitle= "of Loci Genotyped in 70% of Samples, and Samples Genotyped at 70% of Loci")
#plot barchart of genotype rate for loci
ggplot()+geom_bar(data=locusGenoRate_filteredSamples_ranked,aes(x=GenoRate))+ggtitle("Locus Genotype Rate",subtitle= "of Loci Genotyped in 70% of Samples, and Samples Genotyped at 70% of Loci")
#plot histogram of genotype rate for loci
ggplot()+geom_histogram(data=locusGenoRate_filteredSamples,aes(x=GenoRate),binwidth=0.001)+ggtitle("Locus Genotype Rate",subtitle= "of Loci Genotyped in 70% of Samples, and Samples Genotyped at 70% of Loci")
ggsave(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/LocusGenotypeRateFiltered70and70.pdf")


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

#plot histogram of MAF
ggplot()+geom_histogram(data=filteredGenos_filteredSamples_MAF,aes(x=MAF),binwidth=0.01)+xlab("MAF")+ylab("Frequency")+theme_bw()+ggtitle("Minor Allele Frequency",subtitle= "of Loci Genotyped in 70% of Samples, and Samples Genotyped at 70% of Loci")
ggsave(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/GlobalMAF_filtered70and70.pdf")


#plot genotype rate per locus relative to MAF
#combine MAF and genotype rate summaries
filteredLoci_summary<-merge(locusGenoRate_filteredSamples,filteredGenos_filteredSamples_MAF,by="Locus")
head(filteredLoci_summary)
ggplot()+geom_point(data=filteredLoci_summary,aes(x=GenoRate,y=MAF),alpha=0.1)+ggtitle("Genotype Rate per Locus, Relative to MAF")


#filter loci to retain only those with averall MAF>=0.05
MAF05_loci<-filteredGenos_filteredSamples_MAF[filteredGenos_filteredSamples_MAF$MAF>=0.05,]
filteredGenos_filteredSamples_MAF05<-filteredGenos_filteredSamples[,colnames(filteredGenos_filteredSamples)%in%MAF05_loci$Locus]
dim(filteredGenos_filteredSamples_MAF05)


#get unique tags from MAF 0.05 filtered dataset
retainedLoci_MAF_05<-data.frame(str_split_fixed(colnames(filteredGenos_filteredSamples_MAF05),"_",2))
colnames(retainedLoci_MAF_05)<-c("Tag","SNP")
retainedLoci_MAF_05$Locus<-colnames(filteredGenos_filteredSamples_MAF05)
head(retainedLoci_MAF_05)
length(unique(retainedLoci_MAF_05$Tag))

#get number of SNPs per locus
SNPdistribution_MAF_05<-table(retainedLoci_MAF_05$Tag)
SNPdistribution_MAF_05<-as.data.frame(SNPdistribution_MAF_05)
colnames(SNPdistribution_MAF_05)<-c("Tag","SNPnum")
head(SNPdistribution_MAF_05)
#plot histogram of SNP number
ggplot()+geom_histogram(data=SNPdistribution_MAF_05,aes(x=SNPnum),binwidth=1)+xlab("SNPs per 150bp tag")+ylab("Frequency")+theme_bw()
sum(SNPdistribution_MAF_05$SNPnum>1)/dim(SNPdistribution_MAF_05)[1]

#output list of filtered loci, need to open as binary file (using write binary, "wb") to give unix line endings
#format dataset to remove X from locus names
#filteredLociIDs_MAF_05<-MAF02_loci$Locus
#filteredLociIDs_MAF_05<-gsub("X","",filteredLociIDs_MAF_05)
filteredLociIDs_MAF_05<-select(retainedLoci_MAF_05,Tag,SNP)
filteredLociIDs_MAF_05$Tag<-gsub("X","",filteredLociIDs_MAF_05$Tag)
outputFile<-file("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/70perc_genoRate_MAF05_Loci.txt", "wb")
write.table(filteredLociIDs_MAF_05,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",eol="\n")
close(outputFile)



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
