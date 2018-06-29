### Private alleles analysis_PWS for comparison
###   And other analyses
###
### Written by Garrett McKinney
###  edited by Carolyn Tarpey | June 2018
### ---------------------------------------
install.packages("plotrix")
library(RColorBrewer)
library(ggplot2)
library(colorspace)
library(plyr)
library(colorRamps)
library(stringr)
library(lattice)
library(vcfR)
library(dplyr)
library(plotrix)

###############Import the pop info for the populations (IT IS UPDATED!!)
pop_key <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/POPINFO_LS_susitna.txt", header = TRUE, sep = '\t')
head(pop_key)

################Load Pink data, Genepop that has been converted to bases
pop_key <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/POPINFO_LS_susitna.txt", header = TRUE, sep = '\t')

pinkGenepop<-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/Genepop23759_465INDS_bases_edit.txt",header=TRUE)
locusNames<-pinkGenepop[,1]
pinkGenotypes<-pinkGenepop[,2:dim(pinkGenepop)[2]]
row.names(pinkGenotypes)<-locusNames
pinkGenotypes[1:10,1:10]

#transpose dataset
pinkGenotypes<-t(pinkGenotypes)
pinkGenotypes[1:10,1:10]

#get populations for each sample
sampleList<-data.frame(matrix(NA,nrow=length(colnames(pinkGenotypes)),ncol=3))
colnames(sampleList)<-c("sample","population","lineage")
sampleList$sample<-colnames(pinkGenotypes)                    
names(sampleList)
for(i in 1:dim(sampleList)[1]){
  population<-unlist(strsplit(sampleList$sample[i], "_"))[1]
  sampleList$population[i]<-population
  year<-as.numeric(substr(population, nchar(population)-1, nchar(population)))
  if(population == "PLAKEL07"){
    sampleList$lineage[i]<-"even"
  }else if(population == "PLAKEL06"){
    sampleList$lineage[i]<-"odd"
  }else if(year%%2==0){
    sampleList$lineage[i]<-"even"
  }else{
    sampleList$lineage[i]<-"odd"
  }
}
names(sampleList)
sampleList

#assign regions for each sample
beringia<-c("PAMUR10","PAMUR11","PKUSHI06","PKUSHI07","PNOME91","PNOME94","PTAUY09","PTAUY12")
cascadia<-c("PSNOH03","PSNOH96","PLAKEL06","PSPINK14","PKOPE91","PKOPE96")
susitna <-c("PLAKEL07","PDISAP13")
KAM <-c("PHAYLY09","PHAYLY10")

sampleList$region<-ifelse(sampleList$population %in% beringia, "beringia",ifelse(sampleList$population %in% cascadia, "cascadia", ifelse(sampleList$population %in% susitna, "susitna", ifelse(sampleList$population %in% KAM, "KAM",NA))))
sampleList$region
sampleList

#write to file
#write.table(sampleList,file="Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/pink_genepop_sampleInfo.txt",sep="\t",quote=FALSE,row.names=FALSE)

#These dont have the susitna samples in them- take that last part out if you want to include them
evenSamples<-sampleList$sample[sampleList$lineage=="even" & sampleList$region != "susitna" & sampleList$region != "KAM"]
oddSamples<-sampleList$sample[sampleList$lineage=="odd" & sampleList$region != "susitna" & sampleList$region != "KAM"]

alFreqs<-function(pinkGenotypes){
  #make reference table of loci
  locusAlleles<-data.frame(matrix(NA,nrow=dim(pinkGenotypes)[1],ncol=3))
  colnames(locusAlleles)<-c("locus","allele1","allele2")
  locusAlleles$locus<-rownames(pinkGenotypes)
  for(i in 1:dim(locusAlleles)[1]){
    uniqueGenos<-unique(pinkGenotypes[i,]) 
    uniqueGenos<-uniqueGenos[uniqueGenos != "-/-"]
    alleles<-unique(na.omit(unlist(strsplit(uniqueGenos, "[^a-zA-Z]+"))))
    locusAlleles$allele1[i]<-alleles[1]
    locusAlleles$allele2[i]<-alleles[2]
  }
  locusAlleles
  
  #calculate allele frequency per locus
  locusSummary<-data.frame(matrix(NA,nrow=dim(pinkGenotypes)[1],ncol=2))
  colnames(locusSummary)<-c("locus","allele1Freq")
  locusSummary$locus<-rownames(pinkGenotypes)
  for (i in 1:dim(pinkGenotypes)[1]){
    #use ParalogLocusTable to get alleles for each locus
    allele1<-locusAlleles$allele1[i]
    allele2<-locusAlleles$allele2[i]
    allele1Count<-sum(str_count(pinkGenotypes[i,],allele1))
    allele2Count<-sum(str_count(pinkGenotypes[i,],allele2))
    if(is.na(allele1Count)){
      allele1Count=0
    }
    if(is.na(allele2Count)){
      allele2Count=0
    }
    freq<-allele1Count/(allele1Count+allele2Count)
    locusSummary$allele1Freq[i]<-freq
  }
  return(locusSummary)
}

locusSummary<-alFreqs(pinkGenotypes)
names(locusSummary)
hist(locusSummary$allele1Freq)

#write to file
#write.table(locusSummary,file="Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/pink_genepop_locusSummary.txt",sep="\t",quote=FALSE,row.names=FALSE)


#split genotypes by lineage
evenGenotypes<-pinkGenotypes[,colnames(pinkGenotypes) %in% evenSamples]
oddGenotypes<-pinkGenotypes[,colnames(pinkGenotypes) %in% oddSamples]

#get allele frequencies for each lineage
evenAlFreqs<-alFreqs(evenGenotypes)
oddAlFreqs<-alFreqs(oddGenotypes)

#plot histogram of allele frequencies for each lineage
hist(evenAlFreqs$allele1Freq,ylim=c(0,18000),xlab="Allele Frequency",ylab="Number of Loci",main="Even Year Lineage")
hist(oddAlFreqs$allele1Freq,ylim=c(0,18000),xlab="Allele Frequency",ylab="Number of Loci",main="Odd Year Lineage")
#try side by side histogram

#combine allele frequencies for each lineage into one table
evenOddAlFreq<-merge(evenAlFreqs,oddAlFreqs,by="locus")
#plot allele frequencies for each locus with even year allele frequency on the x-axis and
#odd year allele frequency on the y-axis
xyplot(allele1Freq.y ~ allele1Freq.x, data=evenOddAlFreq)
#get the difference in allele frequency between the lineages and plot histogram
evenOddAlFreq$difference = evenOddAlFreq$allele1Freq.x-evenOddAlFreq$allele1Freq.y
hist(evenOddAlFreq$difference)

#get loci that are fixed in one lineage and plot the allele frequency in the other lineage
evenOddAlFreq[1:10,]
#convert allele frequencies to minor allele frequency
for(i in 1:length(evenOddAlFreq$allele1Freq.x)){
  if(evenOddAlFreq$allele1Freq.x[i]>0.5){
    evenOddAlFreq$allele1MAF.x[i]<-1-evenOddAlFreq$allele1Freq.x[i]
  }else{
    evenOddAlFreq$allele1MAF.x[i]<-evenOddAlFreq$allele1Freq.x[i]
  }
  if(evenOddAlFreq$allele1Freq.y[i]>0.5){
    evenOddAlFreq$allele1MAF.y[i]<-1-evenOddAlFreq$allele1Freq.y[i]
  }else{
    evenOddAlFreq$allele1MAF.y[i]<-evenOddAlFreq$allele1Freq.y[i]
  }
}
#write to file
write.table(evenOddAlFreq,file="Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/pink_genepop_alleleFrequencies_NOKAM_NS.txt",sep="\t",quote=FALSE,row.names=FALSE)

#plot histogram of MAF for loci that are fixed in one lineage but not the other
evenFixed<-evenOddAlFreq[(evenOddAlFreq$allele1MAF.x==0) & (evenOddAlFreq$allele1Freq.y !=1) ,]
dim(evenFixed)
hist(evenFixed$allele1MAF.y,xlim=c(0,0.5),ylim=c(0,1500),breaks=seq(0,0.5,by=0.02),main="Even Lineage Fixed Loci",xlab="Odd Lineage Allele Frequency",ylab="Number of Loci")
oddFixed<-evenOddAlFreq[(evenOddAlFreq$allele1MAF.y==0) & (evenOddAlFreq$allele1Freq.x !=1) ,]
dim(oddFixed)
hist(oddFixed$allele1MAF.x,xlim=c(0,0.5),ylim=c(0,1500),breaks=seq(0,0.5,by=0.02),main="Odd Lineage Fixed Loci",xlab="Even Lineage Allele Frequency",ylab="Number of Loci")

#write to file
write.table(evenFixed,file="Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/evenFixed_NOKAM_NS.txt",sep="\t",quote=FALSE,row.names=FALSE)
write.table(oddFixed,file="Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/oddFixed_NOKAM_NS.txt",sep="\t",quote=FALSE,row.names=FALSE)


#plot side by side
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/odd_even_private_allele_frequencies.pdf", width = 9, height = 7)

histlist <- list(evenFixed$allele1MAF.y,oddFixed$allele1MAF.x)
multhist(histlist,xlab="Private Allele Frequency",ylab="Number of Loci")
legend("topright", c("Private Alleles Odd", "Private Alleles Even"), col=c("black", "gray"), lwd=10)

dev.off()
#two sample Kolmogorov-Smirnov
ks.test(evenFixed$allele1MAF.y, oddFixed$allele1MAF.x, alternative = c("two.sided"), exact = NULL)

##############################Comparing to the KAM  pops###################

#These dont have the susitna samples in them- take that last part out if you want to include them
evenSamples_KAM<-sampleList$sample[sampleList$lineage=="even" & sampleList$region == "KAM"]
oddSamples_KAM<-sampleList$sample[sampleList$lineage=="odd" & sampleList$region == "KAM"]

alFreqs<-function(pinkGenotypes){
  #make reference table of loci
  locusAlleles<-data.frame(matrix(NA,nrow=dim(pinkGenotypes)[1],ncol=3))
  colnames(locusAlleles)<-c("locus","allele1","allele2")
  locusAlleles$locus<-rownames(pinkGenotypes)
  for(i in 1:dim(locusAlleles)[1]){
    uniqueGenos<-unique(pinkGenotypes[i,]) 
    uniqueGenos<-uniqueGenos[uniqueGenos != "-/-"]
    alleles<-unique(na.omit(unlist(strsplit(uniqueGenos, "[^a-zA-Z]+"))))
    locusAlleles$allele1[i]<-alleles[1]
    locusAlleles$allele2[i]<-alleles[2]
  }
  locusAlleles
  
  #calculate allele frequency per locus
  locusSummary<-data.frame(matrix(NA,nrow=dim(pinkGenotypes)[1],ncol=2))
  colnames(locusSummary)<-c("locus","allele1Freq")
  locusSummary$locus<-rownames(pinkGenotypes)
  for (i in 1:dim(pinkGenotypes)[1]){
    #use ParalogLocusTable to get alleles for each locus
    allele1<-locusAlleles$allele1[i]
    allele2<-locusAlleles$allele2[i]
    allele1Count<-sum(str_count(pinkGenotypes[i,],allele1))
    allele2Count<-sum(str_count(pinkGenotypes[i,],allele2))
    if(is.na(allele1Count)){
      allele1Count=0
    }
    if(is.na(allele2Count)){
      allele2Count=0
    }
    freq<-allele1Count/(allele1Count+allele2Count)
    locusSummary$allele1Freq[i]<-freq
  }
  return(locusSummary)
}

locusSummary<-alFreqs(pinkGenotypes)
names(locusSummary)
hist(locusSummary$allele1Freq)

#split genotypes by lineage
evenGenotypes_KAM<-pinkGenotypes[,colnames(pinkGenotypes) %in% evenSamples_KAM]
oddGenotypes_KAM<-pinkGenotypes[,colnames(pinkGenotypes) %in% oddSamples_KAM]

#get allele frequencies for each lineage
evenAlFreqs_KAM<-alFreqs(evenGenotypes_KAM)
oddAlFreqs_KAM<-alFreqs(oddGenotypes_KAM)

#plot histogram of allele frequencies for each lineage
hist(evenAlFreqs_KAM$allele1Freq,ylim=c(0,18000),xlab="Allele Frequency",ylab="Number of Loci",main="_KAM Even Year Lineage")
hist(oddAlFreqs_KAM$allele1Freq,ylim=c(0,18000),xlab="Allele Frequency",ylab="Number of Loci",main="_KAM Odd Year Lineage")

#combine allele frequencies for each lineage into one table
evenOddAlFreq_KAM<-merge(evenAlFreqs_KAM,oddAlFreqs_KAM,by="locus")

#plot allele frequencies for each locus with even year allele frequency on the x-axis and
#odd year allele frequency on the y-axis
xyplot(allele1Freq.y ~ allele1Freq.x, data=evenOddAlFreq_KAM)
#get the difference in allele frequency between the lineages and plot histogram
evenOddAlFreq_KAM$difference=evenOddAlFreq_KAM$allele1Freq.x-evenOddAlFreq_KAM$allele1Freq.y
hist(evenOddAlFreq_KAM$difference)

#get loci that are fixed in one lineage and plot the allele frequency in the other lineage
evenOddAlFreq_KAM[1:10,]
#convert allele frequencies to minor allele frequency
for(i in 1:length(evenOddAlFreq_KAM$allele1Freq.x)){
  if(evenOddAlFreq_KAM$allele1Freq.x[i]>0.5){
    evenOddAlFreq_KAM$allele1MAF.x[i]<-1-evenOddAlFreq_KAM$allele1Freq.x[i]
  }else{
    evenOddAlFreq_KAM$allele1MAF.x[i]<-evenOddAlFreq_KAM$allele1Freq.x[i]
  }
  if(evenOddAlFreq_KAM$allele1Freq.y[i]>0.5){
    evenOddAlFreq_KAM$allele1MAF.y[i]<-1-evenOddAlFreq_KAM$allele1Freq.y[i]
  }else{
    evenOddAlFreq_KAM$allele1MAF.y[i]<-evenOddAlFreq_KAM$allele1Freq.y[i]
  }
}
#write to file
write.table(evenOddAlFreq_KAM,file="Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/pink_genepop_alleleFrequencies_KAM.txt",sep="\t",quote=FALSE,row.names=FALSE)

#plot histogram of MAF for loci that are fixed in one lineage but not the other
evenFixed_KAM<-evenOddAlFreq_KAM[(evenOddAlFreq_KAM$allele1MAF.x==0) & (evenOddAlFreq_KAM$allele1Freq.y !=1),]
head(evenOddAlFreq_KAM)
dim(evenFixed_KAM)
hist(evenFixed_KAM$allele1MAF.y,xlim=c(0,0.5),ylim=c(0,1000),breaks=seq(0,0.5,by=0.02),main="_KAM Even Lineage Fixed Loci",xlab="Odd Lineage Allele Frequency",ylab="Number of Loci")
oddFixed_KAM<-evenOddAlFreq_KAM[(evenOddAlFreq_KAM$allele1MAF.y==0) & (evenOddAlFreq_KAM$allele1Freq.x !=1),]
dim(oddFixed_KAM)
hist(oddFixed_KAM$allele1MAF.x,xlim=c(0,0.5),ylim=c(0,1000),breaks=seq(0,0.5,by=0.02),main="_KAMOdd Lineage Fixed Loci",xlab="Even Lineage Allele Frequency",ylab="Number of Loci")

#write to file
write.table(evenFixed_KAM,file="Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/evenFixed_KAM.txt",sep="\t",quote=FALSE,row.names=FALSE)
write.table(oddFixed_KAM,file="Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/oddFixed_KAM.txt",sep="\t",quote=FALSE,row.names=FALSE)


#plot side by side
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/odd_even_private_allele_frequencies_KAM.pdf", width = 9, height = 7)

histlist <- list(evenFixed_KAM$allele1MAF.y,oddFixed_KAM$allele1MAF.x)
multhist(histlist,xlab="Private Allele Frequency",ylab="Number of Loci")
legend("topright", c("Private Alleles Odd", "Private Alleles Even"), col=c("black", "gray"), lwd=10)

dev.off()
#two sample Kolmogorov-Smirnov
ks.test(evenFixed_KAM$allele1MAF.y, oddFixed_KAM$allele1MAF.x, alternative = c("two.sided"), exact = NULL)
