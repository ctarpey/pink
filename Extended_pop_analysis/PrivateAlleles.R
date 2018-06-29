### Private alleles analysis
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
beringia<-c("PAMUR10","PAMUR11","PHAYLY09","PHAYLY10","PKUSHI06","PKUSHI07","PNOME91","PNOME94","PTAUY09","PTAUY12")
cascadia<-c("PKOPE91","PKOPE96","PSNOH03","PSNOH96","PLAKEL06","PSPINK14")
susitna <-c("PLAKEL07","PDISAP13")

sampleList$region<-ifelse(sampleList$population %in% beringia, "beringia",ifelse(sampleList$population %in% cascadia, "cascadia", ifelse(sampleList$population %in% susitna, "susitna",NA)))
sampleList$region
sampleList

#write to file
write.table(sampleList,file="Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/pink_genepop_sampleInfo.txt",sep="\t",quote=FALSE,row.names=FALSE)

#These dont have the susitna samples in them- take that last part out if you want to include them
evenSamples<-sampleList$sample[sampleList$lineage=="even" & sampleList$region != "susitna"]
oddSamples<-sampleList$sample[sampleList$lineage=="odd" & sampleList$region != "susitna"]

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
write.table(locusSummary,file="Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/pink_genepop_locusSummary.txt",sep="\t",quote=FALSE,row.names=FALSE)


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
write.table(evenOddAlFreq,file="Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/pink_genepop_alleleFrequencies.txt",sep="\t",quote=FALSE,row.names=FALSE)

#plot histogram of MAF for loci that are fixed in one lineage but not the other
evenFixed<-evenOddAlFreq[(evenOddAlFreq$allele1MAF.x==0) & (evenOddAlFreq$allele1Freq.y !=1) ,]
dim(evenFixed)
hist(evenFixed$allele1MAF.y,xlim=c(0,0.5),ylim=c(0,1500),breaks=seq(0,0.5,by=0.02),main="Even Lineage Fixed Loci",xlab="Odd Lineage Allele Frequency",ylab="Number of Loci")
oddFixed<-evenOddAlFreq[(evenOddAlFreq$allele1MAF.y==0) & (evenOddAlFreq$allele1Freq.x !=1) ,]
dim(oddFixed)
hist(oddFixed$allele1MAF.x,xlim=c(0,0.5),ylim=c(0,1500),breaks=seq(0,0.5,by=0.02),main="Odd Lineage Fixed Loci",xlab="Even Lineage Allele Frequency",ylab="Number of Loci")

#write to file
write.table(evenFixed,file="Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/evenFixed.txt",sep="\t",quote=FALSE,row.names=FALSE)
write.table(oddFixed,file="Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/oddFixed.txt",sep="\t",quote=FALSE,row.names=FALSE)


#plot side by side
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/odd_even_private_allele_frequencies.pdf", width = 9, height = 7)

histlist <- list(evenFixed$allele1MAF.y,oddFixed$allele1MAF.x)
multhist(histlist, xlab="Private Allele Frequency", ylab="Number of Loci")
legend("topright", c("Private Alleles Odd Lineage", "Private Alleles Even Lineage"), col=c("black", "gray"), lwd=10)


dev.off()
#two sample Kolmogorov-Smirnov
ks.test(evenFixed$allele1MAF.y, oddFixed$allele1MAF.x, alternative = c("two.sided"), exact = NULL)

##############################Comparing to the Susitna pops###################

#These dont have the susitna samples in them- take that last part out if you want to include them
evenSamples_susitna<-sampleList$sample[sampleList$lineage=="even" & sampleList$region == "susitna"]
oddSamples_susitna<-sampleList$sample[sampleList$lineage=="odd" & sampleList$region == "susitna"]

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
evenGenotypes_susitna<-pinkGenotypes[,colnames(pinkGenotypes) %in% evenSamples_susitna]
oddGenotypes_susitna<-pinkGenotypes[,colnames(pinkGenotypes) %in% oddSamples_susitna]

#get allele frequencies for each lineage
evenAlFreqs_susitna<-alFreqs(evenGenotypes_susitna)
oddAlFreqs_susitna<-alFreqs(oddGenotypes_susitna)

#plot histogram of allele frequencies for each lineage
hist(evenAlFreqs_susitna$allele1Freq,ylim=c(0,18000),xlab="Allele Frequency",ylab="Number of Loci",main="_susitna Even Year Lineage")
hist(oddAlFreqs_susitna$allele1Freq,ylim=c(0,18000),xlab="Allele Frequency",ylab="Number of Loci",main="_susitna Odd Year Lineage")

#combine allele frequencies for each lineage into one table
evenOddAlFreq_susitna<-merge(evenAlFreqs_susitna,oddAlFreqs_susitna,by="locus")

#plot allele frequencies for each locus with even year allele frequency on the x-axis and
#odd year allele frequency on the y-axis
xyplot(allele1Freq.y ~ allele1Freq.x, data=evenOddAlFreq_susitna)
#get the difference in allele frequency between the lineages and plot histogram
evenOddAlFreq_susitna$difference=evenOddAlFreq_susitna$allele1Freq.x-evenOddAlFreq_susitna$allele1Freq.y
hist(evenOddAlFreq_susitna$difference)

#get loci that are fixed in one lineage and plot the allele frequency in the other lineage
evenOddAlFreq_susitna[1:10,]
#convert allele frequencies to minor allele frequency
for(i in 1:length(evenOddAlFreq_susitna$allele1Freq.x)){
  if(evenOddAlFreq_susitna$allele1Freq.x[i]>0.5){
    evenOddAlFreq_susitna$allele1MAF.x[i]<-1-evenOddAlFreq_susitna$allele1Freq.x[i]
  }else{
    evenOddAlFreq_susitna$allele1MAF.x[i]<-evenOddAlFreq_susitna$allele1Freq.x[i]
  }
  if(evenOddAlFreq_susitna$allele1Freq.y[i]>0.5){
    evenOddAlFreq_susitna$allele1MAF.y[i]<-1-evenOddAlFreq_susitna$allele1Freq.y[i]
  }else{
    evenOddAlFreq_susitna$allele1MAF.y[i]<-evenOddAlFreq_susitna$allele1Freq.y[i]
  }
}
#write to file
write.table(evenOddAlFreq_susitna,file="Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/pink_genepop_alleleFrequencies_susitna.txt",sep="\t",quote=FALSE,row.names=FALSE)

#plot histogram of MAF for loci that are fixed in one lineage but not the other
evenFixed_susitna<-evenOddAlFreq_susitna[(evenOddAlFreq_susitna$allele1MAF.x==0) & (evenOddAlFreq_susitna$allele1Freq.y !=1),]
head(evenOddAlFreq_susitna)
dim(evenFixed_susitna)
hist(evenFixed_susitna$allele1MAF.y,xlim=c(0,0.5),ylim=c(0,1000),breaks=seq(0,0.5,by=0.02),main="_susitna Even Lineage Fixed Loci",xlab="Odd Lineage Allele Frequency",ylab="Number of Loci")
oddFixed_susitna<-evenOddAlFreq_susitna[(evenOddAlFreq_susitna$allele1MAF.y==0) & (evenOddAlFreq_susitna$allele1Freq.x !=1),]
dim(oddFixed_susitna)
hist(oddFixed_susitna$allele1MAF.x,xlim=c(0,0.5),ylim=c(0,1000),breaks=seq(0,0.5,by=0.02),main="_susitnaOdd Lineage Fixed Loci",xlab="Even Lineage Allele Frequency",ylab="Number of Loci")

#write to file
write.table(evenFixed_susitna,file="Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/evenFixed_susitna.txt",sep="\t",quote=FALSE,row.names=FALSE)
write.table(oddFixed_susitna,file="Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/oddFixed_susitna.txt",sep="\t",quote=FALSE,row.names=FALSE)


#plot side by side
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PrivateAlleles/odd_even_private_allele_frequencies_susitna.pdf", width = 9, height = 7)

histlist <- list(evenFixed_susitna$allele1MAF.y,oddFixed_susitna$allele1MAF.x)
multhist(histlist,xlab="Private Allele Frequency",ylab="Number of Loci")
legend("topright", c("Private Alleles Odd", "Private Alleles Even"), col=c("black", "gray"), lwd=10)

dev.off()
#two sample Kolmogorov-Smirnov
ks.test(evenFixed_susitna$allele1MAF.y, oddFixed_susitna$allele1MAF.x, alternative = c("two.sided"), exact = NULL)



############################## VCF and HDPLOT #################################
#split lineages by region and make histograms as above

#Run HDplot
pinkVCF<-read.vcfR("E:/PinkData/whitelist_vcfCombined/batch_3_filtered.vcf")
pinkVCF_even<-read.vcfR("E:/PinkData/whitelist_vcfCombined/evenSamples.vcf")
pinkVCF_odd<-read.vcfR("E:/PinkData/whitelist_vcfCombined/oddSamples.vcf")

#Loop through VCF data to get information all at once
parseVCF_HDPlot<-function(vcfData){
  HDplotTable<-as.data.frame(matrix(NA,nrow=dim(vcfData@gt)[1],ncol=9))
  colnames(HDplotTable)<-c("Locus_ID","depth_a","depth_b","ratio","num_hets","num_samples","het_perc","std","z")
  HDplotTable$Locus_ID<-vcfData@fix[,3] 
  for (i in 1:dim(vcfData@gt)[1]){ 
    locus_ID<-HDplotTable$Locus_ID[i]
    A_reads<-0
    B_reads<-0
    num_hets<-0
    sampleNum<-(dim(vcfData@gt)[2]-1)
    for (j in 2:dim(vcfData@gt)[2]){
      locusInfo<-unlist(strsplit(vcfData@gt[i,j],":",perl=TRUE))
      geno<-locusInfo[[1]]
      reads<-locusInfo[[3]]
      alleles<-unlist(strsplit(geno,"/",perl=TRUE))
      if (alleles[1]!=alleles[2]){
        num_hets<-num_hets+1
        alleleReads<-as.numeric(unlist(strsplit(reads,",",perl=TRUE)))
        A_reads<-A_reads+alleleReads[1]
        B_reads<-B_reads+alleleReads[2]
      }
    }
    totalReads<-A_reads+B_reads
    ratio<-A_reads/totalReads
    hetPerc<-num_hets/sampleNum
    std<-sqrt(totalReads*0.5*0.5)
    z<- -(totalReads/2-A_reads)/std
    HDplotTable$depth_a[i]<-A_reads
    HDplotTable$depth_b[i]<-B_reads
    HDplotTable$ratio[i]<-ratio
    HDplotTable$num_hets[i]<-num_hets
    HDplotTable$num_samples[i]<-sampleNum
    HDplotTable$het_perc[i]<-hetPerc
    HDplotTable$std[i]<-std
    HDplotTable$z[i]<-z
  }
  return(HDplotTable)
}

#get HDplot results from VCF input
pink_HDplotData<-parseVCF_HDPlot(pinkVCF)
pink_even_HDplotData<-parseVCF_HDPlot(pinkVCF_even)
pink_odd_HDplotData<-parseVCF_HDPlot(pinkVCF_odd)

#plot H and D to visualize results and identify paralogs
xyplot(z~het_perc,data=pink_HDplotData,pch=16,alpha=0.2,xlab="H",ylab="D",main="HDplot Pink Salmon")
xyplot(ratio~het_perc,data=pink_HDplotData,pch=16,alpha=0.2,xlab="H",ylab="D",main="HDplot Pink Salmon")

xyplot(z~het_perc,data=pink_even_HDplotData,pch=16,alpha=0.2,xlab="H",ylab="D",main="HDplot Pink Salmon")
xyplot(ratio~het_perc,data=pink_even_HDplotData,pch=16,alpha=0.2,xlab="H",ylab="D",main="HDplot Pink Salmon")

xyplot(z~het_perc,data=pink_odd_HDplotData,pch=16,alpha=0.2,xlab="H",ylab="D",main="HDplot Pink Salmon")
xyplot(ratio~het_perc,data=pink_odd_HDplotData,pch=16,alpha=0.2,xlab="H",ylab="D",main="HDplot Pink Salmon")

#filter HDplot results to include only loci from population study
#need to convert to tags to match entry from HDplot results since SNP position isn't taken from VCF
locusList<-row.names(pinkGenotypes)
tagList<-str_split_fixed(locusList,"_",2)
tagList<-tagList[,1]
tagList<-gsub("X","",tagList)
tagList<-unique(as.numeric(tagList))

popLoci_HDplotData<-pink_HDplotData[pink_HDplotData$Locus_ID %in% tagList,]

xyplot(z~het_perc,data=popLoci_HDplotData,pch=16,alpha=0.2,xlab="H",ylab="D",main="HDplot Pink Salmon")
xyplot(ratio~het_perc,data=popLoci_HDplotData,pch=16,alpha=0.2,xlab="H",ylab="D",main="HDplot Pink Salmon")

#thresholds set visually at H=0.51 and D=7 and D=-7
thresh_H<-0.51
thresh_D<-7
thresh_Dneg<--7
#add lines for thresholds
xyplot(z ~ het_perc, data=popLoci_HDplotData,xlab="H",ylab="D",main="HDplot Pink Salmon",pch=16,alpha=0.2,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)  ## First call the default panel function for 'xyplot'
         panel.abline(v = thresh_H,lwd=2,lty=2)
         panel.abline(h = thresh_D,lwd=2,lty=2)
         panel.abline(h = thresh_Dneg,lwd=2,lty=2)
       })
#Add paralog status to HDplot results table based on thresholds above
paralogStatus<-function(data,thresh_H,thresh_D,thresh_Dneg){
  #if(data$het_perc>thresh_H|data$z>thresh_D|data$z<thresh_Dneg){
  if(is.na(data[9])){
    paralog<-NA
  }else{
    H<-as.numeric(data[7])
    D<-as.numeric(data[9])
    if(H>thresh_H|D>thresh_D|D<thresh_Dneg){
      paralog<-1
    }else{
      paralog<-0
    }
  }
  return(paralog)
}

popLoci_HDplotData$paralog<-apply(popLoci_HDplotData,1,paralogStatus,thresh_H,thresh_D,thresh_Dneg)

#write to file
write.table(popLoci_HDplotData,file="V:/WORK/MCKINNEY/PinkData/PinkGenepop/pink_genepop_HDplot_popLoci.txt",sep="\t",quote=FALSE,row.names=FALSE)

#replot results color coded by paralog status
xyplot(z~het_perc,data=popLoci_HDplotData,groups=paralog,xlab="H",ylab="D",main="HDplot Pink Salmon",pch=16,alpha=0.2)

#set paralog status for each tag to duplicate (1) if any SNP in that tag was identified as a duplicate
tagParalogStatus<-aggregate(paralog ~ Locus_ID, data=popLoci_HDplotData, FUN=max)
#write to file
write.table(tagParalogStatus,file="V:/WORK/MCKINNEY/PinkData/PinkGenepop/pink_genepop_tagParalogStatus.txt",sep="\t",quote=FALSE,row.names=FALSE)
tagParalogStatus$paralog[tagParalogStatus$Locus_ID=="100029"]

#add paralog status to allele frequency data
evenOddAlFreq_Tag<-str_split_fixed(evenOddAlFreq$locus,"_",2)
evenOddAlFreq_Tag<-evenOddAlFreq_Tag[,1]
evenOddAlFreq_Tag<-gsub("X","",evenOddAlFreq_Tag)
evenOddAlFreq$Locus_ID<-evenOddAlFreq_Tag
evenOddAlFreq<-left_join(evenOddAlFreq,tagParalogStatus)
#convert NA entries to 2 to allow identification, are NA entries paralogs that were previously removed?
evenOddAlFreq$paralog[is.na(evenOddAlFreq$paralog)]<-2

#write to file
write.table(evenOddAlFreq,file="V:/WORK/MCKINNEY/PinkData/PinkGenepop/pink_genepop_evenOddAlFreq.txt",sep="\t",quote=FALSE,row.names=FALSE)

#plot allele frequencies for each locus with even year frequencies on the x-axis and odd-year frequencies
#on the y-axis, color code by paralog status.  Blue is singleton, red is duplicate, green is unknown (not in HDplot data)
xyplot(allele1Freq.y ~ allele1Freq.x, data=evenOddAlFreq, group=paralog)

#identify highly divergent loci as those with an allele frequency difference greater than 0.4
#0.4 threshold is based on histogram of differences
hist(evenOddAlFreq$difference)
evenOddAlFreq$divergent<-ifelse(evenOddAlFreq$difference>=0.5 | evenOddAlFreq$difference<=-0.5, 1, 0)
xyplot(allele1Freq.y ~ allele1Freq.x, data=evenOddAlFreq, group=divergent)
sum(evenOddAlFreq$divergent==1)/length(evenOddAlFreq$divergent)


#write to file
write.table(evenOddAlFreq,file="V:/WORK/MCKINNEY/PinkData/PinkGenepop/pink_genepop_evenOddAlFreq.txt",sep="\t",quote=FALSE,row.names=FALSE)

#load Carolyn's linkage map
linkageMap<-read.delim("V:/WORK/MCKINNEY/PinkData/PinkGenepop/CarolynLinkageMap.txt",header=TRUE)
names(linkageMap)
#add divergent locus status to linkage map
linkageMap$divergent<-evenOddAlFreq$divergent[match(linkageMap$Tag,evenOddAlFreq$Locus_ID)]
#plot linkage map, color coding by divergent status
xyplot(Position ~ LG, group=divergent, data=linkageMap,pch=16)
#create violin plot to show density by divergent status
nonDiverged<-subset(linkageMap,divergent==0)
nonDiverged$LG<-as.character(nonDiverged$LG)
diverged<-subset(linkageMap,divergent==1)
diverged$LG<-as.character(diverged$LG)
ggplot()+geom_violin(data=nonDiverged, aes(x=LG,y=Position),fill="blue",color="blue",alpha=0.5,scale="width")+geom_violin(data=diverged, aes(x=LG,y=Position),fill="red",color="red",alpha=0.5,scale="width")+theme_classic()+labs(title="Distribution of Loci by Divergence",x="Chromosome",y="Position (cM)")+theme(axis.text=element_text(size=12),axis.title=element_text(size=16),plot.title=element_text(size=20),axis.text.x = element_text(angle = 90, hjust = 1))


########################## Private Alleles by Region ##########################
#Separate into Beringia and Cascadia even/odd

evenBeringiaSamples<-sampleList$sample[sampleList$lineage=="even" & sampleList$region=="beringia"]
oddBeringiaSamples<-sampleList$sample[sampleList$lineage=="odd" & sampleList$region=="beringia"]

evenCascadiaSamples<-sampleList$sample[sampleList$lineage=="even" & sampleList$region=="cascadia"]
oddCascadiaSamples<-sampleList$sample[sampleList$lineage=="odd" & sampleList$region=="cascadia"]

#split genotypes by lineage
evenBeringiaGenotypes<-pinkGenotypes[,colnames(pinkGenotypes) %in% evenBeringiaSamples]
oddBeringiaGenotypes<-pinkGenotypes[,colnames(pinkGenotypes) %in% oddBeringiaSamples]

evenCascadiaGenotypes<-pinkGenotypes[,colnames(pinkGenotypes) %in% evenCascadiaSamples]
oddCascadiaGenotypes<-pinkGenotypes[,colnames(pinkGenotypes) %in% oddCascadiaSamples]

#get allele frequencies for each lineage
evenBeringiaAlFreqs<-alFreqs(evenBeringiaGenotypes)
oddBeringiaAlFreqs<-alFreqs(oddBeringiaGenotypes)

evenCascadiaAlFreqs<-alFreqs(evenCascadiaGenotypes)
oddCascadiaAlFreqs<-alFreqs(oddCascadiaGenotypes)

#plot histogram of allele frequencies for each lineage
hist(evenBeringiaAlFreqs$allele1Freq)
hist(oddBeringiaAlFreqs$allele1Freq)

hist(evenCascadiaAlFreqs$allele1Freq)
hist(oddCascadiaAlFreqs$allele1Freq)

#combine allele frequencies for each lineage into one table
evenOddBeringiaAlFreq<-merge(evenBeringiaAlFreqs,oddBeringiaAlFreqs,by="locus")
evenOddCascadiaAlFreq<-merge(evenCascadiaAlFreqs,oddCascadiaAlFreqs,by="locus")

#plot allele frequencies for each locus with even year allele frequency on the x-axis and
#odd year allele frequency on the y-axis
xyplot(allele1Freq.y ~ allele1Freq.x, data=evenOddBeringiaAlFreq)
xyplot(allele1Freq.y ~ allele1Freq.x, data=evenOddCascadiaAlFreq)

#get the difference in allele frequency between the lineages and plot histogram
evenOddBeringiaAlFreq$difference=evenOddBeringiaAlFreq$allele1Freq.x-evenOddBeringiaAlFreq$allele1Freq.y
hist(evenOddBeringiaAlFreq$difference)
evenOddCascadiaAlFreq$difference=evenOddCascadiaAlFreq$allele1Freq.x-evenOddCascadiaAlFreq$allele1Freq.y
hist(evenOddCascadiaAlFreq$difference)

