#Chum data filtering
library(ggplot2)
library(vcfR)

#load genepop files
genos1<-read.table("E:/CHUM/RAD/stacks/whitelist_1/batch_17_whitelist1.genepop",colClasses="factor")
genos1[1:5,1:5]
dim(genos1)
genos2<-read.table("E:/CHUM/RAD/stacks/whitelist_2/batch_17_whitelist2.genepop",colClasses="factor")
genos2[1:5,1:5]
dim(genos2)
genos3<-read.table("E:/CHUM/RAD/stacks/whitelist_3/batch_17_whitelist3.genepop",colClasses="factor")
genos3[1:5,1:5]
dim(genos3)
genos4<-read.table("E:/CHUM/RAD/stacks/whitelist_4/batch_17_whitelist4.genepop",colClasses="factor")
genos4[1:5,1:5]
dim(genos4)
genos5<-read.table("E:/CHUM/RAD/stacks/whitelist_5/batch_17_whitelist5.genepop",colClasses="factor")
genos5[1:5,1:5]
dim(genos5)
genos6<-read.table("E:/CHUM/RAD/stacks/whitelist_6/batch_17_whitelist6.genepop",colClasses="factor")
genos6[1:5,1:5]
dim(genos6)

#use cbind to combine data into single dataset
allGenos<-cbind(genos1,genos2,genos3,genos4,genos5,genos6)
allGenos[1:5,1:5]
dim(allGenos)

#get unique tags from full dataset
allLoci<-data.frame(str_split_fixed(colnames(allGenos),"_",2))
colnames(allLoci)<-c("Tag","SNP")
allLoci$Locus<-colnames(allGenos)
head(allLoci)
length(unique(allLoci$Tag))

#get genotype rate per locus
locusGenoRate<-apply(allGenos,2,function(x) 1-(sum(x=="0000")/dim(allGenos)[1]))
locusGenoRate<-data.frame(keyName=names(locusGenoRate), value=locusGenoRate, row.names=NULL)
colnames(locusGenoRate)<-c("Locus","GenoRate")
dim(locusGenoRate)
head(locusGenoRate)
write.table(locusGenoRate,"E:/CHUM/RAD/Filtering/locusGenoRate.txt",quote=FALSE,row.names=FALSE)

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
ggplot()+geom_point(data=locusGenoRate_ranked,aes(x=rank,y=GenoRate))+xlab("Locus Genotype Rate Rank")+ylab("Genotype Rate")+theme_bw()
#plot barchart of genotype rate for loci
ggplot()+geom_bar(data=locusGenoRate_ranked,aes(x=GenoRate))
#plot histogram of genotype rate for loci
ggplot()+geom_histogram(data=locusGenoRate,aes(x=GenoRate),binwidth=0.001)

#filterLoci with <= 50% genotype rate
filteredLoci<-locusGenoRate[locusGenoRate$GenoRate>=0.5,]
dim(filteredLoci)

#get unique tags from 50% genoRate filtered dataset
retainedLoci_50perc<-data.frame(str_split_fixed(filteredLoci$Locus,"_",2))
colnames(retainedLoci_50perc)<-c("Tag","SNP")
retainedLoci_50perc$Locus<-filteredLoci$Locus
head(retainedLoci_50perc)
length(unique(retainedLoci_50perc$Tag))

#output list of filtered loci, need to open as binary file (using write binary, "wb") to give unix line endings
#format dataset to remove X from locus names
filteredLociIDs<-filteredLoci$Locus
filteredLociIDs<-gsub("X","",filteredLociIDs)
outputFile<-file("E:/CHUM/RAD/Filtering/50perc_genoRateLoci.txt", "wb")
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
ggplot()+geom_point(data=sampleGenoRate_filteredLoci_ranked,aes(x=rank,y=GenoRate))+geom_hline(aes(yintercept=0.75),lty="dashed")
ggplot()+geom_point(data=sampleGenoRate_filteredLoci_ranked,aes(x=rank,y=GenoRate))+xlab("Sample Genotype Rate Rank")+ylab("Genotype Rate")+theme_bw()
#plot histogram of genotype rate filter for samples with filtered loci
ggplot()+geom_histogram(data=sampleGenoRate_filteredLoci,aes(x=GenoRate),binwidth=0.01)

#filter samples with <=75% genotype rate
filteredSamples<-sampleGenoRate_filteredLoci[sampleGenoRate_filteredLoci$GenoRate>=0.75,]
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
ggplot()+geom_point(data=locusGenoRate_filteredSamples_ranked,aes(x=rank,y=GenoRate))
#plot barchart of genotype rate for loci
ggplot()+geom_bar(data=locusGenoRate_filteredSamples_ranked,aes(x=GenoRate))
#plot histogram of genotype rate for loci
ggplot()+geom_histogram(data=locusGenoRate_filteredSamples,aes(x=GenoRate),binwidth=0.001)

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
filteredGenos_filteredSamples_MAF

#plot histogram of MAF
ggplot()+geom_histogram(data=filteredGenos_filteredSamples_MAF,aes(x=MAF),binwidth=0.01)+xlab("MAF")+ylab("Frequency")+theme_bw()

#plot genotype rate per locus relative to MAF
#combine MAF and genotype rate summaries
filteredLoci_summary<-merge(locusGenoRate_filteredSamples,filteredGenos_filteredSamples_MAF,by="Locus")
head(filteredLoci_summary)
ggplot()+geom_point(data=filteredLoci_summary,aes(x=GenoRate,y=MAF),alpha=0.1)

#filter loci to retain only those with averall MAF>=0.02
MAF02_loci<-filteredGenos_filteredSamples_MAF[filteredGenos_filteredSamples_MAF$MAF>=0.02,]
filteredGenos_filteredSamples_MAF02<-filteredGenos_filteredSamples[,colnames(filteredGenos_filteredSamples)%in%MAF02_loci$Locus]
dim(filteredGenos_filteredSamples_MAF02)
#69122 SNPs retained

#get unique tags from MAF 0.02 filtered dataset
retainedLoci_MAF_02<-data.frame(str_split_fixed(colnames(filteredGenos_filteredSamples_MAF02),"_",2))
colnames(retainedLoci_MAF_02)<-c("Tag","SNP")
retainedLoci_MAF_02$Locus<-colnames(filteredGenos_filteredSamples_MAF02)
head(retainedLoci_MAF_02)
length(unique(retainedLoci_MAF_02$Tag))

#get number of SNPs per locus
SNPdistribution_MAF_02<-table(retainedLoci_MAF_02$Tag)
SNPdistribution_MAF_02<-as.data.frame(SNPdistribution_MAF_02)
colnames(SNPdistribution_MAF_02)<-c("Tag","SNPnum")
head(SNPdistribution_MAF_02)
#plot histogram of SNP number
ggplot()+geom_histogram(data=SNPdistribution_MAF_02,aes(x=SNPnum),binwidth=1)+xlab("SNPs per 150bp tag")+ylab("Frequency")+theme_bw()
sum(SNPdistribution_MAF_02$SNPnum>1)/dim(SNPdistribution_MAF_02)[1]

#output list of filtered loci, need to open as binary file (using write binary, "wb") to give unix line endings
#format dataset to remove X from locus names
#filteredLociIDs_MAF_02<-MAF02_loci$Locus
#filteredLociIDs_MAF_02<-gsub("X","",filteredLociIDs_MAF_02)
filteredLociIDs_MAF_02<-select(retainedLoci_MAF_02,Tag,SNP)
filteredLociIDs_MAF_02$Tag<-gsub("X","",filteredLociIDs_MAF_02$Tag)
outputFile<-file("E:/CHUM/RAD/Filtering/50perc_genoRate_MAF02_Loci.txt", "wb")
write.table(filteredLociIDs_MAF_02,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",eol="\n")
close(outputFile)


#plot histogram of MAF 0.02 filtered loci
ggplot()+geom_histogram(data=filteredGenos_filteredSamples_MAF[filteredGenos_filteredSamples_MAF$MAF>=0.02,],aes(x=MAF),binwidth=0.01)
#plot genotype rate per locus relative to MAF for MAF 0.02 filtered loci
ggplot()+geom_point(data=filteredLoci_summary[filteredLoci_summary$MAF>=0.02,],aes(x=GenoRate,y=MAF),alpha=0.1)

#filter loci to retain only those with averall MAF>=0.05
MAF05_loci<-filteredGenos_filteredSamples_MAF[filteredGenos_filteredSamples_MAF$MAF>=0.05,]
filteredGenos_filteredSamples_MAF05<-filteredGenos_filteredSamples[,colnames(filteredGenos_filteredSamples)%in%MAF05_loci$Locus]
dim(filteredGenos_filteredSamples_MAF05)
#54,842 SNPs retained


#get unique tags from filtered dataset
retainedLoci<-data.frame(str_split_fixed(colnames(filteredGenos_filteredSamples_MAF05),"_",2))
colnames(retainedLoci)<-c("Tag","SNP")
retainedLoci$Locus<-colnames(filteredGenos_filteredSamples_MAF05)
head(retainedLoci)
#plot histogram of SNP position
#histogram won't work because it's discrete
ggplot()+geom_histogram(data=retainedLoci,aes(x=SNP),binwidth=1)
#order SNP position for plotting
sort((unique(retainedLoci$SNP)))
#retainedLoci$SNP <- factor(retainedLoci$SNP, levels = c(...))
ggplot()+geom_bar(data=retainedLoci,aes(x=SNP))
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

#print as genepop file
printGenepop<-function(populations,genotypes){
  #print 
}

printGenepop(samplePops,filteredGenos_filteredSamples_MAF05)


#Run vcf data through HDplot
chum_vcf<-read.vcfR("E:/CHUM/RAD/stacks/batch_17_full.vcf")

#create table of locus IDs, reference allele, and alternate allele
chum_LocusTable<-as.data.frame(matrix(NA,nrow=dim(chum_vcf@fix)[1],ncol=3))
colnames(chum_LocusTable)<-c("Locus_ID","refAllele","altAllele")
chum_LocusTable$Locus_ID<-chum_vcf@fix[,3]
chum_LocusTable$refAllele<-chum_vcf@fix[,4]
chum_LocusTable$altAllele<-chum_vcf@fix[,5]
chum_LocusTable$alleles<-paste(chum_LocusTable$refAllele,chum_LocusTable$altAllele,sep=",")
chum_LocusTable$ploidy<-2
dim(chum_LocusTable)
#RUN HDPLOT

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

chum_HDPlotResults<-parseVCF_HDPlot(chum_vcf)
#save results
write.table(chum_HDPlotResults,"E:/CHUM/RAD/Filtering/chum_HDplotResults.txt",quote=FALSE,row.names=FALSE)

#plot results of HDplot
ggplot()+geom_point(data=chum_HDPlotResults,aes(x=het_perc,y=z),alpha=0.5)

#filter samples to exclude individuals with very low genotype rate (<0.75)
chum_HDPlot_filteredSamples<-sampleGenoRate_filteredLoci$Sample[sampleGenoRate_filteredLoci$GenoRate>=0.75]
chum_HDPlot_filteredSamples

#rerun HDplot on filtered samples to see if results differ


#speed up HDplot using apply or by matrix manipulation...
HDPlot<-function(vcfData){
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