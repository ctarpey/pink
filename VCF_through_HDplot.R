### Run vcf# through HDplot
###    runs the chum vcf through HDplot- was part of initial filtering that I removed and retained here. 
### Garrett McKinney | October 2017
### ---------------------------------------

#Run vcf# data through HDplot
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