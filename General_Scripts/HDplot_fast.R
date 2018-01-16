library(vcfR)
library(ggplot2)

vcfInput<-read.vcfR("C:/data.vcf")

HDPlot_fast<-function(vcfData){
  #set up results table
  HDplotTable<-as.data.frame(matrix(NA,nrow=dim(vcfData@gt)[1],ncol=9))
  colnames(HDplotTable)<-c("Locus_ID","depth_a","depth_b","ratio","num_hets","num_samples","het_perc","std","z")
  #format allele reads from vcf data into matrix of comma separated values
  genos<-apply(vcfData@gt[,2:dim(vcfData@gt)[2]], 2, function(x) str_split_fixed(x,":",3)[,1])
  rownames(genos)<-vcfData@fix[,3]
  reads<-apply(vcfData@gt[,2:dim(vcfData@gt)[2]], 2, function(x) str_split_fixed(x,":",4)[,3])
  rownames(reads)<-vcfData@fix[,3] 
  #replace . with 0
  reads<-gsub("\\.","0",reads)
  alleleReads_1<-apply(reads,2,function(x) str_split_fixed(x,",",2)[,1])
  alleleReads_2<-apply(reads,2,function(x) str_split_fixed(x,",",2)[,2])
  #convert to numeric format
  alleleReads_1<-apply(alleleReads_1,2, function(x) as.numeric(x))
  alleleReads_2<-apply(alleleReads_2,2, function(x) as.numeric(x))
  rownames(alleleReads_1)<-vcfData@fix[,3]
  rownames(alleleReads_2)<-vcfData@fix[,3]
  #subset to heterozygous genotypes
  #make genotype matrix where heterozygotes are 1 and other genotypes are 0
  hetMatrix<-genos
  hetMatrix<-apply(hetMatrix,2,function(x) dplyr::recode(x,'0/0'=0,'1/1'=0,'./.'=0,'0/1'=1,'1/0'=1))
  #multiply read count matrices by heterozygote matrix to get read counts for heterozygotes
  alleleReads_1_het<-alleleReads_1*hetMatrix
  alleleReads_2_het<-alleleReads_2*hetMatrix
  #rows are loci and columns are samples
  #sum reads per allele per locus for heterozygous samples
  A_reads<-apply(alleleReads_1_het,1,sum)
  B_reads<-apply(alleleReads_2_het,1,sum)
  totalReads<-A_reads+B_reads
  ratio<-A_reads/totalReads
  std<-sqrt(totalReads*0.5*0.5)
  z<- -(totalReads/2-A_reads)/std
  #get percent heterozygosity for each locus
  numHets<-apply(hetMatrix,1,sum)
  hetPerc<-numHets/dim(hetMatrix)[2]

  #assign results to HDplotTable
  HDplotTable$Locus_ID<-vcfData@fix[,3]
  HDplotTable$depth_a<-A_reads
  HDplotTable$depth_b<-B_reads
  HDplotTable$ratio<-ratio
  HDplotTable$num_hets<-numHets
  HDplotTable$num_samples<-dim(hetMatrix)[2]
  HDplotTable$het_perc<-hetPerc
  HDplotTable$std<-std
  HDplotTable$z<-z
  return(HDplotTable)
}

HDPlotResults<-HDPlot(vcfInput)

ggplot()+geom_point(data=HDPlotResults,aes(x=het_perc,y=z))