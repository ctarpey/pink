### Plot Population Exploratory Statistics for extended pink pops
###    Genotype rate: Locus and Individual, pre and post filtering
###   MAF, HWE
###   Carolyn Tarpey | January 2018
### ---------------------------------------


#Pink data filtering
library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)
library(plotly)
library(gridExtra)
library(scales) 

############ Importing the data that we need to calculate the per population statistics 

#import the original dataset that has only been filtered for 80% genotype rate
OG_80PCfiltered <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/FirstFiltering/filteredGenos_just80PCgenorate_76762.txt",colClasses="factor", header = TRUE)
OG_80PCfiltered[1:5,1:5]
dim(OG_80PCfiltered)

#import the original dataset that has only been filtered for 80% genotype rate
Post_filtered <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/ThirdFiltering/Singleton_One_Tag_Genotypes.txt",colClasses="factor")
Post_filtered[1:5,1:5]
dim(Post_filtered)

#assign each sample to a population, delete the first column (duplicate sample name) and rename the new column Pop
#popmap465<- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/PopMap465_NEWNAMES.txt")
popmap465_NEW<- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/PopMap465_NEWNAMES_NEWSAMPLES.txt")
popmap492<- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/PopMap_NEWNAMES.txt")
head(popmap492)

#add the populations to PostFiltered genotypes
Post_filtered_popz<-cbind(popmap465_NEW,Post_filtered)
Post_filtered_popz$V1<- NULL
colnames(Post_filtered_popz)[1]<-"Pop"
Post_filtered_popz[1:5,1:5]
dim(Post_filtered_popz)

#add the populations to OG_80PCfiltered genotypes
OG_80PCfiltered_popz<-cbind(popmap492$V2,OG_80PCfiltered)
OG_80PCfiltered_popz[,2]<- NULL
colnames(OG_80PCfiltered_popz)[1]<-"Pop"
OG_80PCfiltered_popz[1:5,1:5]
dim(OG_80PCfiltered_popz)



################################## Statistics Per Data Set


### Calculate Individual Genotype Rate per Population 

###OG_filtered
OG_80PCfiltered_IGR_temp<-apply(OG_80PCfiltered_popz[,-1],1,function(x) 1-(sum(x=="0000")/dim(OG_80PCfiltered_popz[,-1])[2]))
OG_80PCfiltered_popz_IGR<-cbind(OG_80PCfiltered_IGR_temp,OG_80PCfiltered_popz)
colnames(OG_80PCfiltered_popz_IGR)[1]<- "SampleGenoRate"
dim(OG_80PCfiltered_popz_IGR)
OG_80PCfiltered_popz_IGR[1:5,1:5]


########
###OG_filtered
Post_filtered_popz_IGR_temp<-apply(Post_filtered_popz[,-1],1,function(x) 1-(sum(x=="0000")/dim(Post_filtered_popz[,-1])[2]))
Post_filtered_popz_IGR<-cbind(Post_filtered_popz_IGR_temp,Post_filtered_popz)
colnames(Post_filtered_popz_IGR)[1]<- "SampleGenoRate"
dim(Post_filtered_popz_IGR)
Post_filtered_popz_IGR[1:5,1:5]


######## Calculate Locus Genotype Rate per Poulation 
#OG_80PCfiltered_popz[1:5,1:5]
###OG_filtered
#create a dataframe for our Locus Geno Rate results
npops<- length(unique(OG_80PCfiltered_popz$Pop))
nloci<- ncol(OG_80PCfiltered_popz[-1])

OG_80PCfiltered_popz_LGR <- matrix(nrow=npops, ncol=nloci) 
rownames(OG_80PCfiltered_popz_LGR)<-as.vector(unique(OG_80PCfiltered_popz$Pop)) #name the rows by the population names
colnames(OG_80PCfiltered_popz_LGR)<-colnames(OG_80PCfiltered_popz[,-1])

for (i in 1:length(unique(OG_80PCfiltered_popz$Pop))){
  tempPop<-(unique(OG_80PCfiltered_popz$Pop)[i])
  tempGeno<-subset(OG_80PCfiltered_popz, Pop == tempPop)
  tempLGR<-apply(tempGeno[,-1],2,function(x) 1-(sum(x=="0000")/dim(tempGeno)[1]))
  OG_80PCfiltered_popz_LGR[i,]<-c(tempLGR)
}

OG_80PCfiltered_popz_LGR[1:5,1:5]
dim(OG_80PCfiltered_popz_LGR)


########

###POST_Filtered_Popz
#create a dataframe for our Locus Geno Rate results
npops<- length(unique(Post_filtered_popz$Pop))
nloci<- length(Post_filtered_popz[,-1])

Post_filtered_popz_LGR <- matrix(nrow=npops, ncol=nloci) 
rownames(Post_filtered_popz_LGR)<-as.vector(unique(Post_filtered_popz$Pop)) #name the rows by the population names
colnames(Post_filtered_popz_LGR)<-colnames(Post_filtered_popz[,-1])

for (i in 1:length(unique(Post_filtered_popz$Pop))){
  tempPop<-(unique(Post_filtered_popz$Pop)[i])
  tempGeno<-subset(Post_filtered_popz, Pop == tempPop)
  tempLGR<-apply(tempGeno[,-1],2,function(x) 1-(sum(x=="0000")/dim(tempGeno)[1]))
  Post_filtered_popz_LGR[i,]<-c(tempLGR)
}

Post_filtered_popz_LGR[1:5,1:5]
dim(Post_filtered_popz_LGR)

#Post_filtered_popz_LGR[,5]

################# Calculate MAF per population 

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

OG_80PCfiltered_popz[1:5,1:5]
Post_filtered_popz[1:5,1:5]
 
######OG_filtered_popz

#create a dataframe for our MAF results
npops<- length(unique(OG_80PCfiltered_popz$Pop))
nloci<- length(OG_80PCfiltered_popz[,-1])

OG_80PCfiltered_popz_MAF<-matrix(nrow=npops, ncol=nloci) 
rownames(OG_80PCfiltered_popz_MAF)<-as.vector(unique(OG_80PCfiltered_popz$Pop)) #name the rows by the population names
colnames(OG_80PCfiltered_popz_MAF)<-colnames(OG_80PCfiltered_popz[,-1])

for (i in 1:length(unique(OG_80PCfiltered_popz$Pop))){
  tempPop<-(unique(OG_80PCfiltered_popz$Pop)[i])
  tempGeno<-subset(OG_80PCfiltered_popz, Pop == tempPop)
  tempMAF<-apply(tempGeno[,-1],2,calculateMAF)
  OG_80PCfiltered_popz_MAF[i,]<-c(tempMAF)
}

OG_80PCfiltered_popz_MAF[1:5,1:5]


###POST_Filtered_Popz

#create a dataframe for our MAF results
npops<- length(unique(Post_filtered_popz$Pop))
nloci<- length(Post_filtered_popz[,-1])

Post_filtered_popz_MAF<-matrix(nrow=npops, ncol=nloci) 
rownames(Post_filtered_popz_MAF)<-as.vector(unique(Post_filtered_popz$Pop)) #name the rows by the population names
colnames(Post_filtered_popz_MAF)<-colnames(Post_filtered_popz[,-1])

for (i in 1:length(unique(Post_filtered_popz$Pop))){
  tempPop<-(unique(Post_filtered_popz$Pop)[i])
  tempGeno<-subset(Post_filtered_popz, Pop == tempPop)
  tempMAF<-apply(tempGeno[,-1],2,calculateMAF)
  Post_filtered_popz_MAF[i,]<-c(tempMAF)
}

Post_filtered_popz_MAF[1:5,1:5]



##################### Run Genepop with the genepop file for HWE and Fis and import the results here ////// This was taken from SecondPinkFiltering.R

HWE_table_Pre_Filtered <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/batch_4_31485LOCI_HWE_results.txt" ,
                                     stringsAsFactors = FALSE, row.names= 1)

colnames(HWE_table_Pre_Filtered)<-c("AMUR_10", "AMUR_11", "SUSIT_13", "HAYLY_09", "HAYLY_10", "KOPE_91", "KOPE_96", "KUSHI_06", "KUSHI_07",
                                    "LAKEL_06", "LAKEL_07", "NOME_91", "NOME_94", "SNOH_03", "SNOH_96", "SUSIT_14", "TAUY_09", "TAUY_12")
head(HWE_table_Pre_Filtered)
dim(HWE_table_Pre_Filtered)

##### Filter the loci by the HWE; retain loci that were at least 0.05 in at least 9 of the populations

loci_HWE_blank_test<-vector()
loci_HWE_blank_test<-apply(HWE_table_Pre_Filtered,1,function(x) sum(x <=0.05, na.rm=TRUE)-sum(x =="-", na.rm=TRUE))
loci_HWE_blank_test<-data.frame(keynames=names(loci_HWE_blank_test), value=loci_HWE_blank_test, row.names = NULL)
colnames(loci_HWE_blank_test)<-c("Locus", "failedHWE-blanks")
head(loci_HWE_blank_test)
dim(loci_HWE_blank_test)

#if the loci does not pass the test, delete it from this group
locipassedHWE_test<-vector()
locipassedHWE_test<-subset(loci_HWE_blank_test, loci_HWE_blank_test[,2]<=9)
dim(locipassedHWE_test)
head(locipassedHWE_test)

#Retain the Loci that Passed HWE
HWE_table_Post_Filtered_test <-  HWE_table_Pre_Filtered[row.names(HWE_table_Pre_Filtered)%in%locipassedHWE_test$Locus, ] #this subsets the matrix by the row names in the list
dim(HWE_table_Post_Filtered_test)
head(HWE_table_Post_Filtered_test)

############ Heterozygosity per Population
## This is only run on the Completely filtered 23,759 loci and it was calculated with a perl script written by garrett: countHets_genepop.pl

HetCounts <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/Heterozygosity/singleton_het_counts.txt" ,
                                     stringsAsFactors = FALSE, row.names= 1,header= TRUE)
head(HetCounts)
dim(HetCounts)

#assign each sample to a population, and rename the new column Pop
popmap465_NEW<- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/PopMap465_NEWNAMES_NEWSAMPLES.txt")

#add the populations to PostFiltered genotypes
HetCounts_popz<-cbind(popmap465_NEW,HetCounts)
HetCounts_popz$V1<- NULL
colnames(HetCounts_popz)[1]<-"Pop"
head(HetCounts_popz)
dim(HetCounts_popz)

outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/Heterozygosity/singleton_het_counts_POPZ.txt", "wb")
write.table(HetCounts_popz,outputFile,quote=FALSE,row.names=TRUE,col.names=TRUE,eol="\n")
close(outputFile)


############################################## PLOTTING The statistics per population  


##Show the colors used for a default palette of 18
show_col(hue_pal()(18))

# this code is from Eleni, she used it to plot something similar in her data. 
#This code requires a dataframe with the individuals and a column that has population designation
#ggplot(mydata, aes(x=loci, fill= population)) + geom_histogram(data= mydata, bins = 20) + facet_wrap(~population) 

####################### Individual Genotype Rate

# plot Individual Genotype Rate Post_Filtered_data
ggplot(Post_filtered_popz_IGR, aes(x=SampleGenoRate, fill= Pop)) + 
  geom_histogram(data= Post_filtered_popz_IGR, bins = 20) + facet_wrap(~Pop) + theme_bw()  + guides(fill= FALSE) +
  ggtitle("Counts of Post-filtered Genotype Rate per Individual by Population")

# plot Individual Genotype Rate OG_data
ggplot(OG_80PCfiltered_popz_IGR, aes(x=SampleGenoRate, fill= Pop)) + 
  geom_histogram(data= OG_80PCfiltered_popz_IGR, bins = 20) + facet_wrap(~Pop) + theme_bw() + guides(fill= FALSE) +
  ggtitle("Counts of Pre-filtered Genotype Rate per Individual by Population")

####################### Locus Genotype Rate


#### OG data Populations 
OG_80PCfiltered_popz_LGR[1:5,1:5]

#Subset the matrix to get each population's numbers
a <-as.data.frame(OG_80PCfiltered_popz_LGR[1,] )
b <-as.data.frame(OG_80PCfiltered_popz_LGR[2,] )
c <-as.data.frame(OG_80PCfiltered_popz_LGR[3,] )
d <-as.data.frame(OG_80PCfiltered_popz_LGR[4,] )
e <-as.data.frame(OG_80PCfiltered_popz_LGR[5,] )
f <-as.data.frame(OG_80PCfiltered_popz_LGR[6,] )
g <-as.data.frame(OG_80PCfiltered_popz_LGR[7,] )
h <-as.data.frame(OG_80PCfiltered_popz_LGR[8,] )
i <-as.data.frame(OG_80PCfiltered_popz_LGR[9,] )
j <-as.data.frame(OG_80PCfiltered_popz_LGR[10,] )
k <-as.data.frame(OG_80PCfiltered_popz_LGR[11,] )
l <-as.data.frame(OG_80PCfiltered_popz_LGR[12,] )
m <-as.data.frame(OG_80PCfiltered_popz_LGR[13,] )
n <-as.data.frame(OG_80PCfiltered_popz_LGR[14,] )
o <-as.data.frame(OG_80PCfiltered_popz_LGR[15,] )
p <-as.data.frame(OG_80PCfiltered_popz_LGR[16,] )
q <-as.data.frame(OG_80PCfiltered_popz_LGR[17,] )
r <-as.data.frame(OG_80PCfiltered_popz_LGR[18,] )

#get the population names in order
row.names(OG_80PCfiltered_popz_LGR)

#assign each of the populations a ggplot to be called later
a1<- ggplot(a, aes(x=a)) + geom_histogram(data= a, bins = 50, fill="#F8766D") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("AMUR_10") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
b1<- ggplot(b, aes(x=b)) + geom_histogram(data= b, bins = 50, fill="#E88526") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("AMUR_11") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
c1<- ggplot(c, aes(x=c)) + geom_histogram(data= c, bins = 50, fill="#D39200") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SUSIT_13") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
d1<- ggplot(d, aes(x=d)) + geom_histogram(data= d, bins = 50, fill="#B79F00") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("HAYLY_09") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
e1<- ggplot(e, aes(x=e)) + geom_histogram(data= e, bins = 50, fill="#93AA00") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("HAYLY_10") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
f1<- ggplot(f, aes(x=f)) + geom_histogram(data= f, bins = 50, fill="#5EB300") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KOPE_91") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g1<- ggplot(g, aes(x=g)) + geom_histogram(data= g, bins = 50, fill="#00BA38") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KOPE_96") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
h1<- ggplot(h, aes(x=h)) + geom_histogram(data= h, bins = 50, fill="#00BF74") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KUSHI_06") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
i1<- ggplot(i, aes(x=i)) + geom_histogram(data= i, bins = 50, fill="#00C19F") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KUSHI_07") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
j1<- ggplot(j, aes(x=j)) + geom_histogram(data= j, bins = 50, fill="#00BFC4") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("LAKEL_06") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
k1<- ggplot(k, aes(x=k)) + geom_histogram(data= k, bins = 50, fill="#00B9C3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("LAKEL07") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
l1<- ggplot(l, aes(x=l)) + geom_histogram(data= l, bins = 50, fill="#00ADFA") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("NOME_91") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m1<- ggplot(m, aes(x=m)) + geom_histogram(data= m, bins = 50, fill="#619CFF") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("NOME_94") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
n1<- ggplot(n, aes(x=n)) + geom_histogram(data= n, bins = 50, fill="#AE87FF") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SNOH_03") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
o1<- ggplot(o, aes(x=o)) + geom_histogram(data= o, bins = 50, fill="#DB72FB") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SNOH_96") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1<- ggplot(p, aes(x=p)) + geom_histogram(data= p, bins = 50, fill="#F564E3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SUSIT_14") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
q1<- ggplot(q, aes(x=q)) + geom_histogram(data= q, bins = 50, fill="#FF61C3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("TAUY_09") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
r1<- ggplot(r, aes(x=r)) + geom_histogram(data= r, bins = 50, fill="#FF699C") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("TAUY_12") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#list the plots and call them in their layout
grid.arrange(a1,b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1,nrow=4, top="Pre Filtered Data; Counts of Locus Genotype Rate by Population")


#### Post Filtered Populations 
Post_filtered_popz_LGR[1:5,1:5]

#Subset the matrix to get each population's numbers
a <-as.data.frame(Post_filtered_popz_LGR[1,] )
b <-as.data.frame(Post_filtered_popz_LGR[2,] )
c <-as.data.frame(Post_filtered_popz_LGR[3,] )
d <-as.data.frame(Post_filtered_popz_LGR[4,] )
e <-as.data.frame(Post_filtered_popz_LGR[5,] )
f <-as.data.frame(Post_filtered_popz_LGR[6,] )
g <-as.data.frame(Post_filtered_popz_LGR[7,] )
h <-as.data.frame(Post_filtered_popz_LGR[8,] )
i <-as.data.frame(Post_filtered_popz_LGR[9,] )
j <-as.data.frame(Post_filtered_popz_LGR[10,] )
k <-as.data.frame(Post_filtered_popz_LGR[11,] )
l <-as.data.frame(Post_filtered_popz_LGR[12,] )
m <-as.data.frame(Post_filtered_popz_LGR[13,] )
n <-as.data.frame(Post_filtered_popz_LGR[14,] )
o <-as.data.frame(Post_filtered_popz_LGR[15,] )
p <-as.data.frame(Post_filtered_popz_LGR[16,] )
q <-as.data.frame(Post_filtered_popz_LGR[17,] )
r <-as.data.frame(Post_filtered_popz_LGR[18,] )

#get the population names in order
row.names(Post_filtered_popz_LGR)

#assign each of the populations a ggplot to be called later
a1<- ggplot(a, aes(x=a)) + geom_histogram(data= a, bins = 50, fill="#F8766D") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("AMUR_10") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
b1<- ggplot(b, aes(x=b)) + geom_histogram(data= b, bins = 50, fill="#E88526") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("AMUR_11") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
c1<- ggplot(c, aes(x=c)) + geom_histogram(data= c, bins = 50, fill="#D39200") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SUSIT_13") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
d1<- ggplot(d, aes(x=d)) + geom_histogram(data= d, bins = 50, fill="#B79F00") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("HAYLY_09") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
e1<- ggplot(e, aes(x=e)) + geom_histogram(data= e, bins = 50, fill="#93AA00") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("HAYLY_10") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
f1<- ggplot(f, aes(x=f)) + geom_histogram(data= f, bins = 50, fill="#5EB300") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KOPE_91") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g1<- ggplot(g, aes(x=g)) + geom_histogram(data= g, bins = 50, fill="#00BA38") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KOPE_96") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
h1<- ggplot(h, aes(x=h)) + geom_histogram(data= h, bins = 50, fill="#00BF74") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KUSHI_06") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
i1<- ggplot(i, aes(x=i)) + geom_histogram(data= i, bins = 50, fill="#00C19F") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KUSHI_07") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
j1<- ggplot(j, aes(x=j)) + geom_histogram(data= j, bins = 50, fill="#00BFC4") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("LAKEL_06") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
k1<- ggplot(k, aes(x=k)) + geom_histogram(data= k, bins = 50, fill="#00B9C3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("LAKEL07") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
l1<- ggplot(l, aes(x=l)) + geom_histogram(data= l, bins = 50, fill="#00ADFA") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("NOME_91") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m1<- ggplot(m, aes(x=m)) + geom_histogram(data= m, bins = 50, fill="#619CFF") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("NOME_94") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
n1<- ggplot(n, aes(x=n)) + geom_histogram(data= n, bins = 50, fill="#AE87FF") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SNOH_03") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
o1<- ggplot(o, aes(x=o)) + geom_histogram(data= o, bins = 50, fill="#DB72FB") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SNOH_96") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1<- ggplot(p, aes(x=p)) + geom_histogram(data= p, bins = 50, fill="#F564E3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SUSIT_14") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
q1<- ggplot(q, aes(x=q)) + geom_histogram(data= q, bins = 50, fill="#FF61C3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("TAUY_09") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
r1<- ggplot(r, aes(x=r)) + geom_histogram(data= r, bins = 50, fill="#FF699C") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("TAUY_12") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#list the plots and call them in their layout
grid.arrange(a1,b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1,nrow=4, top="Post Filtered Data; Counts of Locus Genotype Rate by Population")



####################### Minor Allele Frequency MAF


#### OG data Populations 
OG_80PCfiltered_popz_MAF[1:5,1:5]

#Subset the matrix to get each population's numbers
a <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[1,]))
b <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[2,]))
c <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[3,]))
d <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[4,]))
e <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[5,]))
f <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[6,]))
g <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[7,]))
h <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[8,]))
i <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[9,]))
j <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[10,]))
k <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[11,]))
l <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[12,]))
m <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[13,]))
n <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[14,]))
o <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[15,]))
p <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[16,]))
q <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[17,]))
r <-as.data.frame(as.numeric(OG_80PCfiltered_popz_MAF[18,]))

#get the population names in order
row.names(OG_80PCfiltered_popz_MAF)

#assign each of the populations a ggplot to be called later
a1<- ggplot(a, aes(x=a)) + geom_histogram(data= a, bins = 50, fill="#F8766D") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("AMUR_10") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
b1<- ggplot(b, aes(x=b)) + geom_histogram(data= b, bins = 50, fill="#E88526") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("AMUR_11") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
c1<- ggplot(c, aes(x=c)) + geom_histogram(data= c, bins = 50, fill="#D39200") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SUSIT_13") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
d1<- ggplot(d, aes(x=d)) + geom_histogram(data= d, bins = 50, fill="#B79F00") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("HAYLY_09") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
e1<- ggplot(e, aes(x=e)) + geom_histogram(data= e, bins = 50, fill="#93AA00") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("HAYLY_10") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
f1<- ggplot(f, aes(x=f)) + geom_histogram(data= f, bins = 50, fill="#5EB300") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KOPE_91") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g1<- ggplot(g, aes(x=g)) + geom_histogram(data= g, bins = 50, fill="#00BA38") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KOPE_96") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
h1<- ggplot(h, aes(x=h)) + geom_histogram(data= h, bins = 50, fill="#00BF74") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KUSHI_06") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
i1<- ggplot(i, aes(x=i)) + geom_histogram(data= i, bins = 50, fill="#00C19F") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KUSHI_07") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
j1<- ggplot(j, aes(x=j)) + geom_histogram(data= j, bins = 50, fill="#00BFC4") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("LAKEL_06") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
k1<- ggplot(k, aes(x=k)) + geom_histogram(data= k, bins = 50, fill="#00B9C3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("LAKEL07") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
l1<- ggplot(l, aes(x=l)) + geom_histogram(data= l, bins = 50, fill="#00ADFA") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("NOME_91") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m1<- ggplot(m, aes(x=m)) + geom_histogram(data= m, bins = 50, fill="#619CFF") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("NOME_94") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
n1<- ggplot(n, aes(x=n)) + geom_histogram(data= n, bins = 50, fill="#AE87FF") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SNOH_03") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
o1<- ggplot(o, aes(x=o)) + geom_histogram(data= o, bins = 50, fill="#DB72FB") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SNOH_96") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1<- ggplot(p, aes(x=p)) + geom_histogram(data= p, bins = 50, fill="#F564E3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SUSIT_14") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
q1<- ggplot(q, aes(x=q)) + geom_histogram(data= q, bins = 50, fill="#FF61C3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("TAUY_09") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
r1<- ggplot(r, aes(x=r)) + geom_histogram(data= r, bins = 50, fill="#FF699C") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("TAUY_12") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#list the plots and call them in their layout
grid.arrange(a1,b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1,nrow=4, top="Pre Filtered Data; Counts of Minor Allele Frequencies by Population")


#### Post Filtered Populations 
Post_filtered_popz_MAF[1:5,1:5]

#Subset the matrix to get each population's numbers
a <-as.data.frame(as.numeric(Post_filtered_popz_MAF[1,]))
b <-as.data.frame(as.numeric(Post_filtered_popz_MAF[2,]))
c <-as.data.frame(as.numeric(Post_filtered_popz_MAF[3,]))
d <-as.data.frame(as.numeric(Post_filtered_popz_MAF[4,]))
e <-as.data.frame(as.numeric(Post_filtered_popz_MAF[5,]))
f <-as.data.frame(as.numeric(Post_filtered_popz_MAF[6,]))
g <-as.data.frame(as.numeric(Post_filtered_popz_MAF[7,]))
h <-as.data.frame(as.numeric(Post_filtered_popz_MAF[8,]))
i <-as.data.frame(as.numeric(Post_filtered_popz_MAF[9,]))
j <-as.data.frame(as.numeric(Post_filtered_popz_MAF[10,]))
k <-as.data.frame(as.numeric(Post_filtered_popz_MAF[11,]))
l <-as.data.frame(as.numeric(Post_filtered_popz_MAF[12,]))
m <-as.data.frame(as.numeric(Post_filtered_popz_MAF[13,]))
n <-as.data.frame(as.numeric(Post_filtered_popz_MAF[14,]))
o <-as.data.frame(as.numeric(Post_filtered_popz_MAF[15,]))
p <-as.data.frame(as.numeric(Post_filtered_popz_MAF[16,]))
q <-as.data.frame(as.numeric(Post_filtered_popz_MAF[17,]))
r <-as.data.frame(as.numeric(Post_filtered_popz_MAF[18,]))

#get the population names in order
row.names(Post_filtered_popz_MAF)

#assign each of the populations a ggplot to be called later
a1<- ggplot(a, aes(x=a)) + geom_histogram(data= a, bins = 50, fill="#F8766D") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("AMUR_10") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
b1<- ggplot(b, aes(x=b)) + geom_histogram(data= b, bins = 50, fill="#E88526") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("AMUR_11") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
c1<- ggplot(c, aes(x=c)) + geom_histogram(data= c, bins = 50, fill="#D39200") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SUSIT_13") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
d1<- ggplot(d, aes(x=d)) + geom_histogram(data= d, bins = 50, fill="#B79F00") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("HAYLY_09") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
e1<- ggplot(e, aes(x=e)) + geom_histogram(data= e, bins = 50, fill="#93AA00") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("HAYLY_10") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
f1<- ggplot(f, aes(x=f)) + geom_histogram(data= f, bins = 50, fill="#5EB300") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KOPE_91") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g1<- ggplot(g, aes(x=g)) + geom_histogram(data= g, bins = 50, fill="#00BA38") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KOPE_96") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
h1<- ggplot(h, aes(x=h)) + geom_histogram(data= h, bins = 50, fill="#00BF74") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KUSHI_06") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
i1<- ggplot(i, aes(x=i)) + geom_histogram(data= i, bins = 50, fill="#00C19F") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KUSHI_07") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
j1<- ggplot(j, aes(x=j)) + geom_histogram(data= j, bins = 50, fill="#00BFC4") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("LAKEL_06") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
k1<- ggplot(k, aes(x=k)) + geom_histogram(data= k, bins = 50, fill="#00B9C3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("LAKEL07") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
l1<- ggplot(l, aes(x=l)) + geom_histogram(data= l, bins = 50, fill="#00ADFA") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("NOME_91") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m1<- ggplot(m, aes(x=m)) + geom_histogram(data= m, bins = 50, fill="#619CFF") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("NOME_94") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
n1<- ggplot(n, aes(x=n)) + geom_histogram(data= n, bins = 50, fill="#AE87FF") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SNOH_03") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
o1<- ggplot(o, aes(x=o)) + geom_histogram(data= o, bins = 50, fill="#DB72FB") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SNOH_96") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1<- ggplot(p, aes(x=p)) + geom_histogram(data= p, bins = 50, fill="#F564E3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SUSIT_14") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
q1<- ggplot(q, aes(x=q)) + geom_histogram(data= q, bins = 50, fill="#FF61C3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("TAUY_09") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
r1<- ggplot(r, aes(x=r)) + geom_histogram(data= r, bins = 50, fill="#FF699C") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("TAUY_12") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#list the plots and call them in their layout
grid.arrange(a1,b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1,nrow=4, top="Post Filtered Data; Counts of Minor Allele Frequencies by Population")


####################### HWE 

#### Pre Filtered for HWE: 31485 Loci)
HWE_table_Pre_Filtered[1:5,1:5]

#Subset the matrix to get each population's numbers
a <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,1]))
b <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,2]))
c <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,3]))
d <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,4]))
e <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,5]))
f <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,6]))
g <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,7]))
h <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,8]))
i <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,9]))
j <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,10]))
k <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,11]))
l <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,12]))
m <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,13]))
n <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,14]))
o <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,15]))
p <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,16]))
q <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,17]))
r <-as.data.frame(as.numeric(HWE_table_Pre_Filtered[,18]))

#get the population names in order
colnames(HWE_table_Pre_Filtered)

#assign each of the populations a ggplot to be called later
a1<- ggplot(a, aes(x=a)) + geom_histogram(data= a, bins = 50, fill="#F8766D") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("AMUR_10") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
b1<- ggplot(b, aes(x=b)) + geom_histogram(data= b, bins = 50, fill="#E88526") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("AMUR_11") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
c1<- ggplot(c, aes(x=c)) + geom_histogram(data= c, bins = 50, fill="#D39200") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SUSIT_13") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
d1<- ggplot(d, aes(x=d)) + geom_histogram(data= d, bins = 50, fill="#B79F00") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("HAYLY_09") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
e1<- ggplot(e, aes(x=e)) + geom_histogram(data= e, bins = 50, fill="#93AA00") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("HAYLY_10") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
f1<- ggplot(f, aes(x=f)) + geom_histogram(data= f, bins = 50, fill="#5EB300") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KOPE_91") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g1<- ggplot(g, aes(x=g)) + geom_histogram(data= g, bins = 50, fill="#00BA38") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KOPE_96") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
h1<- ggplot(h, aes(x=h)) + geom_histogram(data= h, bins = 50, fill="#00BF74") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KUSHI_06") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
i1<- ggplot(i, aes(x=i)) + geom_histogram(data= i, bins = 50, fill="#00C19F") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KUSHI_07") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
j1<- ggplot(j, aes(x=j)) + geom_histogram(data= j, bins = 50, fill="#00BFC4") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("LAKEL_06") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
k1<- ggplot(k, aes(x=k)) + geom_histogram(data= k, bins = 50, fill="#00B9C3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("LAKEL07") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
l1<- ggplot(l, aes(x=l)) + geom_histogram(data= l, bins = 50, fill="#00ADFA") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("NOME_91") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m1<- ggplot(m, aes(x=m)) + geom_histogram(data= m, bins = 50, fill="#619CFF") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("NOME_94") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
n1<- ggplot(n, aes(x=n)) + geom_histogram(data= n, bins = 50, fill="#AE87FF") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SNOH_03") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
o1<- ggplot(o, aes(x=o)) + geom_histogram(data= o, bins = 50, fill="#DB72FB") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SNOH_96") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1<- ggplot(p, aes(x=p)) + geom_histogram(data= p, bins = 50, fill="#F564E3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SUSIT_14") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
q1<- ggplot(q, aes(x=q)) + geom_histogram(data= q, bins = 50, fill="#FF61C3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("TAUY_09") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
r1<- ggplot(r, aes(x=r)) + geom_histogram(data= r, bins = 50, fill="#FF699C") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("TAUY_12") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#+geom_point(data=x_ranked,aes(x=rank,y=GenoRate))

#list the plots and call them in their layout
grid.arrange(a1,b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1,nrow=4, top="Before Filtering for HWE; Counts of HWE p-vals by Population")


#### Post Filtered Populations 
HWE_table_Post_Filtered_test[1:5,1:5]

#Subset the matrix to get each population's numbers
a <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,1]))
b <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,2]))
c <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,3]))
d <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,4]))
e <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,5]))
f <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,6]))
g <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,7]))
h <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,8]))
i <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,9]))
j <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,10]))
k <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,11]))
l <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,12]))
m <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,13]))
n <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,14]))
o <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,15]))
p <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,16]))
q <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,17]))
r <-as.data.frame(as.numeric(HWE_table_Post_Filtered_test[,18]))

#get the population names in order
colnames(HWE_table_Post_Filtered_test)

#assign each of the populations a ggplot to be called later
a1<- ggplot(a, aes(x=a)) + geom_histogram(data= a, bins = 50, fill="#F8766D") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("AMUR_10") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
b1<- ggplot(b, aes(x=b)) + geom_histogram(data= b, bins = 50, fill="#E88526") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("AMUR_11") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
c1<- ggplot(c, aes(x=c)) + geom_histogram(data= c, bins = 50, fill="#D39200") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SUSIT_13") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
d1<- ggplot(d, aes(x=d)) + geom_histogram(data= d, bins = 50, fill="#B79F00") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("HAYLY_09") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
e1<- ggplot(e, aes(x=e)) + geom_histogram(data= e, bins = 50, fill="#93AA00") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("HAYLY_10") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
f1<- ggplot(f, aes(x=f)) + geom_histogram(data= f, bins = 50, fill="#5EB300") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KOPE_91") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
g1<- ggplot(g, aes(x=g)) + geom_histogram(data= g, bins = 50, fill="#00BA38") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KOPE_96") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
h1<- ggplot(h, aes(x=h)) + geom_histogram(data= h, bins = 50, fill="#00BF74") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KUSHI_06") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
i1<- ggplot(i, aes(x=i)) + geom_histogram(data= i, bins = 50, fill="#00C19F") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("KUSHI_07") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
j1<- ggplot(j, aes(x=j)) + geom_histogram(data= j, bins = 50, fill="#00BFC4") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("LAKEL_06") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
k1<- ggplot(k, aes(x=k)) + geom_histogram(data= k, bins = 50, fill="#00B9C3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("LAKEL07") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
l1<- ggplot(l, aes(x=l)) + geom_histogram(data= l, bins = 50, fill="#00ADFA") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("NOME_91") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
m1<- ggplot(m, aes(x=m)) + geom_histogram(data= m, bins = 50, fill="#619CFF") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("NOME_94") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
n1<- ggplot(n, aes(x=n)) + geom_histogram(data= n, bins = 50, fill="#AE87FF") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SNOH_03") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
o1<- ggplot(o, aes(x=o)) + geom_histogram(data= o, bins = 50, fill="#DB72FB") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SNOH_96") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1<- ggplot(p, aes(x=p)) + geom_histogram(data= p, bins = 50, fill="#F564E3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("SUSIT_14") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
q1<- ggplot(q, aes(x=q)) + geom_histogram(data= q, bins = 50, fill="#FF61C3") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("TAUY_09") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
r1<- ggplot(r, aes(x=r)) + geom_histogram(data= r, bins = 50, fill="#FF699C") + labs(x= NULL, y=NULL) + theme_bw() + guides(fill= FALSE) + ggtitle("TAUY_12") + theme(plot.title = element_text(size=11), panel.grid.major = element_blank(), panel.grid.minor = element_blank())

#list the plots and call them in their layout
grid.arrange(a1,b1, c1, d1, e1, f1, g1, h1, i1, j1, k1, l1, m1, n1, o1, p1, q1, r1,nrow=4, top="After Filtering for HWE; Counts of HWE p-vals by Population")


# this code is from Eleni, she used it to plot something similar in her data. 
#This code requires a dataframe with the individuals and a column that has population designation
#ggplot(mydata, aes(x=loci, fill= population)) + geom_histogram(data= mydata, bins = 20) + facet_wrap(~population) 

####################### Individual Heterozygosity by Population

# plot Counts of Individual Heterozygosity by Population
ggplot(HetCounts_popz, aes(x=Percent_Het, fill= Pop)) + 
  geom_histogram(data=HetCounts_popz, bins = 20) + facet_wrap(~Pop) + theme_bw()  + guides(fill= FALSE) +
  ggtitle("Counts of Individual Percent Heterozygosity by Population")

#Mean Individual Heterozygosity by Population

# plot Mean Individual Heterozygosity by Population
ggplot(HetCounts_popz, aes(x=factor(Pop), y=Percent_Het, fill=HetCounts_popz$Pop)) + stat_summary(fun.y= "mean", geom="bar") + theme_bw() + 
  ggtitle("Mean Individual Percent Heterozygosity by Population, Post Filtering") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Population") + guides(fill= FALSE)



####################### Number of individuals by population 
popmap465_NEW<- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/PopMap465_NEWNAMES_NEWSAMPLES.txt")
popmap492<- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/PopMap_NEWNAMES.txt")

head(popmap465_NEW)
head(popmap492)

colnames(popmap465_NEW) <- c("Sample","Pop")
colnames(popmap492) <- c("Sample","Pop")

# Post Filtering
ggplot(data=popmap465_NEW) + geom_bar(aes(x=popmap465_NEW$Pop, fill=popmap465_NEW$Pop), stat= "count") + theme_bw() + 
  ggtitle("Number of Samples per Population, After Filtering") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Population") + guides(fill= FALSE)


#Pre-Filtering
ggplot(data=popmap492) + geom_bar(aes(x=popmap492$Pop, fill=popmap492$Pop), stat= "count") + theme_bw() + 
  ggtitle("Number of Samples per Population, Before Filtering") +  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  xlab("Population") + guides(fill= FALSE)



######################################################
######################################################
###################################################### PCA Colors 


# ALL_col <-c(AMUR10 = "#cd5490", AMUR11 = "#cd5490", SUSIT13= "#563e7f", HAYLY09 = "#c9cf4a", HAYLY10 = "#c9cf4a", 
#           KOPPE96 = "#7297ee", KOPPE91 = "#7297ee", KUSHI06 = "#61005e", KUSHI07 = "#61005e", LAKEL06 = "#2171b5",LAKEL07 = "#2171b5",
#            NOME91 = "#677f3e", NOME94 = "#677f3e", SNOH03 = "#00158a", SNOH96 ="#00158a", SUSIT14 ="#563e7f",TAUY09 = "#e0957e", TAUY12 = "#e0957e")

ALL_col_PCA <-c("#cd5490","#cd5490", "#563e7f", "#c9cf4a", "#c9cf4a", 
                "#7297ee", "#7297ee", "#61005e","#61005e","#2171b5","#2171b5",
                "#677f3e", "#677f3e", "#00158a", "#00158a","#563e7f","#e0957e", "#e0957e")

# ALL_col<- c("#cd5490","#cd5490","#c9cf4a","#c9cf4a","#7297ee","#7297ee","#61005e","#61005e", "#677f3e","#677f3e","#00158a","#00158a","#e0957e","#e0957e","#2171b5","#2171b5","#563e7f", "#563e7f")
# ALL_col_l <- rep(1, length=length(ALL_col)) 
# names(ALL_col_l) <- ALL_col
# pie(ALL_col_l, col=ALL_col, cex=.75, main = "ALL_COL")
