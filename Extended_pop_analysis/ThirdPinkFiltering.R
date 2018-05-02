### Filter One SNP per tag dataset for paralogs
###   Caculate the genotype rate of the Haplotype data set 
###   Remove the paralogs from Haplotype data set for Random Forest   
###   Carolyn Tarpey | December 2017
### ---------------------------------------

#This code requires the lists of paralogs from HDPLOT_forPinks.R
#The one SNP per tag genotypes should be filtered through SecondPinkFiltering.R already
#The haplotype file is from STACKS using the 31485 whitelist after initialPinkFiltering.R

#install.packages("vcfR")

#Pink data filtering
library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)

#input the haplotype file, make sure all the headers have no spaces
Haplotype_file<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/Whitelist31485_populations/batch_4.haplotypes.tsv", header=TRUE)
Haplotype_file[1:5,1:5]

#input the single snp filtered genotype file, make sure all the headers have no spaces
genotype_file<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/FilteredGenos_FilteredInds.txt", header=TRUE, colClasses="factor")
genotype_file[1:5,1:5]
dim(genotype_file)

#load the list of loci that were any paralog from HDPLOT_forPinks.R  
paralogs_alldata<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/HDPlot/any_paralog.txt", colClasses="factor")
colnames(paralogs_alldata) <- c("Position", "Tag", "Locus", "Identity", "Identity_2")
head(paralogs_alldata)

#load the list of loci that are singletons from HDPLOT_forPinks.R  
singletons<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/HDPlot/singletons.txt", colClasses="factor")
colnames(singletons) <- c("Position", "Tag", "Locus", "Identity", "Identity_2")
head(singletons)


#load the list of individuals that passed the 80% genotype rate from SecondPinkFiltering.R 
indvs_list<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/SecondFiltering/filteredSample_names.txt", colClasses="character")
#indvs_list <- gsub(",","",paste(indvs_list)) #the names are factors so you can't gsub it w/o paste
#indvs_list <- as.data.frame(indvs_list)
colnames(indvs_list) <-"Sample"
head(indvs_list)
dim(indvs_list)

###################### Get lists of the paralog statuses

#get the unique tags for the paralogs  
paralog_tags<- unique(paralogs_alldata$Tag)
length(paralog_tags)
head(paralog_tags)

#get the unique tags for the singletons  
singletons_tags<- unique(singletons$Tag)
length(singletons_tags)
head(singletons_tags)

#get the tags for singletons that are not in the paralog list
just_singletons_tags<- setdiff(singletons_tags, paralog_tags)
length(just_singletons_tags)
#xx <- which(just_singletons_tags == "153483")
#just_singletons_tags[xx] 

######################## HAPLOTYPE FILTER THE INDIVIDUALS BASED ON 80% genotype rate from SecondPinkFiltering.R

haplo_indvfil <- Haplotype_file[,colnames(Haplotype_file)%in%indvs_list$Sample]
dim(haplo_indvfil)
haplo_indvfil[1:5,1:5]
#Haplotype_file[1:5,1:5]

#add the tag number back to the beginning of the indv filtered haplo file
Catalog_Ids <- Haplotype_file$Catalog_ID
haplo_indvfil <- cbind(Catalog_Ids, haplo_indvfil)
haplo_indvfil[1:5,1:5]


######################## HAPLOTYPE FILTER THE GENOTYPES BASED ON PARALOG RESULTS

#filter to keep only singleton haplotypes
haplo_indv_parafil <- haplo_indvfil[haplo_indvfil$Catalog_Ids%in%just_singletons_tags, ]
haplo_indv_parafil[1:5,1:5]
dim(haplo_indv_parafil)

######################## HAPLOTYPE calculate and plot idv and loci genotype rate to see if there are any more that need to be removed 

#Do the calculation for the individual genotype rate 

#get genotype rate per individual
IND_genorate<-apply(haplo_indv_parafil,2,function(x) 1-(sum(x=="-")/length(haplo_indv_parafil$Catalog_Ids)[1]))
IND_genorate<-data.frame(keyName=names(IND_genorate), value=IND_genorate, row.names=NULL)
IND_genorate <- IND_genorate[-1,] #for the catalogID 
colnames(IND_genorate)<-c("Sample","GenoRate")
dim(IND_genorate)
head(IND_genorate)
min(IND_genorate[,2])

#plot ranked genotype rate for samples with filtered loci
inds_ranked <- IND_genorate[order(IND_genorate$GenoRate),]
inds_ranked$rank<-seq(1,dim(inds_ranked)[1],by=1)
ggplot()+geom_point(data=inds_ranked,aes(x=rank,y=GenoRate))+ggtitle("Individual Genotype Rate for filtered Haplotypes") +theme_bw()
ggsave(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/HDPlot/INDGenotypeRateForHaplotypesFiltered.pdf")


#Do the calculation for the haplotype genotype rate 
#get genotype rate per locus
CAT_IDS <-haplo_indv_parafil$Catalog_Ids
#head(CAT_IDS)
locus_genorate<-apply(haplo_indv_parafil[,-1],1,function(x) 1-(sum(x=="-")/(ncol(haplo_indv_parafil)-1))) #remove one from row count for the CatalogID
locus_genorate<-data.frame(keyName=CAT_IDS, value=locus_genorate)
colnames(locus_genorate)<-c("Locus","GenoRate")
head(locus_genorate)
dim(locus_genorate)
min(locus_genorate[,2])

#plot ranked genotype rate for samples with filtered loci
Loci_ranked <- locus_genorate[order(locus_genorate$GenoRate),]
Loci_ranked$rank<-seq(1,dim(locus_genorate)[1],by=1)
ggplot()+geom_point(data=Loci_ranked,aes(x=rank,y=GenoRate)) +ggtitle("Locus Genotype Rate for filtered Haplotypes")+theme_bw()
ggsave(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/HDPlot/LOCUSGenotypeRateForHaplotypesFiltered.pdf")


###Don't have to do any additional Genotype or Individual filtering
###### export the haplotype file that has been filtered for the individuals and paralogs removed

outputFile<-file("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/HDPlot/Haplotypes_with465ins_singletons.txt", "wb")
write.table(haplo_indv_parafil,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)



##### GENOTYPE filter the loci for paralogs
######################## FILTER THE GENOTYPES BASED ON PARALOG RESULTS
#format dataset to remove X from locus names
genotype_file_t<-gsub("X","",colnames(genotype_file))
head(genotype_file_t)
length(genotype_file_t)
colnames(genotype_file) <- genotype_file_t
genotype_file[1:5,1:5]
dim(genotype_file)

#pull the loci names from the genotype file 
genotype_loci <- colnames(genotype_file)
head(genotype_loci)

#break up the loci names into their components in a table 
loci_table_t <- as.data.frame(str_split_fixed(genotype_loci, "_", 2))
ncol(loci_table_t)
colnames(loci_table_t) <- c("Tag", "Pos")
loci_table_t$Locus <- genotype_loci
head(loci_table_t)
dim(loci_table_t)

#get the tags for singletons that are not in the paralog list
just_singletons_tags <- as.numeric(setdiff(singletons_tags, paralog_tags))
length(just_singletons_tags)
head(just_singletons_tags)

xx <- which(loci_table_t$Tag == "44444")
loci_table_t$Tag[xx]

xx <- which(loci_table_t$Tag == "44444")
loci_table_t$Tag[xx]

singletons_to_keep <- loci_table_t[which(loci_table_t$Tag%in%just_singletons_tags),]
dim(singletons_to_keep)
head(singletons_to_keep)

singletons_locus_keep<- singletons_to_keep$Locus
length(singletons_locus_keep)
head(singletons_locus_keep)

filtered_singleton_genotypes <- genotype_file[,colnames(genotype_file)%in%singletons_locus_keep]
dim(filtered_singleton_genotypes)
filtered_singleton_genotypes[1:5,1:5]

####Write a table of the genotypes/individuals that have passed all the filters so far:
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/ThirdFiltering/Singleton_One_Tag_Genotypes.txt", "wb")
write.table(filtered_singleton_genotypes,outputFile,quote=FALSE,row.names=TRUE,col.names=TRUE,eol="\n")
close(outputFile)

####Write a list of the SNPSthat have passed all the filters so far:
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/ThirdFiltering/Singleton_One_Tag_SNPLIST.txt", "wb")
write.table(singletons_locus_keep,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

#####################Test to see that what we got in the final paralog filtered one snp per tag data set has the 16681 loci

#if we want to see which of our 16681 snps is included at this point
loci_16681 <- readLines("Z:/WORK/TARPEY/Pink_Populations/listof16681LOCI.txt")
head(loci_16681)

#list of all the loci that were singletons that we kept
length(singletons_locus_keep)

#loci that are in  the singleton loci and the remaing 14629 loci of the 16681
length(in_14637_and_30088)
in_14637_and_23759 <- intersect(in_14637_and_30088, singletons_locus_keep)
length(in_14637_and_23759)

diff_14637_and_23759 <- setdiff(in_14637_and_30088, singletons_locus_keep)
length(diff_14637_and_23759)

#####################Test to see that what we got in the final paralog filtered Haplotype data set has the 16681 loci

#if we want to see which of our 16681 snps is included at this point
loci_16681 <- readLines("Z:/WORK/TARPEY/Pink_Populations/listof16681LOCI.txt")
loci_16681_tags <- as.data.frame(str_split_fixed(loci_16681, "_", 2))
head(loci_16681_tags)
colnames(loci_16681_tags)<- c("Tag","Pos")
Just_loci_16681_tags<-loci_16681_tags$Tag 
length(Just_loci_16681_tags)

#list of all the loci that were singletons that we kept
length(just_singletons_tags)

#loci that are in  the singleton loci and the remaing 14629 loci of the 16681

HaplotypeTags_16681Tags <- intersect(just_singletons_tags, Just_loci_16681_tags)
length(HaplotypeTags_16681Tags)

diff_14637_and_23759 <- setdiff(in_14637_and_30088, singletons_locus_keep)
length(diff_14637_and_23759)

###########################################


