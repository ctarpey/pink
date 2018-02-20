### To select markers for two pink Amplicon Panels, one per lineage
###   Using Ranked FST of one snp per tag loci in a hierarchical structure
###   Remove the paralogs from Haplotype data set for Random Forest   
###   Carolyn Tarpey | January 2018
### ---------------------------------------

#This code requires the FST from Genepop. 
#There are three levels to the heirarchy: 
#1. All populations, and All  populations except Susitna 
#2. All Even Populations except Susitna, and all Odd populations except Susitna
#3. Even NA pops w/o Susitna, Odd NA pops w/o Susitna, Even Asia pops, and Odd Asia pops
#Other: Additional level is all the NA pops w/o Susitna and all the Asian pops, 
  #but these are not really for answering out specific marker selection question, just for general interest

#The one SNP per tag genotypes should be filtered through SecondPinkFiltering.R already
#The haplotype file is from STACKS using the 31485 whitelist after initialPinkFiltering.R

#install.packages("vcfR")

#Pink data filtering
library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)
library(hierfstat)
library(adegenet)
library(argparse)
library(stringi)
#install.packages("argparse")


#input the text file that has the results for all the runs of global FST from Genepop
#FST was run on each of the groupings in the three hierarchies listed above, for a total of 10 columns of FST results

FST_file<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Genepop/HWE_Genepop_FST_R.txt", header=TRUE, 
                     stringsAsFactors = FALSE, na.strings = "-" )
FST_file[1:5,1:5]
dim(FST_file)

#input the text file that the modified SNPs from the Catalog file. It has been edited in Excel to have a list of the 
#locus names, as well as the SNP index for each snp at each tag. These are all the SNPS in the whole catalog. 

SNP_positions_all <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/SNPposition/SNP_positions.txt", header=TRUE, 
                                stringsAsFactors = FALSE, na.strings = "-" )
head(SNP_positions_all)
dim(SNP_positions_all)


###################################### SNP Position

##import the haplotypes into R to look for fixed snps at each tag in the lineages
Haplotypes<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/SNPposition/Haplotypes_with465ins_singletons.txt", header=TRUE, 
                          stringsAsFactors = FALSE, na.strings = "-" )
Haplotypes[1:5,1:5]
dim(Haplotypes)


#split the haplotype file by lineage, and remove the Susitna pops by ignoring them
#import the population map that has pop and lineage for each sample. it has "ignore" the susitna lineages
#the lalkelse pops are named their correct lineage, though it doesn't match their population label

popMapLineages<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/SNPposition/PopMapLineages.txt", header=FALSE, 
                          stringsAsFactors = FALSE, na.strings = "-" )
colnames(popMapLineages)<- c("Sample","Pop","Lineage")
head(popMapLineages)
dim(popMapLineages)


##make a dataframe for each lineage with the individuals in that lineages
Even_inds<-popMapLineages[popMapLineages$Lineage=="Even",]
Odd_inds<-popMapLineages[popMapLineages$Lineage=="Odd",]

head(Even_inds)
dim(Even_inds)

dim(Odd_inds)

##filter the haplotype file by lineage

Even_filteredHaplotypes <-Haplotypes[,colnames(Haplotypes)%in%Even_inds$Sample]
Even_filteredHaplotypes[1:5,1:5]
dim(Even_filteredHaplotypes)
names(Even_filteredHaplotypes)

Odd_filteredHaplotypes <-Haplotypes[,colnames(Haplotypes)%in%Odd_inds$Sample]
Odd_filteredHaplotypes[1:5,1:5]
dim(Odd_filteredHaplotypes)
names(Odd_filteredHaplotypes)

# filtered_singleton_genotypes <- genotype_file[,colnames(genotype_file)%in%singletons_locus_keep]
# dim(filtered_singleton_genotypes)
# filtered_singleton_genotypes[1:5,1:5]

is.data.frame(Even_filteredHaplotypes)

#####  THESE WORK BUT TAKE FOREVER:: keep them commented out until you need them again
# First sweep through each list of filtered Haplotypes to remove the genotypes at any individuals with > 2 alleles
# 
#  for (i in 1:nrow(Even_filteredHaplotypes)) {
#    row <- as.character(Even_filteredHaplotypes[i,])
#    for (k in 1:length(row)) {
#      if (!is.na(row[k])) {
#        if (str_count(row[k], pattern = "/") > 1) {
#          Even_filteredHaplotypes[i,k] <- NA
#        } else {
#          Even_filteredHaplotypes[i,k] <- Even_filteredHaplotypes[i,k]
#        }
#      }
#    }
#  }
# 
# 
#  #Even_filteredHaplotypes[92,]
# 
#  for (i in 1:nrow(Odd_filteredHaplotypes)) {
#    row <- as.character(Odd_filteredHaplotypes[i,])
#    for (k in 1:length(row)) {
#      if (!is.na(row[k])) {
#        if (str_count(row[k], pattern = "/") > 1) {
#          Odd_filteredHaplotypes[i,k] <- NA
#        } else {
#          Odd_filteredHaplotypes[i,k] <- Odd_filteredHaplotypes[i,k]
#        }
#      }
#    }
#  }
# 
# 
Even_filteredHaplotypes[1:5,1:5]


# i <- 92
# row <- as.character(test_haplo[154,])
#test_haplo <- Even_filteredHaplotypes[c(153,92),]
# test_haplo <- Even_filteredHaplotypes[1:200,]
# test_haplo[1:5,1:5]


haplo <- Even_filteredHaplotypes
#haplo <- Odd_filteredHaplotypes

base_set <- c("A", "C", "G", "t")
KEEP <- vector("list",100000)
THROW_OUT <- vector("list",100000 )
j = 1 #index for adding things to the list of YEs or NO

#for each row in the datafile 
for (i in 1:nrow(haplo)) { 
  
  row <- as.character(haplo[i,])
  #print(row)
  row_number <- i #this is the row index, not the tag
  tag <- Haplotypes[i,"Catalog_Ids"]
  #geno_count <- sum(!is.na(row)) #Count the number of genotypes at that tag
  #print(geno_count)
  min_allele <- min(nchar(row, type = "chars")) #this doesn't ignore NA, codes them as 2, but we're only looking for 1's
  #print(min_allele)
  #max_allele <- max(nchar(row, type = "chars")) # we dont use this, but it is interesting- here it includes characters in all alleles
  #print(max_allele)
  
  if (min_allele == 1) {  #if single snp per tag, (only one base in the list)
    hets <- sum(grepl("/", row)) # count the number of / in the row
    #are all the letters the same? if so, fixed. 
    if (hets >= 2) { # if there are more than 2 /, it is het
      KEEP[[j]] <- paste(tag, "1", sep="_") #add it to the Keep list 
      j = j + 1 #this is the index for adding things to the list, tags/ position would overlap
      next()
    } else { 
      counts <- table(na.omit(row))
      if (length(counts) == 1 ) {
        THROW_OUT[[j]] <- paste(tag, "1", sep="_") #add it to the throw out list 
        j = j + 1 #this is the index for adding things to the list, tags/ position would overlap  
        next()
      } else if (length(counts)== 2 & any(counts < 2)) { 
        THROW_OUT[[j]] <- paste(tag, "1", sep="_") #add it to the throw out list 
        j = j + 1 #this is the index for adding things to the list, tags/ position would overlap  
        next()
      } else if (length(counts)== 3 & ((counts[1]<2) + (counts[2]<2) + (counts[3]<2)) <= 2) { #boolean logic to choose if two of these are < 2
        THROW_OUT[[j]] <- paste(tag, "1", sep="_") #add it to the throw out list 
        j = j + 1 #this is the index for adding things to the list, tags/ position would overlap  
        next()
      } else if (length(counts)== 4 & ((counts[1]<2) + (counts[2]<2) + (counts[3]<2)+ (counts[4]<2)) <= 3) { #boolean logic to choose if two of these are < 2
        THROW_OUT[[j]] <- paste(tag, "1", sep="_") #add it to the throw out list 
        j = j + 1 #this is the index for adding things to the list, tags/ position would overlap  
        next()
      } else {
        KEEP[[j]] <- paste(tag, "1", sep="_") #add it to the Keep list 
        j = j + 1 #this is the index for adding things to the list, tags/ position would overlap
        next()
      }
    }
  } else if (min_allele != 1) {
    row <- str_split_fixed(row,"/",Inf)
    a_test <- table(row)
      if (length(a_test) == 1){
        for (p in 1:nchar(names(a_test[1]))){ ###<- For tags that have fixed homo variation, names(table[1]) returns the name, not ""
          THROW_OUT[[j]] <- paste(tag, p, sep="_") #add it to the throw out list 
          j = j + 1 #this is the index for adding things to the list, tags/ position would overlap
          #print(paste(tag, p, sep="_"))
        }
      } else{
          alleles <- names(a_test[2:length(a_test)])
          SNPs <- nchar(alleles[1]) # this counts the number of positions in the first allele name, which should be standard for all
          for(m in 1:(SNPs)){
            A <- 0
            C <- 0
            G <- 0
            t <- 0
            for (l in 1:length(alleles)){
              var <- substr(alleles[l], m, m)
              if (var == "A"){
                A <- A+1 
              } else if (var == "C"){
                C <- C+1 
              } else if (var == "G"){
                G <- G+1 
              } else if (var == "T"){
                t <- t+1 }
            }
            if (length(which(c(A, C, G, t) == 0)) == 3){
              THROW_OUT[[j]] <- paste(tag, m, sep="_") #add it to the throw out list 
              j = j + 1 #this is the index for adding things to the list, tags/ position would overlap 
              #print(paste(tag, m, sep="_"))  This works
            } else if (length(which(c(A, C, G, t) >= 2)) >= 2){
              KEEP[[j]] <- paste(tag, m, sep="_") #add it to the Keep list 
              j = j + 1 #this is the index for adding things to the list, tags/ position would overlap 
              #print(paste(tag, m, sep="_")) This works
            } else if ((length(which(c(A, C, G, t) >= 2)) == 1) & (length(which(c(A, C, G, t) == 1)) >= 1)){
                if (length(which(c(A, C, G, t) == 1)) == 1){
                  base <- base_set[which(c(A, C, G, t) == 1)] #gives the base we need to search for ie "A"
                  for (n in 1:length(alleles)){
                    if ((substr(alleles[n], m, m) == base) & (as.vector(a_test[n+1]) >=2)) {
                      KEEP[[j]] <- paste(tag, m, sep="_") #add it to the Keep list 
                      j = j + 1 #this is the index for adding things to the list, tags/ position would overlap
                  } else {
                    THROW_OUT[[j]] <- paste(tag, m, sep="_") #add it to the throw out list 
                    j = j + 1 #this is the index for adding things to the list, tags/ position would overlap
                  }
                  } 
                } else if (length(which(c(A, C, G, t) == 1)) == 2){
                      base <- base_set[which(c(A, C, G, t) == 1)] #gives the base we need to search for ie "A"
                      for (n in 1:length(alleles)){
                        if ((substr(alleles[n], m, m) == base[1]) & (as.vector(a_test[n+1]) >=2)) {
                          KEEP[[j]] <- paste(tag, m, sep="_") #add it to the Keep list 
                          j = j + 1 #this is the index for adding things to the list, tags/ position would overlap
                        } else if ((substr(alleles[n], m, m) == base[2]) & (as.vector(a_test[n+1]) >=2)) {
                          KEEP[[j]] <- paste(tag, m, sep="_") #add it to the Keep list 
                          j = j + 1 #this is the index for adding things to the list, tags/ position would overlap
                      } else {
                        THROW_OUT[[j]] <- paste(tag, m, sep="_") #add it to the throw out list 
                        j = j + 1 #this is the index for adding things to the list, tags/ position would overlap
                      }
                      }
                  } else if ((length(which(c(A, C, G, t) >= 2)) == 1) & (length(which(c(A, C, G, t) == 1)) == 3)){
                      base <- base_set[which(c(A, C, G, t) == 1)] #gives the base we need to search for ie "A"
                      for (n in 1:length(alleles)){
                        if ((substr(alleles[n], m, m) == base[1]) & (as.vector(a_test[n+1]) >=2)) {
                          KEEP[[j]] <- paste(tag, m, sep="_") #add it to the Keep list 
                          j = j + 1 #this is the index for adding things to the list, tags/ position would overlap
                        } else if ((substr(alleles[n], m, m) == base[2]) & (as.vector(a_test[n+1]) >=2)) {
                          KEEP[[j]] <- paste(tag, m, sep="_") #add it to the Keep list 
                          j = j + 1 #this is the index for adding things to the list, tags/ position would overlap
                        } else if ((substr(alleles[n], m, m) == base[3]) & (as.vector(a_test[n+1]) >=2)) {
                          KEEP[[j]] <- paste(tag, m, sep="_") #add it to the Keep list 
                          j = j + 1 #this is the index for adding things to the list, tags/ position would overlap
                        } else {
                          THROW_OUT[[j]] <- paste(tag, m, sep="_") #add it to the throw out list 
                          j = j + 1 #this is the index for adding things to the list, tags/ position would overlap
                        }
                      }
                    }
                }
            }
            }
          }
    }

KEEP[1:50] 
THROW_OUT[1:50]



### EVEN ~~~~~~~~~~~~~~~~~~take the Keep and throw out lists and remove most of the NULLs 
# 
# KEEP_even_raw <- KEEP
# THROW_OUT_even_raw <- THROW_OUT
# 
# KEEP_Even_filt_Hap <- KEEP[-which(sapply(KEEP, is.null))] #remove the NULLS
# length(KEEP_Even_filt_Hap)
# 
# INT_THROW_Even_filt_Hap <- THROW_OUT[-which(sapply(THROW_OUT, is.null))] # remove the NULLS
# length(INT_THROW_Even_filt_Hap)
# 
# INT_THROW_Even_filt_Hap <- as.vector(unique(INT_THROW_Even_filt_Hap))# remove the duplicates
# length(INT_THROW_Even_filt_Hap)
# 
# THROW_Even_filt_Hap <- setdiff(INT_THROW_Even_filt_Hap,KEEP_Even_filt_Hap) #in INT_throw that were not in KEEP
# length(THROW_Even_filt_Hap)


# ### ODD ~~~~~~~~~~~~~~~~~~take the Keep and throw out lists and  NULLs 
# # KEEP_odd_raw
# KEEP_odd_raw <- KEEP
# THROW_OUT_odd_raw <- THROW_OUT
# 
# KEEP_ODD_filt_Hap <- KEEP[-which(sapply(KEEP, is.null))] #remove the NULLS
# length(KEEP_ODD_filt_Hap)
# 
# INT_THROW_ODD_filt_Hap <- THROW_OUT[-which(sapply(THROW_OUT, is.null))] # remove the NULLS
# length(INT_THROW_ODD_filt_Hap)
# 
# INT_THROW_ODD_filt_Hap <- as.vector(unique(INT_THROW_ODD_filt_Hap))# remove the duplicates
# length(INT_THROW_ODD_filt_Hap)
# 
# THROW_ODD_filt_Hap <- setdiff(INT_THROW_ODD_filt_Hap,KEEP_ODD_filt_Hap) #in INT_throw that were not in KEEP
# length(THROW_ODD_filt_Hap)

# # KEEP_odd_raw
# KEEP_odd_raw <- KEEP
# THROW_OUT_odd_raw <- THROW_OUT
# 
# KEEP_ODD_filt_Hap <- KEEP_odd_raw[-which(sapply(KEEP_odd_raw, is.null))] #remove the NULLS
# length(KEEP_ODD_filt_Hap)
# 
# INT_THROW_ODD_filt_Hap <- THROW_OUT_odd_raw[-which(sapply(THROW_OUT_odd_raw, is.null))] # remove the NULLS
# length(INT_THROW_ODD_filt_Hap)
# 
# INT_THROW_ODD_filt_Hap <- as.vector(unique(INT_THROW_ODD_filt_Hap))# remove the duplicates
# length(INT_THROW_ODD_filt_Hap)
# 
# THROW_ODD_filt_Hap <- setdiff(INT_THROW_ODD_filt_Hap,KEEP_ODD_filt_Hap) #in INT_throw that were not in KEEP
# length(THROW_ODD_filt_Hap)




## convert the vector of filtered real SNP indexes to a Data Frame
Keep_even <-KEEP_Even_filt_Hap
Even_filt_snps <-data.frame(stri_list2matrix(Keep_even, byrow=TRUE))
colnames(Even_filt_snps)<-"Even_SNP_Index"
head(Even_filt_snps)

Keep_odd <- KEEP_ODD_filt_Hap
Odd_filt_snps <-data.frame(stri_list2matrix(Keep_odd, byrow=TRUE))
colnames(Odd_filt_snps)<-"Odd_SNP_Index"
head(Odd_filt_snps)



### Find the true SNP position for the variation at each tag. 

head(SNP_positions_all)
head(Even_filt_snps)
head(Odd_filt_snps)

Even_SNPS <- SNP_positions_all[SNP_positions_all$SNP_tag_index %in% Even_filt_snps$Even_SNP_Index,]
Odd_SNPS <- SNP_positions_all[SNP_positions_all$SNP_tag_index %in% Odd_filt_snps$Odd_SNP_Index,]

head(Even_SNPS)
head(Odd_SNPS)

########################### More filtering of SNPS based on position in the tag. 

###Flag SNP Positions >= 16 and <= 74, (we want tags 17-73 to be able to make primers)

dim(Even_SNPS)

#add new empty column to the end of each data frame for the flags
Even_SNPS$Btw_17_73 <- NA
Odd_SNPS$Btw_17_73 <- NA

#make the position numeric so you can compare with numbers
Even_SNPS$Position <- as.numeric(Even_SNPS$Position)
Odd_SNPS$Position <- as.numeric(Odd_SNPS$Position)

for (r in 1:dim(Even_SNPS)[1]){
  if ((Even_SNPS[r,3] >= 17) & (Even_SNPS[r,3] <= 73)) {
    Even_SNPS[r,10] <- "TRUE"
  } else {
    Even_SNPS[r,10] <- "FALSE"
  }
}

for (r in 1:dim(Odd_SNPS)[1]){
  if ((Odd_SNPS[r,3] >= 17) & (Odd_SNPS[r,3] <= 73)) {
    Odd_SNPS[r,10] <- "TRUE"
  } else {
    Odd_SNPS[r,10] <- "FALSE"
  }
}

head(Even_SNPS)
head(Odd_SNPS)

###Flag TAGS that have a snp with Positions >= 16 and <= 74- we don't want them

#add new empty column to the end of each data frame for the flags
head(test)

test<- Even_SNPS
test_tags<- unique(Even_SNPS$Tag)
test$SNPs_InRange <- NA


for (s in 1:dim(test)[1]){
  for (u in 1:length(test_tags)) {
    if (any(test[]))  ####<------------------------------------------------START HERE, testing the positions in each tag to make sure none have FALSE in btw
  }
}
  
  
#Even_SNPS$SNPs_InRange <- NA
#Odd_SNPS$SNPs_InRange <- NA







###################################### Marker selection based on ranked FSt of loci 


##############EVEN LINEAGE PANEL 
names(FST_file)

##Find the top FST for the NA_even
FST_NA_even<- FST_file[,c("Locus","NA_even")]
head(FST_NA_even)

FST_file_NAeven_sort <- FST_NA_even[order(FST_NA_even$NA_even, decreasing=TRUE),]
head(FST_file_NAeven_sort)

#top FST loci for NA_even
NAeven_400 <- FST_file_NAeven_sort[1:400,]
NAeven_500 <- FST_file_NAeven_sort[1:500,]


##Find the top FST for the Asia_even
FST_ASIA_even<- FST_file[,c("Locus","ASIA_even")]
head(FST_ASIA_even)

FST_file_ASIAeven_sort <- FST_ASIA_even[order(FST_ASIA_even$ASIA_even, decreasing=TRUE),]
head(FST_file_ASIAeven_sort)

#top FST loci for Asia_even
ASIAeven_400 <- FST_file_ASIAeven_sort[1:400,]
ASIAeven_500 <- FST_file_ASIAeven_sort[1:500,]


###set operations to see how many overlap: 

NA_ASIA_even_400_int <- union(NAeven_400$Locus, ASIAeven_400$Locus)
NA_ASIA_even_500_int <- union(NAeven_500$Locus, ASIAeven_500$Locus)

length(NA_ASIA_even_400_int)

length(NA_ASIA_even_500_int)

# ####Write a list of the loci that were found in the 400 choice for Evens:
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/NA_ASIA_even_400.txt", "wb")
# write.table(NA_ASIA_even_400_int,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
# close(outputFile)


##############ODD LINEAGE PANEL 
names(FST_file)

##Find the top FST for the NA_odd
FST_NA_odd<- FST_file[,c("Locus","NA_odd")]
head(FST_NA_odd)

FST_file_NAodd_sort <- FST_NA_odd[order(FST_NA_odd$NA_odd, decreasing=TRUE),]
head(FST_file_NAodd_sort)

#top FST loci for NA_even
NAodd_400 <- FST_file_NAodd_sort[1:400,]
NAodd_500 <- FST_file_NAodd_sort[1:500,]
NAodd_600 <- FST_file_NAodd_sort[1:600,]


#Find the top FST for the Asia_even
FST_ASIA_odd<- FST_file[,c("Locus","ASIA_odd")]
head(FST_ASIA_odd)

FST_file_ASIAodd_sort <- FST_ASIA_odd[order(FST_ASIA_odd$ASIA_odd, decreasing=TRUE),]
head(FST_file_ASIAodd_sort)

#top FST loci for Asia_even
ASIAodd_400 <- FST_file_ASIAodd_sort[1:400,]
ASIAodd_500 <- FST_file_ASIAodd_sort[1:500,]
ASIAodd_600 <- FST_file_ASIAodd_sort[1:600,]

###set operations to see how many overlap: 
NA_ASIA_odd_400_int <- union(NAodd_400$Locus, ASIAodd_400$Locus)
NA_ASIA_odd_500_int <- union(NAodd_500$Locus, ASIAodd_500$Locus)
NA_ASIA_odd_600_int <- union(NAodd_600$Locus, ASIAodd_600$Locus)


length(NA_ASIA_odd_400_int)
length(NA_ASIA_odd_500_int)
length(NA_ASIA_odd_600_int)


#########overlap between two panels

overlap <- intersect(NA_ASIA_even_400_int, NA_ASIA_odd_600_int)
length(overlap)

# 
# ####Write a list of the loci that were found in the 600 choice for Odds:
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/NA_ASIA_odd_600.txt", "wb")
# write.table(NA_ASIA_odd_600_int,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
# close(outputFile)


############################################# USE the dataframes from above that have the real snps for each lineage to  













#####################################################################################################################################
##################### scrap code to use for identifying true snps in each lineage

Even_tag_pos <- loci_table_t[loci_table_t$Tag%in%just_singletons_tags,]
dim(singletons_to_keep)
head(singletons_to_keep)


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

genotype_file_t<-gsub("X","",colnames(genotype_file))

#split the tags and snp positions
allGenos_oneSNP_temp_tags<-data.frame(str_split_fixed(allGenos_oneSNP_temp$Locus,"_",2))
colnames(allGenos_oneSNP_temp_tags)<-c("Tag","SNP")
head(allGenos_oneSNP_temp_tags)
allGenos_oneSNP_temp$Tag<-allGenos_oneSNP_temp_tags$Tag
allGenos_oneSNP_temp$SNP<-allGenos_oneSNP_temp_tags$SNP
head(allGenos_oneSNP_temp)

