### To select markers for two pink Amplicon Panels, one per lineage
###   Using Ranked FST of one snp per tag loci in a hierarchical structure
###   Remove the paralogs from Haplotype data set for Random Forest   
###   Carolyn Tarpey | January 2018
### ---------------------------------------

#This code requires the FST from Genepop, which was run with the One SNP per TAG data set, 23759 loci. 
#There are three levels to the heirarchy: 
#1. All populations, and All  populations except Susitna 
#2. All Even Populations except Susitna, and all Odd populations except Susitna
#3. Even NA pops w/o Susitna, Odd NA pops w/o Susitna, Even Asia pops, and Odd Asia pops
#Other: Additional level is all the NA pops w/o Susitna and all the Asian pops, 
  #but these are not really for answering out specific marker selection question, just for general interest

#The one SNP per tag genotypes should be filtered through SecondPinkFiltering.R already
#The haplotype file is from STACKS using the 31485 whitelist after initialPinkFiltering.R

###This code also requires a table of all the SNP positions for each of the tags. 
# This is made from the SNP file in the Catalog from STACKS. It has been edited in Excel to have a list of the 
#locus names, as well as the SNP index for each snp at each tag. These are all the SNPS in the whole catalog. 

###This code also requires a population map, which each individual of the analysis in the rows, and their population and lineage named in columns. 
### The lakelse population individuals are named incorrectly, the lineage/year for each of the individuals is switched, 
##The correct population/lineage is named in the pop map


library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)
library(hierfstat)
library(adegenet)
library(argparse)
library(stringi)
library(reshape2)
#install.packages("reshape2")

#______________________________________________________________________________________________
################################IMPORT YOUR DATA FILES
#input the text file that has the results for all the runs of global FST from Genepop
#FST was run on each of the groupings in the three hierarchies listed above, for a total of 10 columns of FST results
FST_file<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Genepop/HWE_Genepop_FST_R.txt", header=TRUE, 
                     stringsAsFactors = FALSE, na.strings = "-" )
FST_file[1:5,1:5]
dim(FST_file)

#input the text file that has the modified SNPs from the Catalog file. 
SNP_positions_all <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/SNPposition/SNP_positions.txt", header=TRUE, 
                                stringsAsFactors = FALSE, na.strings = "-" )
head(SNP_positions_all)
dim(SNP_positions_all)

##import the haplotypes into R to look for fixed snps at each tag in the lineages
Haplotypes<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/SNPposition/Haplotypes_with465ins_singletons.txt", header=TRUE, 
                          stringsAsFactors = FALSE, na.strings = "-" )
Haplotypes[1:5,1:5]
dim(Haplotypes)

#import the population map that has pop and lineage for each sample. it has "ignore" for the susitna lineages
#the lalkelse pops are named their correct lineage, though their sample name doesn't match their population label

popMapLineages<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/SNPposition/PopMapLineages.txt", header=FALSE, 
                          stringsAsFactors = FALSE, na.strings = "-" )
colnames(popMapLineages)<- c("Sample","Pop","Lineage")
head(popMapLineages)
dim(popMapLineages)

#______________________________________________________________________________________________
###################################### SNP Position
#split the haplotype file by lineage, and remove the Susitna pops by ignoring them
##make a dataframe for each lineage with the individuals in that lineage
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
# 
# #####  THESE WORK BUT TAKE FOREVER:: keep them commented out until you need them again
# # First sweep through each list of filtered Haplotypes to remove the genotypes at any individuals with > 2 alleles
# 
# for (i in 1:nrow(Even_filteredHaplotypes)) {
#   row <- as.character(Even_filteredHaplotypes[i,])
#   for (k in 1:length(row)) {
#     if (!is.na(row[k])) {
#       if (str_count(row[k], pattern = "/") > 1) {
#         Even_filteredHaplotypes[i,k] <- NA
#       } else {
#         Even_filteredHaplotypes[i,k] <- Even_filteredHaplotypes[i,k]
#       }
#     }
#   }
# }
# 
# ####Write the even filtered for multiple allele haplotype file in case we need to load it back in:
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Even_filtered_noMultAlleles_Haplo.txt", "wb")
# write.table(Even_filteredHaplotypes,outputFile,quote=FALSE,row.names=TRUE,col.names=TRUE,eol="\n")
# close(outputFile)
# 
# 
# for (i in 1:nrow(Odd_filteredHaplotypes)) {
#   row <- as.character(Odd_filteredHaplotypes[i,])
#   for (k in 1:length(row)) {
#     if (!is.na(row[k])) {
#       if (str_count(row[k], pattern = "/") > 1) {
#         Odd_filteredHaplotypes[i,k] <- NA
#       } else {
#         Odd_filteredHaplotypes[i,k] <- Odd_filteredHaplotypes[i,k]
#       }
#     }
#   }
# }
# 
# ####Write the odd filtered for multiple allele haplotype file in case we need to load it back in:
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Odd_filtered_noMultAlleles_Haplo.txt", "wb")
# write.table(Odd_filteredHaplotypes,outputFile,quote=FALSE,row.names=TRUE,col.names=TRUE,eol="\n")
# close(outputFile)


############# This section takes the haplotypes for each lineage and parses them to see which SNPS are "real" at each tag
############# it does this with a series of if/else statements to judge if the SNP is fixed in the lineage, so it is labled false,
############# or if there is more than 1 example of that same variant it is labeled as a true snp for that lineage, and retained. 

############################# For the EVEN set
haplo <- Even_filteredHaplotypes

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

KEEP[1:20] 
THROW_OUT[1:20]


### EVEN ~~~~~~~~~~take the Keep snp list and throw out snp list and remove the NULLs and DUPLICATES 
KEEP_even_raw <- KEEP
THROW_OUT_even_raw <- THROW_OUT
KEEP_Even_filt_Hap <- KEEP_even_raw[-which(sapply(KEEP_even_raw, is.null))] #remove the NULLS
length(KEEP_Even_filt_Hap)
INT_THROW_Even_filt_Hap <- THROW_OUT_even_raw[-which(sapply(THROW_OUT_even_raw, is.null))] # remove the NULLS
length(INT_THROW_Even_filt_Hap)
INT_THROW_Even_filt_Hap <- as.vector(unique(INT_THROW_Even_filt_Hap))# remove the duplicates
length(INT_THROW_Even_filt_Hap)
THROW_Even_filt_Hap <- setdiff(INT_THROW_Even_filt_Hap,KEEP_Even_filt_Hap) #in INT_throw that were not in KEEP
length(THROW_Even_filt_Hap)


############################# For the Odd set
haplo <- Odd_filteredHaplotypes

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

KEEP[1:20] 
THROW_OUT[1:20]

#### ODD ~~~~~take the Keep snp list and throw out snp list and remove the NULLs and DUPLICATES 
KEEP_odd_raw <- KEEP
THROW_OUT_odd_raw <- THROW_OUT
KEEP_ODD_filt_Hap <- KEEP_odd_raw[-which(sapply(KEEP_odd_raw, is.null))] #remove the NULLS
length(KEEP_ODD_filt_Hap)
INT_THROW_ODD_filt_Hap <- THROW_OUT_odd_raw[-which(sapply(THROW_OUT_odd_raw, is.null))] # remove the NULLS
length(INT_THROW_ODD_filt_Hap)
INT_THROW_ODD_filt_Hap <- as.vector(unique(INT_THROW_ODD_filt_Hap))# remove the duplicates
length(INT_THROW_ODD_filt_Hap)
THROW_ODD_filt_Hap <- setdiff(INT_THROW_ODD_filt_Hap,KEEP_ODD_filt_Hap) #in INT_throw that were not in KEEP
length(THROW_ODD_filt_Hap)

## convert the vector of filtered real SNP indexes to a Data Frame
Keep_even <-KEEP_Even_filt_Hap
Even_filt_snps <-data.frame(stri_list2matrix(Keep_even, byrow=TRUE))
colnames(Even_filt_snps)<-"Even_SNP_Index"
head(Even_filt_snps)

Keep_odd <- KEEP_ODD_filt_Hap
Odd_filt_snps <-data.frame(stri_list2matrix(Keep_odd, byrow=TRUE))
colnames(Odd_filt_snps)<-"Odd_SNP_Index"
head(Odd_filt_snps)


### Find the true SNP position for the variation at each tag. This requires that the SNP_positions_all table be loaded into R, 
#see very top of code for that loading command
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

#make the position numeric so you can compare numbers to numbers
Even_SNPS$Position <- as.numeric(Even_SNPS$Position)
Odd_SNPS$Position <- as.numeric(Odd_SNPS$Position)

##For loops that look over the snp position and flag TRUE for those that are in the range and FALSE for those outside the range
##This uses the absolute column numbers for indexing position in the dataframe, instead of column names
#if the column positions change these will be wrong. 

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

##############Flag TAGS that have a snp with Positions >= 16 and <= 74- we don't want them
#make a new data frame that has the tag and the flag of whether it fits the criteria of SNPS between 17 and 73


##EVEN
head(Even_SNPS)
Even_tags <- unique(Even_SNPS$Tag)
Even_SNPs_InRange <- vector()
head(Even_tags)

for (s in 1:length(Even_tags)) {
  tested_tag <- Even_tags[s]
  trial_set <- Even_SNPS[which(Even_SNPS$Tag %in% tested_tag),]
  if (any(trial_set$Btw_17_73 == FALSE)){
    Even_SNPs_InRange[[s]]<- "FALSE"
  } else {
    Even_SNPs_InRange[[s]]<- "TRUE"
  }
}

###Combine the output into a matrix that has the tag and the flag
Even_Tag_Flags <-as.data.frame(cbind(Even_tags, Even_SNPs_InRange))
colnames(Even_Tag_Flags) <- c("Tag","Even_SNPs_InRange")
head(Even_Tag_Flags)

##ODD
head(Odd_SNPS)
Odd_tags <- unique(Odd_SNPS$Tag)
Odd_SNPs_InRange <- vector()
head(Odd_tags)

for (s in 1:length(Odd_tags)) {
  tested_tag <- Odd_tags[s]
  trial_set <- Odd_SNPS[which(Odd_SNPS$Tag %in% tested_tag),]
  if (any(trial_set$Btw_17_73 == FALSE)){
    Odd_SNPs_InRange[[s]]<- "FALSE"
  } else {
    Odd_SNPs_InRange[[s]]<- "TRUE"
  }
}

###Combine the output into a matrix that has the tag and the flag
#the TRUE are the tags that have snps that we can use to make primers- the GOOD Ones
Odd_Tag_Flags <-as.data.frame(cbind(Odd_tags, Odd_SNPs_InRange))
colnames(Odd_Tag_Flags) <- c("Tag","Odd_SNPs_InRange")
head(Odd_Tag_Flags)

###################################### Marker selection based on ranked FSt of loci 
## we now have a matrix for each lineage that shows whether the snps at a tag are within the range of primer design
## next we want to look at the FST of the snps, and choose 

names(FST_file)
dim(FST_file)

### Add the SNP_inRange designation for the EVEn and ODD lineages to the FST_file 
# split the LOCUS column so that we have a column of Tags to merge with 
newColNames <- c("Tag", "Pos")
newCols <- colsplit(FST_file$Locus, "_", newColNames)
FST_file <- cbind(newCols, FST_file)
head(FST_file)

# merge the SNP_inRange designations
FST_file <- merge(FST_file, Even_Tag_Flags, by = "Tag")
FST_file <- merge(FST_file, Odd_Tag_Flags, by = "Tag")

##########################
##### Create data frames for each of the lineages that has the FST and the SNPs_inRange flag
##Find the top FST for the NA_even
FST_NA_even<- FST_file[,c("Locus","NA_even","Even_SNPs_InRange")]
head(FST_NA_even)
FST_file_NAeven_sort <- FST_NA_even[order(FST_NA_even$NA_even, decreasing=TRUE),]
head(FST_file_NAeven_sort)

##Find the top FST for the NA_even that have SNPS_inRange
FST_file_NAeven_sort_inRange <- FST_file_NAeven_sort[which(FST_file_NAeven_sort$Even_SNPs_InRange == "TRUE"),]
head(FST_file_NAeven_sort_inRange)

##Find the top FST for the NA_odd
FST_NA_odd<- FST_file[,c("Locus","NA_odd","Odd_SNPs_InRange")]
head(FST_NA_odd)
FST_file_NAodd_sort <- FST_NA_odd[order(FST_NA_odd$NA_odd, decreasing=TRUE),]
head(FST_file_NAodd_sort)

##Find the top FST for the NA_odd that have SNPS_inRange
FST_file_NAodd_sort_inRange <- FST_file_NAodd_sort[which(FST_file_NAodd_sort$Odd_SNPs_InRange == "TRUE"),]
head(FST_file_NAodd_sort_inRange)

##Find the top FST for the Asia_even
FST_ASIA_even<- FST_file[,c("Locus","ASIA_even","Even_SNPs_InRange")]
head(FST_ASIA_even)
FST_file_ASIAeven_sort <- FST_ASIA_even[order(FST_ASIA_even$ASIA_even, decreasing=TRUE),]
head(FST_file_ASIAeven_sort)

##Find the top FST for the Asia_even that have SNPS_inRange
FST_file_ASIAeven_sort_inRange <- FST_file_ASIAeven_sort[which(FST_file_ASIAeven_sort$Even_SNPs_InRange == "TRUE"),]
head(FST_file_ASIAeven_sort_inRange)

##Find the top FST for the Asia_odd
FST_ASIA_odd<- FST_file[,c("Locus","ASIA_odd","Odd_SNPs_InRange")]
head(FST_ASIA_odd)
FST_file_ASIAodd_sort <- FST_ASIA_odd[order(FST_ASIA_odd$ASIA_odd, decreasing=TRUE),]
head(FST_file_ASIAodd_sort)

##Find the top FST for the Asia_odd that have SNPS_inRange
FST_file_ASIAodd_sort_inRange <- FST_file_ASIAodd_sort[which(FST_file_ASIAodd_sort$Odd_SNPs_InRange == "TRUE"),]
head(FST_file_ASIAodd_sort_inRange)

######################################
################# EVEN 
#top FST loci for NA_even, 
#change the variable to the number you want 
NAeven_var <- 500 #<--------------------------------------------dictates the number of loci
NAeven_set <- FST_file_NAeven_sort_inRange[1:NAeven_var,]
NAeven_set_avgFST <- mean(NAeven_set$NA_even)
NAeven_set_maxFST <- max(NAeven_set$NA_even)
NAeven_set_minFST <- min(NAeven_set$NA_even)
cat("NAeven FST max, min, average: ", NAeven_set_avgFST, NAeven_set_maxFST, NAeven_set_minFST)

#top FST loci for ASIAeven, 
#change the variable to the number you want 
ASIAeven_var <- 500 #<--------------------------------------------dictates the number of loci
ASIAeven_set <- FST_file_ASIAeven_sort_inRange[1:ASIAeven_var,]
ASIAeven_set_avgFST <- mean(ASIAeven_set$ASIA_even)
ASIAeven_set_maxFST <- max(ASIAeven_set$ASIA_even)
ASIAeven_set_minFST <- min(ASIAeven_set$ASIA_even)
cat("ASIAeven FST max, min, average: ", ASIAeven_set_avgFST, ASIAeven_set_maxFST, ASIAeven_set_minFST)

###set operations to see how many overlap: 
NA_ASIA_even_overlap <- union(NAeven_set$Locus, ASIAeven_set$Locus)
length(NA_ASIA_even_overlap)

## make a new dataframe that has the FST info for the overlapping set 
NA_ASIA_even_overlap_FST <- FST_file[FST_file$Locus %in% NA_ASIA_even_overlap,]  
NA_ASIA_even_overlap_FST_sort <-NA_ASIA_even_overlap_FST[order(NA_ASIA_even_overlap_FST$EVEN, decreasing = TRUE),]
head( NA_ASIA_even_overlap_FST_sort)
dim(NA_ASIA_even_overlap_FST_sort)

# #### Write a list of the loci that were found in the overlapping set for Evens:
# ### Change the name to reflect the number choice
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/NA_ASIA_even_overlap_FST_sort_500_955.txt", "wb")
# write.table(NA_ASIA_even_overlap_FST_sort,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)

################# ODD
#top FST loci for NA_odd,
#change the variable to the number you want 
NAodd_var <- 500 #<--------------------------------------------dictates the number of loci
NAodd_set <- FST_file_NAodd_sort_inRange[1:NAodd_var,]
NAodd_set_avgFST <- mean(NAodd_set$NA_odd)
NAodd_set_maxFST <- max(NAodd_set$NA_odd)
NAodd_set_minFST <- min(NAodd_set$NA_odd)
cat("NAodd FST max, min, average: ", NAodd_set_avgFST, NAodd_set_maxFST, NAodd_set_minFST)

#top FST loci for ASIAodd, 
#change the variable to the number you want 
ASIAodd_var <- 500 #<--------------------------------------------dictates the number of loci
ASIAodd_set <- FST_file_ASIAodd_sort_inRange[1:ASIAodd_var,]
ASIAodd_set_avgFST <- mean(ASIAodd_set$ASIA_odd)
ASIAodd_set_maxFST <- max(ASIAodd_set$ASIA_odd)
ASIAodd_set_minFST <- min(ASIAodd_set$ASIA_odd)
cat("ASIAodd FST max, min, average: ", ASIAodd_set_avgFST, ASIAodd_set_maxFST, ASIAodd_set_minFST)

###set operations to see how many overlap: 
NA_ASIA_odd_overlap <- union(NAodd_set$Locus, ASIAodd_set$Locus)
length(NA_ASIA_odd_overlap)

## make a new dataframe that has the FST info for the overlapping set 
NA_ASIA_odd_overlap_FST <- FST_file[FST_file$Locus %in% NA_ASIA_odd_overlap,]  
NA_ASIA_odd_overlap_FST_sort <-NA_ASIA_odd_overlap_FST[order(NA_ASIA_odd_overlap_FST$ODD, decreasing = TRUE),]
head( NA_ASIA_odd_overlap_FST_sort)
dim(NA_ASIA_odd_overlap_FST_sort)

# #### Write a list of the loci that were found in the overlapping set for Odds:
# ### Change the name to reflect the number choice
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/NA_ASIA_odd_overlap_FST_sort_500_924.txt", "wb")
# write.table(NA_ASIA_odd_overlap_FST_sort,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)

#####################################################
#########overlap between two panels
overlap <- intersect(NA_ASIA_even_overlap, NA_ASIA_odd_overlap)
length(overlap)
head(overlap)

## make a new dataframe that has the FST info for the overlapping set 
Lineage_overlap_FST <- FST_file[FST_file$Locus %in% overlap,]  
Lineage_overlap_FST_sort <-Lineage_overlap_FST[order(Lineage_overlap_FST$ALL_noSusit, decreasing = TRUE),]
head(Lineage_overlap_FST_sort)
dim(Lineage_overlap_FST_sort)

# ####Write a list of the loci that were found to be overlapping between the even and the odd panels:
# ### Change the name of the file to reflect the choice of loci and the final number
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/OverlapBetweenEvenOdd500_197.txt", "wb")
# write.table(Lineage_overlap_FST_sort,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)

##############################
########## top X FST loci that were excluded from the combined analysis because of the SNPS_inRange flag

##pull out the Failed in the SNPS_inRange flag
FST_NAeven_failed <- FST_file_NAeven_sort[which(FST_file_NAeven_sort$Even_SNPs_InRange == "FALSE"),] 
FST_ASIAeven_failed <- FST_file_ASIAeven_sort[which(FST_file_ASIAeven_sort$Even_SNPs_InRange == "FALSE"),]
FST_ASIAodd_failed <- FST_file_ASIAodd_sort[which(FST_file_ASIAodd_sort$Odd_SNPs_InRange == "FALSE"),]
FST_NAodd_failed <- FST_file_NAodd_sort[which(FST_file_NAodd_sort$Odd_SNPs_InRange == "FALSE"),]

##Keep X number of loci that have the highest FST #<--------------------------------------------dictates the number of loci
NAeven_top_failed <- FST_NAeven_failed[1:50,]
ASIAeven_top_failed <- FST_ASIAeven_failed[1:50,]
NAodd_top_failed <- FST_NAodd_failed[1:50,]
ASIAodd_top_failed <- FST_ASIAodd_failed[1:50,]

head(NAeven_top_failed)
head(ASIAeven_top_failed)
head(NAodd_top_failed)
head(ASIAodd_top_failed)

##merge the FST_file the SNP_positions_all to see what the SNP positions are for the loci that failed. 
masterFST_file <- merge(SNP_positions_all, FST_file, by= "Tag")
masterFST_file <- rename(masterFST_file, Locus = Locus.x, FST_Locus = Locus.y)
head(masterFST_file)

#trim columns we dont need
masterFST_file$Pos <-NULL
masterFST_file$SNP_Index <-NULL
masterFST_file$SNP_tag_index <-NULL
head(masterFST_file)

###Pull out the positions of the SNPS in the top x number of loci with the highest FST that failed the SNP postion flag
FST_NAeven_top_failed_snps <- masterFST_file[masterFST_file$FST_Locus %in% NAeven_top_failed$Locus,] 
FST_ASIAeven_top_failed_snps <- masterFST_file[masterFST_file$FST_Locus %in% ASIAeven_top_failed$Locus,]
FST_ASIAodd_top_failed_snps <- masterFST_file[masterFST_file$FST_Locus %in% ASIAodd_top_failed$Locus,]
FST_NAodd_top_failed_snps <- masterFST_file[masterFST_file$FST_Locus %in% NAodd_top_failed$Locus,]

head(FST_NAeven_top_failed_snps)
dim(FST_NAeven_top_failed_snps)
# ### Change the name of the file to reflect the choice of loci 
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/FST_NAeven_top100_failed_snps.txt", "wb")
# write.table(FST_NAeven_top_failed_snps,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)

head(FST_ASIAeven_top_failed_snps)
dim(FST_ASIAeven_top_failed_snps)
# ### Change the name of the file to reflect the choice of loci 
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/FST_ASIAeven_top100_failed_snps.txt", "wb")
# write.table(FST_ASIAeven_top_failed_snps,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)

head(FST_ASIAodd_top_failed_snps)
dim(FST_ASIAodd_top_failed_snps)
# ### Change the name of the file to reflect the choice of loci 
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/FST_ASIAodd_top100_failed_snps.txt", "wb")
# write.table(FST_ASIAodd_top_failed_snps,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)

head(FST_NAodd_top_failed_snps)
dim(FST_NAodd_top_failed_snps)
# ### Change the name of the file to reflect the choice of loci 
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/FST_NAodd_top100_failed_snps.txt", "wb")
# write.table(FST_NAodd_top_failed_snps,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)

####### overlap of the Even NA and ASian failed 
Even_failed <- intersect(NAeven_top_failed$Locus, ASIAeven_top_failed$Locus)
length(Even_failed)
head(Even_failed)

####### overlap of the Odd NA and ASian failed 
Odd_failed<- intersect(NAodd_top_failed$Locus, ASIAodd_top_failed$Locus)
length(Odd_failed)
head(Odd_failed)

#########################################################################
#### Use the combination of the 50 top failed from each of the groups for the lineage. So take the top 50 from each and then do a union
####### overlap of the Even NA and ASian failed 
Even_failed_union <- union(NAeven_top_failed$Locus, ASIAeven_top_failed$Locus)
length(Even_failed_union)
head(Even_failed_union)

####### overlap of the Odd NA and Asian failed 
Odd_failed_union<- union(NAodd_top_failed$Locus, ASIAodd_top_failed$Locus)
length(Odd_failed_union)
head(Odd_failed_union)

## Add these list of failed to the list of passed to get a list to pull the sequences to make a FASTA file 
Even_failed_passed_list <- append(NA_ASIA_even_overlap_FST_sort$Locus, Even_failed_union)
head(Even_failed_passed_list)
length(Even_failed_passed_list)
### Change the name of the file 
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Even_failed_passed_list.txt", "wb")
# write.table(Even_failed_passed_list,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
# close(outputFile)

Odd_failed_passed_list <- append(NA_ASIA_odd_overlap_FST_sort$Locus, Odd_failed_union)
head(Odd_failed_passed_list)
length(Odd_failed_passed_list)
### Change the name of the file 
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Odd_failed_passed_list.txt", "wb")
# write.table(Odd_failed_passed_list,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
# close(outputFile)

###########Make a table that has the FST, Locus_FST, SNP Pos and SNP Flag for all the pass and failed loci tog. 
#combine the two failed lists, so both lineages in one 
Both_failed_union <- union(Odd_failed_union, Even_failed_union)
length(Both_failed_union)

#combine ALL THE LISTS, failed and passed, so both lineages in one 
masterCombined_BothLineages <- union(Even_failed_passed_list, Odd_failed_passed_list)
length(masterCombined_BothLineages)


#Make the COMBINED loci list table 
FST_file_COMBINED_table <- FST_file[FST_file$Locus %in% masterCombined_BothLineages,] 
head(FST_file_COMBINED_table)

### Change the name of the file 
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/FST_file_COMBINED_table_1809.txt", "wb")
# write.table(FST_file_COMBINED_table,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)


#Make the failed loci their own table that has the many SNPS that they have. 
masterFST_file_failed_union <- masterFST_file[masterFST_file$FST_Locus %in% Both_failed_union,] 
head(masterFST_file_failed_union)

### Change the name of the file 
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/masterFST_file_failed_union.txt", "wb")
# write.table(masterFST_file_failed_union,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)
