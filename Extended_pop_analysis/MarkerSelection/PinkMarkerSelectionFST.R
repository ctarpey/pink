### To select markers for two pink Amplicon Panels, one per lineage
###   Using Ranked FST of one snp per tag loci in a hierarchical structure
###    Identify which SNPs are real variation in each lineage using Haplotype file 
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
library(RColorBrewer)
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
# 
# ##import the test haplotypes into R to look for fixed snps at each tag in the lineages <--------------------------------------TEST SET
# Test_Haplotypes<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Haplotypes_R_test.txt", header=TRUE, 
#                        stringsAsFactors = FALSE, na.strings = "-" )
# Test_Haplotypes[1:5,1:5]
# dim(Test_Haplotypes)

#import the population map that has pop and lineage for each sample. it has "ignore" for the susitna lineages
#the lalkelse pops are named their correct lineage, though their sample name doesn't match their population label

popMapLineages<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/SNPposition/PopMapLineages.txt", header=FALSE, 
                          stringsAsFactors = FALSE, na.strings = "-" )
colnames(popMapLineages)<- c("Sample","Pop","Lineage")
head(popMapLineages)
dim(popMapLineages)

#_______________________________ The below are not required to run the script the first time
# 
# # ##Input the even filtered for multiple allele haplotype file (this is made further down in this script, but is a pain to make):
# Even_filteredHaplotypes <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Even_filtered_noMultAlleles_Haplo.txt", header=TRUE,
#                                       stringsAsFactors = FALSE, row.names=1)
# Even_filteredHaplotypes[1:5,1:5]
# 
# # ##Input the even filtered for multiple allele haplotype file (this is made further down in this script, but is a pain to make):
# Odd_filteredHaplotypes <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Odd_filtered_noMultAlleles_Haplo.txt", header=TRUE,
#                                       stringsAsFactors = FALSE, row.names=1)
# Odd_filteredHaplotypes[1:5,1:5]
#####_________________________________________________________________________________ 

hap_ids <- Haplotypes$Catalog_Ids
fst_ids <- FST_file$Tag

length(hap_ids)
length(fst_ids)

hap_fst_intersect <- intersect(hap_ids, fst_ids)
length(hap_fst_intersect)

inhap_notfst <- setdiff(hap_ids, fst_ids)
length(inhap_notfst)

inhap_notfst

infst_nothap <- setdiff(fst_ids, hap_ids)
length(infst_nothap)


#####_________________________________________________________________________________ 
###################################### Prepping Haplotypes 
# The following section takes a Haplotype file from populations in STACKs and splits it up by lineage
# Then identifies any haplotype that is not biallelic and removes that haplotype at that individual. 

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

Odd_filteredHaplotypes <-Haplotypes[,colnames(Haplotypes)%in%Odd_inds$Sample]
Odd_filteredHaplotypes[1:5,1:5]
dim(Odd_filteredHaplotypes)

is.data.frame(Even_filteredHaplotypes)

# ################### For the test section    <---------------------------------------------------------------------TEST SET
# ##filter the test haplotype file by lineage
# Test_Even_filteredHaplotypes <-Test_Haplotypes[,colnames(Test_Haplotypes)%in%Even_inds$Sample]
# Test_Even_filteredHaplotypes[1:5,1:5]
# dim(Test_Even_filteredHaplotypes)
# names(Test_Even_filteredHaplotypes)
# 
# Test_Odd_filteredHaplotypes <-Test_Haplotypes[,colnames(Test_Haplotypes)%in%Odd_inds$Sample]
# Test_Odd_filteredHaplotypes[1:5,1:5]
# dim(Test_Odd_filteredHaplotypes)
# names(Test_Odd_filteredHaplotypes)

####  THESE WORK BUT TAKE FOREVER:: keep them commented out until you need them again
# First sweep through each list of filtered Haplotypes to remove the genotypes at any individuals with > 2 alleles
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
# #
# 
# ####Write a table to a file
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Even_filtered_noMultAlleles_Haplo.txt", "wb")
# write.table(Even_filteredHaplotypes,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)
# 
# 
# ####Write a table to a file
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Odd_filtered_noMultAlleles_Haplo.txt", "wb")
# write.table(Odd_filteredHaplotypes,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)


###################################### SNP Position
### This section takes the haplotypes for each lineage at each individual and identifies which SNPS are "real" 
### at each position of the haplotype it does this with a series of if/else statements to judge if the SNP is 
### fixed in the lineage, if the snps is fixed, or has no variation in any individual, the name of the tag 
### and the position of the snp is added to a list of snps to throw out.If there is more than 1 example of that 
### same variant at that position in any individual, it is labeled as a true snp for that lineage, and the name 
### of the tag and the position of the snp is added to a list of snps to retained. 
### The Even Haplotype file and the Odd Haplotype file are run through the for loop seperately to get the real snps
### in each lineage. The for loop only has generic names in it, and to change which Haplotype file you are running, 
### assign it to haplo below. 

############################# For the EVEN set  #Uncomment the haplo below to run the Even , and stop at line 
haplo <- Even_filteredHaplotypes
############################# For the Odd set
haplo <- Odd_filteredHaplotypes
haplo[1:5,1:5]
#test
#haplo <- Odd_filteredHaplotypes[1:10,]
#looking for a specific haplotype in the haplotype file
#which(Haplotypes$Catalog_Ids == 12823) #this is the catalog ID #
#Haplotypes[3806,] #put the row number here and it will output the haplotype for that catalog id

base_set <- c("A", "C", "G", "T")
KEEP <- vector("list",200000) #these are the snps (tag_position) that are showing <2 of the same variation in the lineage
THROW_OUT <- vector("list",200000 ) #these are the snps that are not showing variation, (or only one instance of of that variation) or are fixed- no variation
FIXXED <- vector("list",200000) #These are the snps that have no variation
j = 1 #index for adding things to the list of keep, Throw out, or Fixxed


for (i in 1:nrow(haplo)) { #for each row in the haplotype file
  row <- as.character(haplo[i,])
  #print(row)
  row_number <- i #this is the row index of the file, not the tag
  tag <- Haplotypes[i,"Catalog_Ids"] #this pulls the tag number from the original haplotype file 
  #tag <- row_number #<-------------------------------------- for the test it just uses the row number of the test file
  #geno_count <- sum(!is.na(row)) #Count the number of genotypes at that tag
  #print(geno_count)
  min_allele <- min(nchar(row, type = "chars")) #this doesn't ignore NA, but it they have 2 chars and we're only looking for haplotypes with 1 char
  #print(min_allele)
  #max_allele <- max(nchar(row, type = "chars")) #  dont use this, but it is interesting- here it includes characters in all alleles
  #print(max_allele)
  
  if (min_allele == 1) {  #if the haplotype has one character it is a single snp per tag 
    hets <- sum(grepl("/", row)) # count the number of / in the row
    #are all the letters the same? if so, that snp is fixed in all the idividuals and we want to throw it out and add it to the fixxed list. 
    #<-- Add a section that counts instances of the het bases, if they are all different then it is not useful variation and we dont want it
    if (hets >= 2) { # if there are more than 2 /, it is het and we want to keep it
      KEEP[[j]] <- paste(tag, "1", sep="_") #add it to the Keep list with a # of 1 because there is only one position in these haplotypes
      j = j + 1 #this is the index for adding things to the list, tags/ position would overlap
      #next() #<---------------------------------- DO I NEED THIS?? WHen it finishes this loop I want it to go to the next row of the haplotype file
    } else {  #each individual has one allele in that row of the haplotype file, but we don't know if it is variabl between individuals
      counts <- table(na.omit(row)) #Counts the number of each base 
      if (length(counts) == 1 ) { #If there is only one base, no variation, throw out and add to fixed list
        THROW_OUT[[j]] <- paste(tag, "1", sep="_") #add it to the throw out list 
        FIXXED[[j]]<- paste(tag, "1", sep="_") # add it to a list of tags that are fixed in the lineage
        j = j + 1 
        #next() #<---------------------------------- DO I NEED THIS?? WHen it finishes this loop I want it to go to the next row of the haplotype file
      } else if (length(counts)== 2 & any(counts < 2)) { #this catches cases where there is one Het- it will be its own count on the table and would be a count of one
        THROW_OUT[[j]] <- paste(tag, "1", sep="_") #add it to the throw out list 
        j = j + 1
        #next() #<---------------------------------- DO I NEED THIS?? WHen it finishes this loop I want it to go to the next row of the haplotype file
      } else if (length(counts)== 3 & ((counts[1]<2) + (counts[2]<2) + (counts[3]<2)) == 2) { #if two of these bases occur only once, not enough variation, throw out. 
        THROW_OUT[[j]] <- paste(tag, "1", sep="_") #add it to the throw out list 
        j = j + 1 
        #next() #<---------------------------------- DO I NEED THIS?? WHen it finishes this loop I want it to go to the next row of the haplotype file
      } else if (length(counts)== 4 & ((counts[1]<2) + (counts[2]<2) + (counts[3]<2)+ (counts[4]<2)) == 3) { #if three of these bases occur only once, not enough variation, throw out.
        THROW_OUT[[j]] <- paste(tag, "1", sep="_") #add it to the throw out list 
        j = j + 1 
        #next() #<---------------------------------- DO I NEED THIS?? WHen it finishes this loop I want it to go to the next row of the haplotype file
      } else { #if it doesn't fit the criteria above then it should be kept, it has enough variation. 
        KEEP[[j]] <- paste(tag, "1", sep="_") #add it to the Keep list 
        j = j + 1 
        #next() #<---------------------------------- DO I NEED THIS?? WHen it finishes this loop I want it to go to the next row of the haplotype file
      }
    }
  } else if (min_allele != 1) { # If there is not just one base at the haplotype
    row <- str_split_fixed(row,"/",Inf)
    a_test <- table(row)
    if (length(a_test) == 1){
      for (p in 1:nchar(names(a_test[1]))){ ###<- For tags that have fixed homozygous variation, ie AA in all individuals, names(table[1]) returns the name, not ""
        THROW_OUT[[j]] <- paste(tag, p, sep="_") #add every position of the haplotype to the throw out list
        FIXXED[[j]]<- paste(tag, p, sep="_")# add every position of the haplotype to a list of tags that are fixed in the lineage
        j = j + 1 #index to prevent tags_positions from overwriting in the lists
        #print(paste(tag, p, sep="_"))
      } #m<- 1
    } else { 
      alleles <- names(a_test[2:length(a_test)]) #Gives each of the variants ie: "AGGAA" "AGAAT" "GGAAA"
      SNPs <- nchar(alleles[1]) # this counts the number of positions in the first allele name, which should be standard for all
      for(m in 1:(SNPs)){ #for every position in the number of positions in the haplotype 
        A <- 0        #restart the counter for the bases
        C <- 0
        G <- 0
        T <- 0 
        for (l in 1:length(alleles)){ # for every variant, ie:"AGGAA" "AGAAT" "GGAAA" 
          var <- substr(alleles[l], m, m) #pull out the base at the m position in the l# variant, like the bases at the second position of AGGAA = G
          if (var == "A"){  #Count that base if it is an A, (other wise count it where it should be counted, then go to the next allele and count the base at m position in the same way)
            A <- A+1 
          } else if (var == "C"){
            C <- C+1 
          } else if (var == "G"){
            G <- G+1 
          } else if (var == "T"){
            T <- T+1 
          } 
        }# <----------------------------Here I want to make sure that it has counted every base at a particular position (m) for all the alleles before going on to the next loop part
        if (length(which(c(A, C, G, T) == 0)) == 3){ #If there is only one type of base at this position in all alleles, it is fixed, and we dont want it. 
          THROW_OUT[[j]] <- paste(tag, m ,sep="_") #add it to the throw out list 
          FIXXED[[j]]<- paste(tag, m, sep="_")# add it to a list of tags that are fixed in the lineage 
          j = j + 1 
          #next() #<---------------I'm worried that when the first position is fixed, it is writing it to a list but isnt going to the next position. 
        } else if (length(which(c(A, C, G, T) >= 2)) >= 2){ # If there are more than 2 with 2 or more instances of that base, we want it. 
          KEEP[[j]] <- paste(tag, m, sep="_") #add it to the Keep list 
          j = j + 1 
          ##the next line is the ambiguous situation, if a base has a count of 1 in the table above, say A is 1, for position 1 that means it was only in one 
          #allele (ie AGGAA)at that position but there might be 32 of that allele,in that case we want to keep that snp 
          # or if there is only one AGGAA, all the other alleles have G in position 1 then it isn't enough variation and we dont want it
        } else if ((length(which(c(A, C, G, T) >= 2)) == 1) & (length(which(c(A, C, G, T) == 1)) >= 1)){
          if (length(which(c(A, C, G, T) == 1)) == 1){ 
            base <- base_set[which(c(A, C, G, T) == 1)] #gives the base we need to search for the counts for ie "A"
            for (n in 1:length(alleles)){ #go to each allele in turn and look for the position we want 
              if ((substr(alleles[n], m, m) == base) & (as.vector(a_test[n+1]) >=2)) { # and return the base at that position (m), and then go look for the count of the allele that the base came from in the table a_test, if there are more than 2 of it, keep the snp
                KEEP[[j]] <- paste(tag, m, sep="_") #add it to the Keep list 
                j = j + 1 
              } else if ((substr(alleles[n], m, m) == base) & (as.vector(a_test[n+1]) <2)){ #if there was not more than one instance of that allele, throw out the snp
                THROW_OUT[[j]] <- paste(tag, m, sep="_") #add it to the throw out list 
                j = j + 1 
              } else {
                next() # <- I want this to take us to the top of this loop, line 261, to evaluate the next position in the alleles, but <----I DONT THINK THIS WORKS
              }
            }
          } else if ((length(which(c(A, C, G, T) >= 2)) == 1) & (length(which(c(A, C, G, T) == 1)) == 3)){ # this is the same ambiguous situation as above, but with 3 ambiguous snp counts, like A and G and t are both 1, but they might come from alleles with counts 2 or more
            base <- base_set[which(c(A, C, G, T) == 1)] #gives the base we need to search for ie "A"
            for (n in 1:length(alleles)){
              if ((substr(alleles[n], m, m) == base[1]) & (as.vector(a_test[n+1]) >=2)) {
                KEEP[[j]] <- paste(tag, m, sep="_") #add it to the Keep list 
                j = j + 1 
              } else if ((substr(alleles[n], m, m) == base[2]) & (as.vector(a_test[n+1]) >=2)) {
                KEEP[[j]] <- paste(tag, m, sep="_") #add it to the Keep list 
                j = j + 1 
              } else if ((substr(alleles[n], m, m) == base[3]) & (as.vector(a_test[n+1]) >=2)) {
                KEEP[[j]] <- paste(tag, m, sep="_") #add it to the Keep list 
                j = j + 1 
              } else {
                THROW_OUT[[j]] <- paste(tag, m, sep="_") #add it to the throw out list  # if it doesnt fit any of these instances throw it out. 
                j = j + 1 
              }
            }
          }
        } else if ((length(which(c(A, C, G, T) == 1)) == 2)) { # this is the same ambiguous situation as above, but with 2 ambiguous snp counts, like A and G are both 1, but they might come from alleles with counts 2 or more
          base <- base_set[which(c(A, C, G, T) == 1)] #gives the base we need to search for ie "A"
          base_1 <- which(substr(alleles, m, m) == base[1])
          base_2 <- which(substr(alleles, m, m) == base[2])
          if ((as.vector(a_test[base_1+1]) >=2) & (as.vector(a_test[base_2+1]) >=2)) { # If the counts for the allele that the first ambiguous snp came from are greater than 2, keep it
            KEEP[[j]] <- paste(tag, m, sep="_")
            j = j + 1 
            next() #<- if we kept it we want to go back to line 261 and start the next position in the allele
          } else {
            THROW_OUT[[j]] <- paste(tag, m, sep="_") #add it to the throw out list 
            j = j + 1 
          }
        }
      }
    }
  }
}

KEEP[1:20] 
THROW_OUT[1:20]
FIXXED[1:20]


##EVEN ~~~~~~~~~~take the Keep snp list and throw out snp list and remove the NULLs and DUPLICATES
# KEEP_even_raw <- KEEP #<---------------------Uncomment this and the next 2 lines when you are running the function for Even
# THROW_OUT_even_raw <- THROW_OUT #saves the output as lineage specific
# FIXXED_even_raw <- FIXXED

KEEP_Even_filt_Hap <- KEEP_even_raw[-which(sapply(KEEP_even_raw, is.null))] #remove the NULLS
length(KEEP_Even_filt_Hap)
KEEP_Even_filt_Hap <- as.vector(unique(KEEP_Even_filt_Hap)) # test to remove Duplicates 
length(KEEP_Even_filt_Hap)
INT_THROW_Even_filt_Hap <- THROW_OUT_even_raw[-which(sapply(THROW_OUT_even_raw, is.null))] # remove the NULLS
length(INT_THROW_Even_filt_Hap)
INT_THROW_Even_filt_Hap <- as.vector(unique(INT_THROW_Even_filt_Hap)) # test to remove Duplicates 
length(INT_THROW_Even_filt_Hap)
THROW_Even_filt_Hap <- setdiff(INT_THROW_Even_filt_Hap,KEEP_Even_filt_Hap) #test to see if there is overlap between lists
length(THROW_Even_filt_Hap)
FIXXED_Even_filt_Hap <- FIXXED_even_raw[-which(sapply(FIXXED_even_raw, is.null))] #remove the NULLS
length(FIXXED_Even_filt_Hap)
FIXXED_Even_filt_Hap <- as.vector(unique(FIXXED_Even_filt_Hap))# test to remove Duplicates 
length(FIXXED_Even_filt_Hap)
#head(FIXXED_Even_filt_Hap)

#__________________________________________________________________ Stop Here and Return to line 174 and re-run the Function with the ODD before continuing

# ## ODD ~~~~~take the Keep snp list and throw out snp list and remove the NULLs and DUPLICATES
# KEEP_odd_raw <- KEEP #<--------------------- #Uncomment this and the next 2 lines when you are running the function for ODD
# THROW_OUT_odd_raw <- THROW_OUT #saves the output as lineage specific
# FIXXED_ODD_raw <- FIXXED

KEEP_ODD_filt_Hap <- KEEP_odd_raw[-which(sapply(KEEP_odd_raw, is.null))] #remove the NULLS
length(KEEP_ODD_filt_Hap)
KEEP_ODD_filt_Hap <- as.vector(unique(KEEP_ODD_filt_Hap)) # test to remove Duplicates 
length(KEEP_ODD_filt_Hap)
INT_THROW_ODD_filt_Hap <- THROW_OUT_odd_raw[-which(sapply(THROW_OUT_odd_raw, is.null))] # remove the NULLS
length(INT_THROW_ODD_filt_Hap)
INT_THROW_ODD_filt_Hap <- as.vector(unique(INT_THROW_ODD_filt_Hap))# test to remove Duplicates 
length(INT_THROW_ODD_filt_Hap)
THROW_ODD_filt_Hap <- setdiff(INT_THROW_ODD_filt_Hap,KEEP_ODD_filt_Hap)  #test to see if there is overlap between lists
length(THROW_ODD_filt_Hap) 
FIXXED_ODD_filt_Hap <- FIXXED_ODD_raw[-which(sapply(FIXXED_ODD_raw, is.null))] #remove the NULLS
length(FIXXED_ODD_filt_Hap)
FIXXED_ODD_filt_Hap <- as.vector(unique(FIXXED_ODD_filt_Hap))# test to remove Duplicates 
length(FIXXED_ODD_filt_Hap)
#head(FIXXED_ODD_filt_Hap)


###+++++++++++++++++++++++++++++++++++++++++++ ONLY After the above has been run 
## Convert the above vectors for KEEP, THROW OUT and FIXXEd to Data Frames

## convert the vector of filtered SNP indexes to a Data Frame
Keep_even <-KEEP_Even_filt_Hap
length(Keep_even)
Even_filt_snps <-data.frame(stri_list2matrix(Keep_even, byrow=TRUE))
colnames(Even_filt_snps)<-"Even_SNP_Index"
head(Even_filt_snps)
length(Even_filt_snps$Even_SNP_Index)
# 

Keep_odd <- KEEP_ODD_filt_Hap
length(Keep_odd)
Odd_filt_snps <-data.frame(stri_list2matrix(Keep_odd, byrow=TRUE))
colnames(Odd_filt_snps)<-"Odd_SNP_Index"
head(Odd_filt_snps)
dim(Odd_filt_snps)


## convert the vector of filtered Thrown out SNP indexes to a Data Frame
Throw_even <-THROW_Even_filt_Hap
Even_filt_throw_snps <-data.frame(stri_list2matrix(Throw_even, byrow=TRUE))
colnames(Even_filt_throw_snps)<-"Even_SNP_throw_Index"
head(Even_filt_throw_snps)

Throw_odd <- THROW_ODD_filt_Hap
Odd_filt_throw_snps <-data.frame(stri_list2matrix(Throw_odd, byrow=TRUE))
colnames(Odd_filt_throw_snps)<-"Odd_SNP_throw_Index"
head(Odd_filt_throw_snps)

## convert the vector of filtered fixxed SNP indexes to a Data Frame
FIXXED_even <-FIXXED_Even_filt_Hap
Even_fixxed_snps <-data.frame(stri_list2matrix(FIXXED_even, byrow=TRUE))
colnames(Even_fixxed_snps)<-"Even_SNP_Index"
head(Even_fixxed_snps)

FIXXED_odd <- FIXXED_ODD_filt_Hap
Odd_fixxed_snps <-data.frame(stri_list2matrix(FIXXED_odd, byrow=TRUE))
colnames(Odd_fixxed_snps)<-"Odd_SNP_Index"
head(Odd_fixxed_snps)

### Find the true SNP position for the variation at each tag. This requires that the SNP_positions_all table be loaded into R, 
#see very top of code for that loading command
head(SNP_positions_all) 
dim(SNP_positions_all)

dim(Even_filt_snps)
dim(Odd_filt_snps)
head(Odd_filt_snps)

Even_SNPS <- SNP_positions_all[which(SNP_positions_all$SNP_tag_index %in% Even_filt_snps$Even_SNP_Index),]
Odd_SNPS <- SNP_positions_all[which(SNP_positions_all$SNP_tag_index %in% Odd_filt_snps$Odd_SNP_Index),]

dim(Even_SNPS)
dim(Odd_SNPS)
head(Even_SNPS)
head(Odd_SNPS)

# ####Write the Even_SNPS
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Even_SNPS.txt", "wb")
write.table(Even_SNPS,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

####Write the Odd_SNPS
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Odd_SNPS.txt", "wb")
write.table(Odd_SNPS,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

Even_throw_SNPS <- SNP_positions_all[which(SNP_positions_all$SNP_tag_index %in% Even_filt_throw_snps$Even_SNP_throw_Index),]
Odd_throw_SNPS <- SNP_positions_all[which(SNP_positions_all$SNP_tag_index %in% Odd_filt_throw_snps$Odd_SNP_throw_Index),]
dim(Even_throw_SNPS)
dim(Odd_throw_SNPS)
head(Even_throw_SNPS)
head(Odd_filt_throw_snps)

# ####Write the Even_throw_SNPS
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Even_throw_SNPS.txt", "wb")
write.table(Even_throw_SNPS,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

####Write the Odd_throw_SNPS
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Odd_throw_SNPS.txt", "wb")
write.table(Odd_throw_SNPS,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

Even_fixxed_SNPS <- SNP_positions_all[which(SNP_positions_all$SNP_tag_index %in% Even_fixxed_snps$Even_SNP_Index),]
Odd_fixxed_SNPS <- SNP_positions_all[which(SNP_positions_all$SNP_tag_index %in% Odd_fixxed_snps$Odd_SNP_Index),]
dim(Even_fixxed_SNPS)
dim(Odd_fixxed_SNPS)
head(Even_fixxed_SNPS)
head(Odd_fixxed_SNPS)

# ####Write the Even_fixxed_SNPS
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Even_fixxed_SNPS.txt", "wb")
write.table(Even_fixxed_SNPS,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

####Write the Odd_fixxed_SNPS
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Odd_fixxed_SNPS.txt", "wb")
write.table(Odd_fixxed_SNPS,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

#####################################  Stats on the differences between the lineages: number of variable sites per tag, and # of tags. 
head(Even_SNPS)
head(Odd_SNPS)
Histo_7_colors<- c("#6e016b","#88419d","#8c6bb1","#8c96c6","#7bccc4",'#4eb3d3','#2b8cbe')
Histo_14_colors<- c("#6e016b","#6e016b","#88419d","#88419d","#8c6bb1","#8c6bb1","#8c96c6","#8c96c6","#7bccc4","#7bccc4",'#4eb3d3','#4eb3d3','#2b8cbe','#2b8cbe')
E_O_Colors <- c("#00158a", "#c9cf4a")
  
even_snps_counts <- as.data.frame(table(Even_SNPS$Tag))
even_snps_counts <- as.data.frame(table(even_snps_counts$Freq))
colnames(even_snps_counts) <-c("SNPs","Even_Count")
even_snps_counts

odd_snps_counts <- as.data.frame(table(Odd_SNPS$Tag))
odd_snps_counts <- as.data.frame(table(odd_snps_counts$Freq))
colnames(odd_snps_counts) <-c("SNPs","Odd_Count")
odd_snps_counts

snp_counts<-merge(odd_snps_counts,even_snps_counts, by = "SNPs")
snp_counts
snp_counts_melt <- melt(snp_counts)
snp_counts_melt

#pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLOTs_Snps_per_Haplotype_lineage.pdf", width = 9, height = 7)

#individual histograms of even and odd
ggplot(data = even_snps_counts, aes(x=SNPs, y= Even_Count)) + theme_bw() + 
  geom_bar(stat= "identity", fill = Histo_7_colors)  + xlab("Number of SNPs per Haplotype") + ylab("Count") +
  ggtitle("Variable SNPs per Haplotype in the Even Lineage")+ scale_y_continuous(limits=c(0, 15000)) 

ggplot(odd_snps_counts, aes(x=SNPs, y= Odd_Count)) + theme_bw() + 
  geom_bar(stat= "identity", fill=Histo_7_colors) +xlab("Number of SNPs per Haplotype") + ylab("Count") +
  ggtitle("Variable SNPs per Haplotype in the Odd Lineage") + scale_y_continuous(limits=c(0, 15000))

# histogram with the even and odd next to eachother in one plot 
ggplot(data=snp_counts_melt, aes(snp_counts_melt$SNPs, snp_counts_melt$value)) + theme_bw() +
  geom_bar(aes(fill= snp_counts_melt$variable), alpha = 0.5, position = "dodge", stat="identity") +
   xlab("Number of SNPs per Haplotype") + ylab("Count") + ggtitle("Variable SNPs per Haplotype ") +
  scale_fill_manual(values = E_O_Colors, labels = c("Odd", "Even")) + labs(fill="Lineage")

#dev.off()  

########################### More filtering of SNPS based on position in the tag. 
###Flag SNP Positions >= 16 and <= 74, (we want tags 17-73 to be able to make primers)
dim(Even_SNPS)

#add new empty column to the end of each data frame for the flags
Even_SNPS$Btw_17_73 <- NA
Odd_SNPS$Btw_17_73 <- NA

#make the position numeric so you can compare numbers to numbers
Even_SNPS$Position <- as.numeric(Even_SNPS$Position)
Odd_SNPS$Position <- as.numeric(Odd_SNPS$Position)

##For loops that look over the snp position for the SNPs that have been identified as variable in the lineage
# and flag TRUE for those that are in the range and FALSE for those outside the range
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
dim(Even_SNPS)
head(Odd_SNPS)
dim(Odd_SNPS)


##############Flag TAGS that have a snp with Positions >= 16 and <= 74- we don't want them
#make a new data frame that has the tag and the flag of whether it fits the criteria of SNPS between 17 and 73


##EVEN
head(Even_SNPS)
Even_tags <- unique(Even_SNPS$Tag)
length(Even_tags)
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

#Combine the output into a matrix that has the tag and the flag
Even_Tag_Flags <-as.data.frame(cbind(Even_tags, Even_SNPs_InRange))
colnames(Even_Tag_Flags) <- c("Tag","Even_SNPs_InRange")
head(Even_Tag_Flags)

dim(Even_Tag_Flags)
length(which(Even_Tag_Flags$Even_SNPs_InRange=="TRUE"))
length(which(Even_Tag_Flags$Even_SNPs_InRange=="FALSE"))


# #### Write the Even tag flags:
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Even_tag_flags.txt", "wb")
# write.table(Even_Tag_Flags,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)


##ODD
head(Odd_SNPS)
Odd_tags <- unique(Odd_SNPS$Tag)
length(Odd_tags)
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

#Combine the output into a matrix that has the tag and the flag
#the TRUE are the tags that have snps that we can use to make primers- the GOOD Ones
Odd_Tag_Flags <-as.data.frame(cbind(Odd_tags, Odd_SNPs_InRange))
colnames(Odd_Tag_Flags) <- c("Tag","Odd_SNPs_InRange")
head(Odd_Tag_Flags)

dim(Odd_Tag_Flags)
length(which(Odd_Tag_Flags$Odd_SNPs_InRange=="TRUE"))
length(which(Odd_Tag_Flags$Odd_SNPs_InRange=="FALSE"))

# #### Write the Odd tag flags:
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Odd_Tag_Flags.txt", "wb")
# write.table(Odd_Tag_Flags,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)


################### This section is similar to the above, but it is just for the benefit of Garrett's primer design pipeline. 

Even_SNPS_PD <- Even_SNPS
Odd_SNPS_PD <- Odd_SNPS

#add two new empty columns to the end of each data frame for the flags of the left and right limits- needed for Garrett's primer pipeline
Even_SNPS_PD$LeftLimit<- NA
Even_SNPS_PD$RightLimit<- NA

Odd_SNPS_PD$LeftLimit <- NA
Odd_SNPS_PD$RightLimit <- NA

##For loops that look over the snp position for the SNPs that have been identified as variable in the lineage
# and look at whether the snp is at a position smaller than 17 and flags it 0 is yes, 1 is no
# Then it looks at whether the snps is at a position that is larger than 73 and flags it 0 is yes, 1 is no
##This uses the absolute column numbers for indexing position in the dataframe, instead of column names
#if the column positions change these will be wrong. 

for (r in 1:dim(Even_SNPS_PD)[1]){
  if (Even_SNPS_PD[r,3] < 17){
    Even_SNPS_PD[r,11] <- 0
  } else {
    Even_SNPS_PD[r,11] <- 1
  }
  if (Even_SNPS_PD[r,3] > 73){
    Even_SNPS_PD[r,12] <- 0
  } else {
    Even_SNPS_PD[r,12] <- 1
  }
}

for (r in 1:dim(Odd_SNPS_PD)[1]){
  if (Odd_SNPS_PD[r,3] < 17){
    Odd_SNPS_PD[r,11] <- 0
  } else {
    Odd_SNPS_PD[r,11] <- 1
  }
  if (Odd_SNPS_PD[r,3] > 73){
    Odd_SNPS_PD[r,12] <- 0
  } else {
    Odd_SNPS_PD[r,12] <- 1
  }
}
head(Even_SNPS_PD)
dim(Even_SNPS_PD)
head(Odd_SNPS_PD)
dim(Odd_SNPS_PD)


####Combine those flags for each SNP into flags for each tag, This is for Garrett's primer design pipeline

##EVEN
head(Even_SNPS_PD)
Even_tags_PD <- unique(Even_SNPS_PD$Tag)
head(Even_tags_PD)
length(Even_tags_PD)

#create the data frame to put the results into 
Even_SNPs_InRange_PD <-data.frame(Even_tags_PD)
Even_SNPs_InRange_PD$LeftLimit <-NA
Even_SNPs_InRange_PD$RightLimit <-NA
head(Even_SNPs_InRange_PD)

for (s in 1:length(Even_tags_PD)) {
  tested_tag <- Even_tags_PD[s]
  trial_set <- Even_SNPS_PD[which(Even_SNPS_PD$Tag %in% tested_tag),]
  if (any(trial_set$LeftLimit == 0)){
    Even_SNPs_InRange_PD$LeftLimit[s]<- 0
  } else {
    Even_SNPs_InRange_PD$LeftLimit[s]<- 1
  }
  if (any(trial_set$RightLimit == 0)){
    Even_SNPs_InRange_PD$RightLimit[s]<- 0
  } else {
    Even_SNPs_InRange_PD$RightLimit[s]<- 1
  }
}

head(Even_SNPs_InRange_PD)


#### Write the Even even flags for the primer pipeline:
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Even_SNPs_InRange_PD.txt", "wb")
# write.table(Even_SNPs_InRange_PD,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)

##EVEN
head(Odd_SNPS_PD)
Odd_tags_PD <- unique(Odd_SNPS_PD$Tag)
head(Odd_tags_PD)
length(Odd_tags_PD)

#create the data frame to put the results into 
Odd_SNPs_InRange_PD <-data.frame(Odd_tags_PD)
Odd_SNPs_InRange_PD$LeftLimit <-NA
Odd_SNPs_InRange_PD$RightLimit <-NA
head(Odd_SNPs_InRange_PD)

for (s in 1:length(Odd_tags_PD)) {
  tested_tag <- Odd_tags_PD[s]
  trial_set <- Odd_SNPS_PD[which(Odd_SNPS_PD$Tag %in% tested_tag),]
  if (any(trial_set$LeftLimit == 0)){
    Odd_SNPs_InRange_PD$LeftLimit[s]<- 0
  } else {
    Odd_SNPs_InRange_PD$LeftLimit[s]<- 1
  }
  if (any(trial_set$RightLimit == 0)){
    Odd_SNPs_InRange_PD$RightLimit[s]<- 0
  } else {
    Odd_SNPs_InRange_PD$RightLimit[s]<- 1
  }
}

head(Odd_SNPs_InRange_PD)


### Write the Odd flags for the primer pipeline:
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Odd_SNPs_InRange_PD", "wb")
write.table(Odd_SNPs_InRange_PD,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)




###################################### Marker selection based on ranked FSt of loci 
## we now have a matrix for each lineage that shows whether the snps at a tag are within the range of primer design
## next we want to look at the FST of the snps, and choose 
#FST_file_raw <- FST_file
#FST_file <- FST_file_raw

names(FST_file)
dim(FST_file)
head(FST_file)

Even_FST_file <- FST_file
Odd_FST_file <- FST_file

### Add the SNP_inRange designation for the EVEn and ODD lineages to the FST_file 
# split the LOCUS column so that we have a column of Tags to merge with 
#Even
Even_newColNames <- c("Tag", "Pos")
Even_newCols <- colsplit(Even_FST_file$Locus, "_", Even_newColNames)
Even_FST_file <- cbind(Even_newCols, Even_FST_file)
Even_FST_file[,3] <- NULL
head(Even_FST_file)
dim(Even_FST_file)


# split the LOCUS column so that we have a column of Tags to merge with 
#Odd
Odd_newColNames <- c("Tag", "Pos")
Odd_newCols <- colsplit(Odd_FST_file$Locus, "_", Odd_newColNames)
Odd_FST_file <- cbind(Odd_newCols, Odd_FST_file)
Odd_FST_file[,3] <- NULL
head(Odd_FST_file)
dim(Odd_FST_file)

dim(Even_Tag_Flags)
head(Even_Tag_Flags)

#these are the tags that are in the haplotype file that are not in the FST file, that is OK- they were likely filtered out in the ONE SNP per tag data set 
bb <- setdiff(Even_Tag_Flags$Tag, Even_FST_file$Tag)
length(bb)

cc <- setdiff(Odd_Tag_Flags$Tag, Odd_FST_file$Tag)
length(cc)

#### This part is unecessary:
# ### Subset the FST files for the tags that we have in the Even and Odd tag flag files 
# Even_FST_file <- Even_FST_file[which(Even_Tag_Flags$Tag %in% Even_FST_file$Tag),]
# head(Even_FST_file)
# dim(Even_FST_file)
# 
# Odd_FST_file <- Odd_FST_file[which(Odd_Tag_Flags$Tag %in% Odd_FST_file$Tag),]
# head(Odd_FST_file)
# dim(Odd_FST_file)

# merge the SNP_inRange designations
Even_FST_file <- merge(Even_FST_file, Even_Tag_Flags, by = "Tag")
head(Even_FST_file)
dim(Even_FST_file)

Odd_FST_file <- merge(Odd_FST_file, Odd_Tag_Flags, by = "Tag")
head(Odd_FST_file)
dim(Odd_FST_file)


##########################
##### Create data frames for each of the lineages that has the FST and the SNPs_inRange flag
##Find the top FST for the NA_even
FST_NA_even<- Even_FST_file[,c("Locus","NA_even","Even_SNPs_InRange")]
head(FST_NA_even)
FST_file_NAeven_sort <- FST_NA_even[order(FST_NA_even$NA_even, decreasing=TRUE),]
head(FST_file_NAeven_sort)

##Find the top FST for the NA_even that have SNPS_inRange
FST_file_NAeven_sort_inRange <- FST_file_NAeven_sort[which(FST_file_NAeven_sort$Even_SNPs_InRange == "TRUE"),]
head(FST_file_NAeven_sort_inRange)
dim(FST_file_NAeven_sort_inRange)

##Find the top FST for the Asia_even
FST_ASIA_even<- Even_FST_file[,c("Locus","ASIA_even","Even_SNPs_InRange")]
head(FST_ASIA_even)
FST_file_ASIAeven_sort <- FST_ASIA_even[order(FST_ASIA_even$ASIA_even, decreasing=TRUE),]
head(FST_file_ASIAeven_sort)

##Find the top FST for the Asia_even that have SNPS_inRange
FST_file_ASIAeven_sort_inRange <- FST_file_ASIAeven_sort[which(FST_file_ASIAeven_sort$Even_SNPs_InRange == "TRUE"),]
head(FST_file_ASIAeven_sort_inRange)
dim(FST_file_ASIAeven_sort_inRange)

##Find the top FST for the NA_odd
FST_NA_odd<- Odd_FST_file[,c("Locus","NA_odd","Odd_SNPs_InRange")]
head(FST_NA_odd)
FST_file_NAodd_sort <- FST_NA_odd[order(FST_NA_odd$NA_odd, decreasing=TRUE),]
head(FST_file_NAodd_sort)

##Find the top FST for the NA_odd that have SNPS_inRange
FST_file_NAodd_sort_inRange <- FST_file_NAodd_sort[which(FST_file_NAodd_sort$Odd_SNPs_InRange == "TRUE"),]
head(FST_file_NAodd_sort_inRange)
dim(FST_file_NAodd_sort_inRange)

##Find the top FST for the Asia_odd
FST_ASIA_odd<- Odd_FST_file[,c("Locus","ASIA_odd","Odd_SNPs_InRange")]
head(FST_ASIA_odd)
FST_file_ASIAodd_sort <- FST_ASIA_odd[order(FST_ASIA_odd$ASIA_odd, decreasing=TRUE),]
head(FST_file_ASIAodd_sort)

##Find the top FST for the Asia_odd that have SNPS_inRange
FST_file_ASIAodd_sort_inRange <- FST_file_ASIAodd_sort[which(FST_file_ASIAodd_sort$Odd_SNPs_InRange == "TRUE"),]
head(FST_file_ASIAodd_sort_inRange)
dim(FST_file_NAodd_sort_inRange)

######################################
################# EVEN 
#top FST loci for NA_even, 
#change the variable to the number you want 
NAeven_var <- 500 #<--------------------------------------------dictates the number of loci
NAeven_set <- FST_file_NAeven_sort_inRange[1:NAeven_var,]
NAeven_set_avgFST <- mean(NAeven_set$NA_even)
NAeven_set_maxFST <- max(NAeven_set$NA_even)
NAeven_set_minFST <- min(NAeven_set$NA_even)
cat("NAeven FST max, min, average: ",NAeven_set_maxFST, NAeven_set_minFST, NAeven_set_avgFST )

#top FST loci for ASIAeven, 
#change the variable to the number you want 
ASIAeven_var <- 500 #<--------------------------------------------dictates the number of loci
ASIAeven_set <- FST_file_ASIAeven_sort_inRange[1:ASIAeven_var,]
dim(ASIAeven_set)
ASIAeven_set_avgFST <- mean(ASIAeven_set$ASIA_even)
ASIAeven_set_maxFST <- max(ASIAeven_set$ASIA_even)
ASIAeven_set_minFST <- min(ASIAeven_set$ASIA_even)
cat("ASIAeven FST max, min, average: ", ASIAeven_set_maxFST, ASIAeven_set_minFST, ASIAeven_set_avgFST )

###set operations to see how many overlap: 
NA_ASIA_even_union <- union(NAeven_set$Locus, ASIAeven_set$Locus)
length(NA_ASIA_even_union)

##In Asia Even, not in NA Even
In_ASIA_even_Not_NA <- setdiff(ASIAeven_set$Locus, NAeven_set$Locus)
length(In_ASIA_even_Not_NA)

Intersection_NA_even_ASIA <- intersect(NAeven_set$Locus, ASIAeven_set$Locus)
length(Intersection_NA_even_ASIA)

## make a new dataframe that has the FST info for the union of the sets 
NA_ASIA_even_union_FST <- FST_file[FST_file$Locus %in% NA_ASIA_even_union,]  
NA_ASIA_even_union_FST_sort <-NA_ASIA_even_union_FST[order(NA_ASIA_even_union_FST$EVEN, decreasing = TRUE),]
head( NA_ASIA_even_union_FST_sort)
dim(NA_ASIA_even_union_FST_sort)

# #### Write a list of the loci that were found in the overlapping set for Evens:
# ### Change the name to reflect the number choice
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/ASIAeven_set", "wb")
# write.table(ASIAeven_set,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)

################# ODD
#top FST loci for NA_odd,
#change the variable to the number you want 
NAodd_var <- 500 #<--------------------------------------------dictates the number of loci
NAodd_set <- FST_file_NAodd_sort_inRange[1:NAodd_var,]
NAodd_set_avgFST <- mean(NAodd_set$NA_odd)
NAodd_set_maxFST <- max(NAodd_set$NA_odd)
NAodd_set_minFST <- min(NAodd_set$NA_odd)
cat("NAodd FST max, min, average: ", NAodd_set_maxFST, NAodd_set_minFST, NAodd_set_avgFST)

#top FST loci for ASIAodd, 
#change the variable to the number you want 
ASIAodd_var <- 500 #<--------------------------------------------dictates the number of loci
ASIAodd_set <- FST_file_ASIAodd_sort_inRange[1:ASIAodd_var,]
ASIAodd_set_avgFST <- mean(ASIAodd_set$ASIA_odd)
ASIAodd_set_maxFST <- max(ASIAodd_set$ASIA_odd)
ASIAodd_set_minFST <- min(ASIAodd_set$ASIA_odd)
cat("ASIAodd FST max, min, average: ", ASIAodd_set_maxFST, ASIAodd_set_minFST, ASIAodd_set_avgFST )

###set operations to see how many overlap: 
NA_ASIA_odd_union <- union(NAodd_set$Locus, ASIAodd_set$Locus)
length(NA_ASIA_odd_union)

##In Asia Odd, not in NA Odd
In_ASIA_ODD_Not_NA <- setdiff(ASIAodd_set$Locus, NAodd_set$Locus)
length(In_ASIA_ODD_Not_NA)

##In NA Odd, not in Asia Odd
In_NA_ODD_Not_ASIA <- setdiff(NAodd_set$Locus, ASIAodd_set$Locus)
length(In_NA_ODD_Not_ASIA)

Intersection_NA_odd_ASIA <- intersect(NAodd_set$Locus, ASIAodd_set$Locus)
length(Intersection_NA_odd_ASIA)

## make a new dataframe that has the FST info for the overlapping set 
NA_ASIA_odd_union_FST <- FST_file[FST_file$Locus %in% NA_ASIA_odd_union,]  
NA_ASIA_odd_union_FST_sort <-NA_ASIA_odd_union_FST[order(NA_ASIA_odd_union_FST$ODD, decreasing = TRUE),]
head( NA_ASIA_odd_union_FST_sort)
dim(NA_ASIA_odd_union_FST_sort)

# #### Write a list of the loci that were found in the overlapping set for Odds:
# ### Change the name to reflect the number choice
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/ASIAodd_set.txt", "wb")
# write.table(ASIAodd_set,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)

#####################################################
#########overlap between two panels
overlap <- intersect(NA_ASIA_even_union, NA_ASIA_odd_union)
length(overlap)
head(overlap)

## make a new dataframe that has the FST info for the overlapping set 
Lineage_overlap_FST <- FST_file[FST_file$Locus %in% overlap,]  
Lineage_overlap_FST_sort <-Lineage_overlap_FST[order(Lineage_overlap_FST$ALL_noSusit, decreasing = TRUE),]
head(Lineage_overlap_FST_sort)
dim(Lineage_overlap_FST_sort)

# ####Write a list of the loci that were found to be overlapping between the even and the odd panels:
# ### Change the name of the file to reflect the choice of loci and the final number
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/OverlapBetweenEvenOdd500_126.txt", "wb")
# write.table(Lineage_overlap_FST_sort,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)

##############################
########## top X FST loci that were excluded from the combined analysis because of the SNPS_inRange flag
FST_file_NAeven_sort

##pull out the Failed in the SNPS_inRange flag
FST_NAeven_failed <- FST_file_NAeven_sort[which(FST_file_NAeven_sort$Even_SNPs_InRange == "FALSE"),] 
FST_ASIAeven_failed <- FST_file_ASIAeven_sort[which(FST_file_ASIAeven_sort$Even_SNPs_InRange == "FALSE"),]
FST_ASIAodd_failed <- FST_file_ASIAodd_sort[which(FST_file_ASIAodd_sort$Odd_SNPs_InRange == "FALSE"),]
FST_NAodd_failed <- FST_file_NAodd_sort[which(FST_file_NAodd_sort$Odd_SNPs_InRange == "FALSE"),]

##Keep X number of loci that have the highest FST #<--------------------------------------------dictates the number of loci
NAeven_top_failed <- FST_NAeven_failed[1:250,]
ASIAeven_top_failed <- FST_ASIAeven_failed[1:250,]
NAodd_top_failed <- FST_NAodd_failed[1:250,]
ASIAodd_top_failed <- FST_ASIAodd_failed[1:250,]

head(NAeven_top_failed)
head(ASIAeven_top_failed)
head(NAodd_top_failed)
head(ASIAodd_top_failed)

### Export the odd and even lists of top failed loci
Even_top_failed_tags <- unique(append(NAeven_top_failed$Locus,ASIAeven_top_failed$Locus))
length(Even_top_failed_tags)

Odd_top_failed_tags <- unique(append(NAodd_top_failed$Locus,ASIAodd_top_failed$Locus))
length(Odd_top_failed_tags)

# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Odd_top_failed_tags.txt", "wb")
# write.table(Odd_top_failed_tags,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
# close(outputFile)
# 
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Even_top_failed_tags.txt", "wb")
# write.table(Even_top_failed_tags,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
# close(outputFile)


################## Make the Master_FST 
#FST_file_raw <- FST_file
#FST_file <- FST_file_raw

head(FST_file)
dim(FST_file)
head(SNP_positions_all)
dim(SNP_positions_all)

newColNames <- c("Tag", "Pos")
newCols <- colsplit(FST_file$Locus, "_", newColNames)
FST_file <- cbind(newCols, FST_file)
FST_file[,3] <- NULL
head(FST_file)
dim(FST_file)
#FST_file <- within(FST_file, rm("Pos","Index"))
head(FST_file)


##merge the FST_file the SNP_positions_all to see what the SNP positions are for the loci that failed. 
masterFST_file <- merge(SNP_positions_all, FST_file, by= "Tag")
masterFST_file <- rename(masterFST_file, Locus = Locus.x, FST_Locus = Locus.y)
head(masterFST_file)

#trim columns we dont need
#masterFST_file$Pos <-NULL
masterFST_file$SNP_Index <-NULL
masterFST_file$SNP_tag_index <-NULL
head(masterFST_file)

###Pull out the positions of the SNPS in the top x number of loci with the highest FST that failed the SNP postion flag
FST_NAeven_top_failed_snps <- masterFST_file[masterFST_file$FST_Locus %in% NAeven_top_failed$Locus,]
FST_NAeven_top_failed_snps_avg <- FST_NAeven_top_failed_snps[FST_NAeven_top_failed_snps$Locus %in% unique(FST_NAeven_top_failed_snps$FST_Locus),]
dim(FST_NAeven_top_failed_snps_avg)
head(FST_NAeven_top_failed_snps_avg)

FST_ASIAeven_top_failed_snps <- masterFST_file[masterFST_file$FST_Locus %in% ASIAeven_top_failed$Locus,]
FST_ASIAeven_top_failed_snps_avg <- FST_ASIAeven_top_failed_snps[FST_ASIAeven_top_failed_snps$Locus %in% unique(FST_ASIAeven_top_failed_snps$FST_Locus),]
dim(FST_ASIAeven_top_failed_snps_avg)
head(FST_ASIAeven_top_failed_snps_avg)

FST_ASIAodd_top_failed_snps <- masterFST_file[masterFST_file$FST_Locus %in% ASIAodd_top_failed$Locus,]
FST_ASIAodd_top_failed_snps_avg <- FST_ASIAodd_top_failed_snps[FST_ASIAodd_top_failed_snps$Locus %in% unique(FST_ASIAodd_top_failed_snps$FST_Locus),]
dim(FST_ASIAodd_top_failed_snps_avg)
head(FST_ASIAodd_top_failed_snps_avg)

FST_NAodd_top_failed_snps <- masterFST_file[masterFST_file$FST_Locus %in% NAodd_top_failed$Locus,]
FST_NAodd_top_failed_snps_avg <- FST_NAodd_top_failed_snps[FST_NAodd_top_failed_snps$Locus %in% unique(FST_NAodd_top_failed_snps$FST_Locus),]
dim(FST_NAodd_top_failed_snps_avg)
head(FST_NAodd_top_failed_snps_avg)

FST_NAeven_top_failed_snps_avgFST <- mean(FST_NAeven_top_failed_snps_avg$NA_even)
FST_NAeven_top_failed_snps_maxFST <- max(FST_NAeven_top_failed_snps_avg$NA_even)
FST_NAeven_top_failed_snps_minFST <- min(FST_NAeven_top_failed_snps_avg$NA_even)
cat("NAeven FST max, min, average: ", FST_NAeven_top_failed_snps_maxFST, FST_NAeven_top_failed_snps_minFST,
    FST_NAeven_top_failed_snps_avgFST )
### Change the name of the file to reflect the choice of loci
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/FST_NAeven_top250_failed_snps.txt", "wb")
write.table(FST_NAeven_top_failed_snps_avg,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

head(FST_ASIAeven_top_failed_snps_avg)
dim(FST_ASIAeven_top_failed_snps_avg)
FST_ASIAeven_top_failed_snps_avgFST <- mean(FST_ASIAeven_top_failed_snps_avg$ASIA_even)
FST_ASIAeven_top_failed_snps_maxFST <- max(FST_ASIAeven_top_failed_snps_avg$ASIA_even)
FST_ASIAeven_top_failed_snps_minFST <- min(FST_ASIAeven_top_failed_snps_avg$ASIA_even)
cat("ASIAeven FST max, min, average: ", FST_ASIAeven_top_failed_snps_maxFST, FST_ASIAeven_top_failed_snps_minFST,
    FST_ASIAeven_top_failed_snps_avgFST )
## Change the name of the file to reflect the choice of loci
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/FST_ASIAeven_top250_failed_snps.txt", "wb")
write.table(FST_ASIAeven_top_failed_snps_avg,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

head(FST_NAodd_top_failed_snps_avg)
dim(FST_NAodd_top_failed_snps_avg)
FST_NAodd_top_failed_snps_avgFST <- mean(FST_NAodd_top_failed_snps_avg$NA_odd)
FST_NAodd_top_failed_snps_maxFST <- max(FST_NAodd_top_failed_snps_avg$NA_odd)
FST_NAodd_top_failed_snps_minFST <- min(FST_NAodd_top_failed_snps_avg$NA_odd)
cat("NAodd FST max, min, average: ", FST_NAodd_top_failed_snps_maxFST, FST_NAodd_top_failed_snps_minFST,
    FST_NAodd_top_failed_snps_avgFST )
### Change the name of the file to reflect the choice of loci
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/FST_NAodd_top250_failed_snps.txt", "wb")
write.table(FST_NAodd_top_failed_snps_avg,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

head(FST_ASIAodd_top_failed_snps_avg)
dim(FST_ASIAodd_top_failed_snps_avg)
FST_ASIAodd_top_failed_snps_avgFST <- mean(FST_ASIAodd_top_failed_snps_avg$ASIA_odd)
FST_ASIAodd_top_failed_snps_maxFST <- max(FST_ASIAodd_top_failed_snps_avg$ASIA_odd)
FST_ASIAodd_top_failed_snps_minFST <- min(FST_ASIAodd_top_failed_snps_avg$ASIA_odd)
cat("ASIAodd FST max, min, average: ", FST_ASIAodd_top_failed_snps_maxFST, FST_ASIAodd_top_failed_snps_minFST,
    FST_ASIAodd_top_failed_snps_avgFST )
### Change the name of the file to reflect the choice of loci
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/FST_ASIAodd_top250_failed_snps.txt", "wb")
write.table(FST_ASIAodd_top_failed_snps_avg,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
close(outputFile)

#########################################################################
#### Use the combination of the 250 top failed from each of the groups for the lineage. So take the top 250 from each group and then do a union
####### overlap of the Even NA and ASian failed 
Even_failed_union <- union(NAeven_top_failed$Locus, ASIAeven_top_failed$Locus)
length(Even_failed_union)
head(Even_failed_union)

##In Asia Even, not in NA Even
In_ASIA_Even_Not_NA_failed <- setdiff(ASIAeven_top_failed$Locus, NAeven_top_failed$Locus)
length(In_ASIA_Even_Not_NA_failed)

##In NA Odd, not in Asia Odd
In_NA_Even_Not_ASIA_failed <- setdiff(NAeven_top_failed$Locus, ASIAeven_top_failed$Locus)
length(In_NA_Even_Not_ASIA_failed)

Intersection_NA_Even_ASIA_failed <- intersect(NAeven_top_failed$Locus, ASIAeven_top_failed$Locus)
length(Intersection_NA_Even_ASIA_failed)


####### overlap of the Odd NA and ASian failed 
Odd_failed_union<- union(NAodd_top_failed$Locus, ASIAodd_top_failed$Locus)
length(Odd_failed_union)
head(Odd_failed_union)

##In Asia Odd, not in NA Odd
In_ASIA_ODD_Not_NA_failed <- setdiff(ASIAodd_top_failed$Locus, NAodd_top_failed$Locus)
length(In_ASIA_ODD_Not_NA_failed)

##In NA Odd, not in Asia Odd
In_NA_ODD_Not_ASIA_failed <- setdiff(NAodd_top_failed$Locus, ASIAodd_top_failed$Locus)
length(In_NA_ODD_Not_ASIA_failed)

Intersection_NA_odd_ASIA_failed <- intersect(NAodd_top_failed$Locus, ASIAodd_top_failed$Locus)
length(Intersection_NA_odd_ASIA_failed)

## Add these list of failed to the list of passed to get a list to pull the sequences to make a FASTA file 
Even_failed_passed_list <- append(NA_ASIA_even_union_FST_sort$Locus, Even_failed_union)
head(Even_failed_passed_list)
length(Even_failed_passed_list)
# ## Change the name of the file
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Even_failed_passed_list.txt", "wb")
# write.table(Even_failed_passed_list,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
# close(outputFile)

Odd_failed_passed_list <- append(NA_ASIA_odd_union_FST_sort$Locus, Odd_failed_union)
head(Odd_failed_passed_list)
length(Odd_failed_passed_list)
#  ## Change the name of the file
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Odd_failed_passed_list.txt", "wb")
# write.table(Odd_failed_passed_list,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
# close(outputFile)

###########Make a table that has the FST, Locus_FST, SNP Pos and SNP Flag for all the pass and failed loci tog. 
#combine the two failed lists, so both lineages in one 
Both_failed_union <- union(Odd_failed_union, Even_failed_union)
length(Both_failed_union)
head(Both_failed_union)

#combine ALL THE LISTS, failed and passed, so both lineages in one 
masterCombined_BothLineages <- union(Even_failed_passed_list, Odd_failed_passed_list) #<- THIS IS THE LIST OF ALL THE LOCI TO MAKE FASTA FILE 
length(masterCombined_BothLineages)
head(masterCombined_BothLineages)
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/masterCombined_BothLineages_2312.txt", "wb")
# write.table(masterCombined_BothLineages,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
# close(outputFile)

#Make the COMBINED loci list table 
FST_file_COMBINED_table <- FST_file[FST_file$Locus %in% masterCombined_BothLineages,] 
head(FST_file_COMBINED_table)
dim(FST_file_COMBINED_table)

# ## Change the name of the file
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/FST_file_COMBINED_table_2312.txt", "wb")
# write.table(FST_file_COMBINED_table,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)

#Make the failed loci their own table that has the many SNPS that they have. 
masterFST_file_failed_union <- masterFST_file[masterFST_file$FST_Locus %in% Both_failed_union,] 
head(masterFST_file_failed_union)
dim(masterFST_file_failed_union)

# ## Change the name of the file
# outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/masterFST_file_failed_union.txt", "wb")
# write.table(masterFST_file_failed_union,outputFile,quote=FALSE,row.names=FALSE,col.names=TRUE,eol="\n")
# close(outputFile)
