### To select markers for two pink Amplicon Panels, one per lineage
###   Using Ranked FST of one snp per tag loci in a hierarchical structure
###    Identify which SNPs are real variation in each lineage using Haplotype file 
###   Carolyn Tarpey | May 2018
### ---------------------------------------

#This code is designed to help make Genepop files of the loci that are in the panels for each of the lineages
#The loci have not yet had primers made, so they are the 1197 - low complexity and 14661 - low complexity for odd and even respectively.
#The genepop file we are using to trim down to size is the general batch_4.genepop, which was created with the 31485 whitelist.  

###This code also requires a population map, with one column that has each individual and and another column for their population number. 
#These were made one for each lineage in R. in the script Pink_prep_PLINK_LD_popPlots.R
### The lakelse population individuals are named incorrectly, the lineage/year for each of the individuals is switched, 
##The correct population/lineage is named in the pop map


library(vcfR)
library(stringr)
library(dplyr)
library(gdata)
library(hierfstat)
library(adegenet)
library(argparse)
library(stringi)
library(reshape2)

#______________________________________________________________________________________________
################################IMPORT The population maps These have to have a comma after each name, to match the genepop names. 

EVEN_pops <-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Genepop/EvenPanelMarkers/EVEN_inds_R.txt", header=FALSE, 
                     stringsAsFactors = FALSE)
names(EVEN_pops) <- c("POPnum", "Indv_name")
head(EVEN_pops)
dim(EVEN_pops)

ODD_pops <-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Genepop/OddPanelMarkers/Odd_inds_R.txt", header=FALSE, 
                       stringsAsFactors = FALSE )
names(ODD_pops) <- c("POPnum", "Indv_name")
head(ODD_pops)
dim(ODD_pops)

#add comma to name if not already there!
EVEN_pops$Indv_name <- paste0(EVEN_pops$Indv_name, ",")
ODD_pops$Indv_name <- paste0(ODD_pops$Indv_name, ",")

#####################Import the two genepop files to be edited-  these not one line per snp, but they have been extensively edited. 
#They have had the "pop" removed, the , between each snp name replaced with a tab, a tab to start that snp name line, the header line deleted, and the , in the pop names deleted. 
#The col as factor is to keep the 4 digit format, 0303, and 0000 etc

EVEN_r_genepop <-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Genepop/EvenPanelMarkers/batch_4edited_forR.genepop", header=TRUE, 
                        colClasses="factor", na.strings = "0000" )
EVEN_r_genepop[1:5,1:5]
dim(EVEN_r_genepop)

ODD_r_genepop <-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Genepop/OddPanelMarkers/batch_4edited_forR.genepop", header=TRUE, 
                           colClasses="factor", na.strings = "0000" )
ODD_r_genepop[1:5,1:5]
dim(ODD_r_genepop)

################################IMPORT the list of the loci in each panel

EVEN_paneltags <-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Even_failed_passed_list.txt", header=FALSE, 
                       stringsAsFactors = FALSE)
dim(EVEN_paneltags)

ODD_paneltags <-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Odd_failed_passed_list.txt", header=FALSE, 
                      stringsAsFactors = FALSE )
dim(ODD_paneltags)


#______________________________________________________________________________________________
##################### #SPlit the locus from each panel into two columns, tags and SNPS

EVEN_panelSNPS <- data.frame(str_split_fixed(EVEN_paneltags$V1,"_",2))
names(EVEN_panelSNPS) <-c("Tag","SNP")
head(EVEN_panelSNPS)

ODD_panelSNPS <- data.frame(str_split_fixed(ODD_paneltags$V1,"_",2))
names(ODD_panelSNPS) <-c("Tag","SNP")
head(ODD_panelSNPS)

#Get a list of the tags for each of the panels. 
EVEN_paneltags<- EVEN_panelSNPS$Tag
length(EVEN_paneltags)
EVEN_paneltags<- sort(as.numeric(as.character(EVEN_paneltags)), decreasing= FALSE)

ODD_paneltags <- ODD_panelSNPS$Tag
length(ODD_paneltags)
ODD_paneltags<- sort(as.numeric(as.character(ODD_paneltags)), decreasing= FALSE)

#Remove the X that gets added by R for some reason
names(EVEN_r_genepop) <- gsub("X","",names(EVEN_r_genepop))
names(ODD_r_genepop) <- gsub("X","",names(ODD_r_genepop))

#get a list of each of the snps in the OG genepop files, to be edited down to the few that we want.
Even_genepop_names <- names(EVEN_r_genepop)
Even_genepop_names_split <- data.frame(str_split_fixed(Even_genepop_names,"_",2))
names(Even_genepop_names_split) <-c("Tag","SNP")
Even_genepop_names_split$Locus <- Even_genepop_names
head(Even_genepop_names_split)

Odd_genepop_names <- names(ODD_r_genepop)
Odd_genepop_names_split <- data.frame(str_split_fixed(Odd_genepop_names,"_",2))
names(Odd_genepop_names_split) <-c("Tag","SNP")
Odd_genepop_names_split$Locus <- Odd_genepop_names
head(Odd_genepop_names_split)

Even_genepop_keep_loci <- Even_genepop_names_split[ Even_genepop_names_split$Tag%in% EVEN_paneltags, ]
head(Even_genepop_keep_loci)
dim(Even_genepop_keep_loci)

Odd_genepop_keep_loci <- Odd_genepop_names_split[ Odd_genepop_names_split$Tag %in% ODD_paneltags, ]
head(Odd_genepop_keep_loci)
dim(Odd_genepop_keep_loci)

##Filter down the raw genepop files so that they only have the columns of the snps that we want to keep
Even_filt_snps_genepop <- EVEN_r_genepop[, names(EVEN_r_genepop) %in% Even_genepop_keep_loci$Locus ]
Even_filt_snps_genepop[1:5,1:5]
dim(Even_filt_snps_genepop)

Odd_filt_snps_genepop <- ODD_r_genepop[, names(ODD_r_genepop) %in% Odd_genepop_keep_loci$Locus ]
Odd_filt_snps_genepop[1:5,1:5]
dim(Odd_filt_snps_genepop)

##Filter down the genotype files that have just the snps we want for each panel, to have only the individuals in each panel
Even_filt_snps_inds_genepop <- Even_filt_snps_genepop[ rownames(Even_filt_snps_genepop)  %in% EVEN_pops$Indv_name ,]
Even_filt_snps_inds_genepop[1:5,1:5]
dim(Even_filt_snps_inds_genepop)

Odd_filt_snps_inds_genepop <- Odd_filt_snps_genepop[ rownames(Odd_filt_snps_genepop) %in% ODD_pops$Indv_name ,]
Odd_filt_snps_inds_genepop[1:5,1:5]
dim(Odd_filt_snps_inds_genepop)

###THIS DOESNT WORK> AT THIS POINT< EXPORT THE ABOVE AND FIND REPLACE NA WITH 0000 in notebook++
#Replace the NA that has been coerced with "0000' which is the GENEPOP coding for NA
# Even_filt_snps_inds_genepop <- as.data.frame(Even_filt_snps_inds_genepop)
# Even_filt_snps_inds_genepop[is.na(Even_filt_snps_inds_genepop)] <-"0000"
# Even_filt_snps_inds_genepop[1:5,1:5]

####Write a table to a file
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Genepop/OddPanelMarkers/Odd_filt_snps_inds_genepop.txt", "wb")
write.table(Odd_filt_snps_inds_genepop,outputFile,quote=FALSE,row.names=TRUE,col.names=TRUE,eol="\n")
close(outputFile)


####Write a table to a file
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Genepop/EvenPanelMarkers/Even_filt_snps_inds_genepop.txt", "wb")
write.table(Even_filt_snps_inds_genepop,outputFile,quote=FALSE,row.names=TRUE,col.names=TRUE,eol="\n")
close(outputFile)
