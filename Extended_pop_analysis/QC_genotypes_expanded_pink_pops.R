### Genotype QC on expanded pink populations 
###   Testing the loci genotyped the same between previous study and current efforts.
###   QC individual based PCA plots for each population
###  Carolyn Tarpey | Nov 2017
### ---------------------------------------

#This is R code requires a genepop file for each of the genotype sets and should be stripped of the header line that Stacks puts in.  
#The second line should start with a tab, the loci names should be tab deliminated, and there should be no "pop" between the samples. 
#This formatting ensures R imports the genepop file as a table. 

#use the R code initialPinkFiltering_edited.R to figure out what the loci overlap is between the old and new data sets. 
#Use that to subset the old and the new genepop files. 

#Pink data filtering
library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)


#####################Test to see that what we got in our whitelist has the 16681 loci that we want

pink16681 <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/pink16681.txt")
pink16681 <- unlist(pink16681, use.names = FALSE)
length(pink16681)
head(pink16681)

locipassedMAF_temp<-locipassedMAF_temp[,1]
head(locipassedMAF_temp)
length(locipassedMAF_temp)
#locipassedMAF_temp<-locipassedMAF #to re-set
intersect_filteredWhitelist_16681 <- intersect(locipassedMAF_temp, pink16681)
diff_filteredWhitelist_16681 <-setdiff(pink16681,locipassedMAF_temp)
diff_16681_filteredWhitelist <-setdiff(locipassedMAF_temp, pink16681)

length(intersect_filteredWhitelist_16681)
head(locipassedMAF_temp)

length(diff_filteredWhitelist_16681)
head(diff_filteredWhitelist_16681)

length(diff_16681_filteredWhitelist)
head(diff_16681_filteredWhitelist)

#########################################

#output file to use to reduce the original genepop file and the new genepop, this list is the intersect of the two. 

#output list of filtered loci, need to open as binary file (using write binary, "wb") to give unix line endings
#format dataset to remove X from locus names
intersect_filteredWhitelist_16681_IDS<-gsub("X","",intersect_filteredWhitelist_16681)
head(intersect_filteredWhitelist_16681_IDS)
outputFile<-file("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/intersect_filteredWhitelist_16681_IDS", "wb")
write.table(intersect_filteredWhitelist_16681_IDS,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
close(outputFile)

#######################################

#import the two new genepop files 
#make sure they have been formatted, missing their title line, tab separated loci, 
#tab at the start of the loci line and the same number of inds
OLD_genepop<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/QC_test_16681_new_genotypes/OLD_genepop_16662.txt",colClasses="character")
OLD_genepop[1:5,1:5]
dim(OLD_genepop)
#head(OLD_genepop)

NEW_genepop<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/QC_test_16681_new_genotypes/NEW_genepop_16662.txt",colClasses="character")
NEW_genepop[1:5,1:5]
dim(NEW_genepop)

# the old genepop file has been filtered for individuals, and IDK which so figure that our below and remove them from the new file
OLD_inds<-rownames(OLD_genepop)
length(OLD_inds)
NEW_inds<-rownames(NEW_genepop)
length(NEW_inds)

# rows2remove<-setdiff(NEW_inds, OLD_inds)
# length(rows2remove)
# head(rows2remove)

rows2keep<-intersect(NEW_inds, OLD_inds)
length(rows2keep)
head(rows2keep)

filtered_inds <- setdiff(NEW_inds, OLD_inds)

NEW_genepop_ed<-subset(NEW_genepop[rows2keep,])
dim(NEW_genepop_ed)
NEW_genepop_ed[1:5,1:5]

#####################Compare the matrices

z<-NEW_genepop_ed==OLD_genepop

z <- NEW_genepop_ed==OLD_genepop # avoid to later compute this twice
#testComparemat<- ifelse(z, 0, OLD_genepop) # get the desired matrix

Percent_diff_MAt<- round(sum(!z)/length(z)*100, 2) # get the percentage of non equal values
Percent_diff_MAt<- sum(!z)/length(z)*100 # get the percentage of non equal values unrounded
Percent_diff_MAt

# row by row, col by col

percentDif <- function(x)   (round((sum(!x)/length(x)*100),2))
#percentDif <- function(x)   (sum(!x)/length(x)*100)
#percentDif <- function(x)   sum(!x)

MatCompareRow<- apply(z,1,percentDif)
which.max(MatCompareRow)
max(MatCompareRow)

MatCompareCol<- apply(z,2,percentDif)
which.max(MatCompareCol)
max(MatCompareCol)

#########plot the percent difference by locus
#plot histogram of SNP position
#histogram won't work because it's discrete
MatCompareCol_df<-data.frame(keynames=names(MatCompareCol), value=MatCompareCol, row.names = NULL)
colnames(MatCompareCol_df)<-c("Locus", "PerCDiff")

LociPerC<-as.data.frame(table(MatCompareCol_df$PerCDiff))
colnames(LociPerC)<-c("Percent_Difference","Count_of_Loci")
LociPerC_trim<-LociPerC[-1,]

ggplot()+geom_point(data=LociPerC_trim,aes(y=Count_of_Loci,x=Percent_Difference))+xlab("Percent Difference")+ylab("Count of Loci")+theme_bw()+ggtitle("Percent Difference in Loci Genotypes Between Old and New Datasets")
ggsave(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Filtering/Percent Difference.pdf")
#plot barchart of genotype rate for loci
ggplot()+geom_bar(data=LociPerC,aes(x=Percent_Difference))+ggtitle("Count of Loci per Genotype Rate ")
