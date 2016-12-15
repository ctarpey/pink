### R script for Fish 444 pop structure lab
### Charlie Waters 1/30/2014

#########################################
# 1. Set up R session
#########################################
#A. #Install each of the following packages; 
#To minimize the chance of errors, 
#hit Ctrl+Enter for each line separately
install.packages("BDgraph")
install.packages("ape")   
install.packages("ade4")
install.packages("adegenet")
install.packages("diveRsity")
install.packages("doParallel")
install.packages("foreach")
install.packages("genetics")
install.packages("hierfstat")
install.packages("iterators")
install.packages("parallel")### This may not be available but it's ok.
install.packages("sendplot")
install.packages("xlsx")

#B. Load the installed packages in to R for use
#Again, hit Ctrl+Enter for each line separately
#Or use R studio panel

library(ape)  
library(ade4)
library(adegenet)
library(diveRsity)
library(doParallel)
library(foreach)
library(genetics)
library(hierfstat)
library(iterators)
library(parallel)
library(sendplot)
library(xlsx)
library(BDgraph)

# C. Set the working directory to a specific folder 
#just change my UWNetID to yours 
#and change the folder to the one you created on the Desktop

setwd("C:/users/jheare/Desktop/researchproject") 

#D. import a Genepop file and saves it as "cod_data_genepop" 
#in the form of a "genind" object
salmon <- read.genepop("Class_data_genepop_new.gen", missing=NA)  

#E. Gives a summary of the data, 
#including number of individuals, alleles per locus, etc.
summary(salmon)  



######################################################################################################
# CODE FOR PART A - TEST WHETHER THE POPULATIONS DEVIATE FROM HARDY-WEINBERG EQUILIBRIUM
######################################################################################################

### Break up the entire data set into genind objects for each population
pop_labels <- c(rep("JA",90),rep("KO",86),
                rep("AD",45),rep("AI",92),rep("AT",45),rep("UP",87),
                rep("KI",106),rep("HS",89),rep("WA",69),rep("SG",94),rep("PS",20))  
## Creates a vector containing the population assignments of each individual
#numbers in code are numbers of individuals in each population


## Creates a list of genind objects for each population
pops_separated <- seppop(cod_data_genepop,pop=pop_labels)
names(pops_separated)


#Creates a genind object comprising only the AD individuals
data_AD <-pops_separated$AD  
#Verify that the genind object has the correct number of individuals and loci
data_AD                      
###Repeat for all populations
data_AI <-pops_separated$AI  
data_AT <-pops_separated$AT
data_HS <-pops_separated$HS
data_JA <-pops_separated$JA
data_KI <-pops_separated$KI
data_KO <-pops_separated$KO
data_PS <-pops_separated$PS
data_SG <-pops_separated$SG
data_UP <-pops_separated$UP
data_WA <-pops_separated$WA

#### Compute observed and expected heterzygosity for 
#each population over all loci
summary_AD <- summary(data_AD)
mean(summary_AD$Hexp)
mean(summary_AD$Hobs)
### Test whether Hobs and Hexp are significantly different
t.test(summary_AD$Hobs,summary_AD$Hexp,paired=TRUE)  

summary_AI <- summary(data_AI)
mean(summary_AI$Hexp)
mean(summary_AI$Hobs)
t.test(summary_AI$Hobs,summary_AI$Hexp,paired=TRUE)

summary_AT <- summary(data_AT)
mean(summary_AT$Hexp)
mean(summary_AT$Hobs)
t.test(summary_AT$Hobs,summary_AT$Hexp,paired=TRUE)

summary_HS <- summary(data_HS)
mean(summary_HS$Hexp)
mean(summary_HS$Hobs)
t.test(summary_HS$Hobs,summary_HS$Hexp,paired=TRUE)

summary_JA <- summary(data_JA)
mean(summary_JA$Hexp)
mean(summary_JA$Hobs)
t.test(summary_JA$Hobs,summary_JA$Hexp,paired=TRUE)

summary_KI <- summary(data_KI)
mean(summary_KI$Hexp)
mean(summary_KI$Hobs)
t.test(summary_JA$Hobs,summary_JA$Hexp,paired=TRUE)

summary_KO <- summary(data_KO)
mean(summary_KO$Hexp)
mean(summary_KO$Hobs)
t.test(summary_KO$Hobs,summary_KO$Hexp,paired=TRUE)

summary_PS <- summary(data_PS)
mean(summary_PS$Hexp)
mean(summary_PS$Hobs)
t.test(summary_PS$Hobs,summary_PS$Hexp,paired=TRUE)

summary_SG <- summary(data_SG)
mean(summary_SG$Hexp)
mean(summary_SG$Hobs)
t.test(summary_SG$Hobs,summary_SG$Hexp,paired=TRUE)

summary_UP <- summary(data_UP)
mean(summary_UP$Hexp)
mean(summary_UP$Hobs)
t.test(summary_UP$Hobs,summary_UP$Hexp,paired=TRUE)

summary_WA <- summary(data_WA)
mean(summary_WA$Hexp)
mean(summary_WA$Hobs)
t.test(summary_WA$Hobs,summary_WA$Hexp,paired=TRUE)


##################################################################################################################################
### QUESTION:
### ARE THERE ANY POPULATIONS THAT SIGNIFICANTLY DEVIATE FROM HWE? IF SO, WHICH ONES?
##################################################################################################################################



### Now test for deviations from HWE at each locus within each population

#This test uses simulation to compute a pvalue for HWE
HWE_test_results <- HWE.test(salmon,pop=NULL,permut=TRUE,
                             nsim=10000,res.type="matrix") 

#We want to correct for multiple tests, even though the Bonferroni is conservative
corrected_pval <- 0.05/(11*11) 
corrected_pval
#Creates a table of True/False for loci that are out of HWE (TRUE=out of HWE)
HWE_logical <- HWE_test_results<corrected_pval  
#Identify which loci are out of HWE in the various populations
HWE_logical   


#################################################################################################################################
### QUESTION: 
### IS THERE ANY EVIDENCE FOR DEVIATIONS FROM HWE AT ANY OF THE LOCI? 
### EXPORT THE HW_TEST_RESULTS TABLE, BOLD THE LOCI.
write.table(HWE_test_results,file="HWresults.csv",sep=",",row.names=F)
##################################################################################################################################


######################################################################################################
# CODE FOR PART B - TEST FOR GENETIC DIFFERENTIATION BETWEEN PAIRS OF POPULATIONS 
######################################################################################################


#####Test to see if population has a significant effect on genetic structure
#(A test for homogeneity across all populations)

###converts data for use in package "hierfstat"
data_hierfstat <- genind2hierfstat(salmon)  

#Create a vector that assigns each individual (total 823) to one panmictic population
one_pop <- rep(1,823)
#Creates a vector of the actual population (or subpopulation) assignments 
levels <- data_hierfstat[,1]  
#Create a data frame only of locus data (formatted for hierfstat)
locus_data <- data_hierfstat[,2:12]  
#Permutes the populations and estimates G statistics 
testwithin <- test.within(locus_data,within=one_pop,test.lev=levels,nperm=1000) 
testwithin$p.val    

#################################################################################################################################
### QUESTION: 
###  Is there evidence that the samples comprise one single population?
##################################################################################################################################


######################################################################################################
# CODE FOR PART C - DERIVE GENETIC DISTANCES BETWEEN POPULATIONS (FST) 
######################################################################################################

#### Calculate pairwise Fst per population
##Function in package "diveRsity" that calculates pairwise W&C's Fst; 
#Samples individuals with replacement to create a new dataset and recalculates Fst 
#repeats 100 times to generate a mean and 95% CI's 

pop_stats <- fstOnly(infile="class_data_genepop_new.gen",outfile="Salmon_Pairwise_Fst",
                     gp=3,bs_pairwise=TRUE,bootstraps=100,parallel=TRUE) 

### Open the output from the function above 
#and examine the Pairwise Fst values with confidence intervals

#################################################################################################################################
### QUESTION: 
### EXAMINE THE CONFIDENCE INTERVALS FOR ALL PAIRWISE COMPARISIONS. 
# WHICH POPULATION PAIRWISE COMPARISONS ARE SIGNIFICANTLY DIFFERENT FROM ZERO? 
#HINT: THEY ARE DIFFERENT IF THEIR CI'S DO NOT INCLUDE ZERO
write.table(pop_stats,file="FSTresults.csv",sep=",",row.names=T)
##################################################################################################################################

#### Before lab, I took the list of actual pairwise Fst values 
#and transformed them into a simple pairwise matrix (file=Fst.csv). 
#Import this matrix now. 
pairwise_fst <- read.csv("Fst.csv",header=TRUE,row.names=1)
pairwise_fst
fst <- as.dist(pairwise_fst) ## Convert pairwise_fst to a distance object
fst  ### we'll use this matrix later in the analysis


######################################################################################################
# CODE FOR PART D - VISUALIZE THE DATA BY CONSTRUCTING A NEIGHBOR-JOINING TREE 
######################################################################################################

### Now let's make a phylogenetic tree!!
### We have to import actual genotypic data to construct 
#population distances and a phylogenetic tree
genotype_data <- read.csv("Genotype_data.csv",header=TRUE) 

### Remove the first column of the data file, 
#which contains the sample ID, so that the matrix only has genotypes
genotypes <- genotype_data[,2:12] 
## Creates a vector containing the population assignments of each individual
pop_labels <- c(rep("JA",90),rep("KO",86),rep("AD",45),
                rep("AI",92),rep("AT",45),rep("UP",87),rep("KI",106),rep("HS",89),
                rep("WA",69),rep("SG",94),rep("PS",20))  
## Converts labels to "factor"
pops <- as.factor(pop_labels) 

###Convert to format compatible for use with phylogenetic tree 
#packages in R; stores data as gene frequencies
genet_file <- char2genet(genotypes,pops,complete=TRUE)  

#### Computes Nei's distances between populations
Nei_dist <- dist.genet(genet_file,method=1)    
## View the distance matrix
Nei_dist 
###Construct a Neighbor-joining tree from the distance matrix
tree <- nj(Nei_dist) 
### Plot the Neighbor-joining tree
plot.phylo(tree)  


#### To bootstrap the tree
## Write a function to bootstap the NJ tree
func <- function(x) nj(dist.genet(char2genet(x,pops)))  
bootstraps <- boot.phylo(tree,genotypes,func,B=100,block=1,rooted=TRUE)

#### Plot your original tree
plot.phylo(tree)    
###Add the bootstrap values from 
#the 100 trees made from randomly sampling your data
nodelabels(bootstraps)  
bootstraps


#################################################################################################################################
### SAVE THE IMAGE OF YOUR TREE FOR YOUR LAB WRITE-UP!!
##################################################################################################################################


###################################################################################################################################
# CODE FOR PART E - VISUALIZE THE DATA THROUGH MULTIDIMENSIONAL SCALING ANALYSIS (ALSO KNOWN AS A PRINCIPAL COORDINATES ANALYSIS)
##################################################################################################################################

######## Conduct a Principal Coordinates Analysis on the Nei_dist matrix generated previously

#### Conducts the PCoA; choose to retain 2 axes
pcoa<-dudi.pco(Nei_dist,scannf=FALSE,nf=2)   
### Lists the eigenvalues for the axes; 
#represents the variation explained by each axis
pcoa$eig  
### Plots the populations on the first two axes
scatter(pcoa,xax=1,yax=2,clab.row=1,posieig="topright") 

### This code yields the proportion of 
#total variation explained by the x-axis in your plot
variance_explained1<-pcoa$eig[1]/sum(pcoa$eig)   
variance_explained1 
## This code yields the proportion of 
#total variation explained by the y-axis in your plot
variance_explained2<-pcoa$eig[2]/sum(pcoa$eig)   
variance_explained2

#################################################################################################################################
### SAVE THE IMAGE OF YOUR TREE FOR YOUR LAB WRITE-UP!!

### QUESTION: HOW MUCH VARIATION DOES THE X-AXIs EXPLAIN? THE Y-AXIS? 
#What population relationships are described by these axes?
##################################################################################################################################


###################################################################################################################################
# CODE FOR PART F - VISUALIZE THE DATA: 
#THE RELATIONSHIP BETWEEN GENETIC AND GEOGRAPHIC DISTANCE
##################################################################################################################################

###Read in the matix of geographic distances
geo_dists <- read.csv("geographic_distances.csv",header=TRUE,row.names=1)  
##convert to a distance object
geo_dist_object <- as.dist(geo_dists)
geo_dist_object
par(mar=c(5,5,5,5))
plot(geo_dist_object,fst,type='p',pch=19,col="blue",
     main="Genetic vs. Geographic Distance", xlab="Geographic Distance (km)",
     ylab="Genetic Distance (W&C Fst)",cex.main=2.5,cex.lab=2,cex.axis=2)

### Fits a linear model to the data
relationship <- lm(fst~geo_dist_object)  
### Get the results of the linear regression model
relationship   
### Plot the regression line to the graph of genetic vs. geographic distance
abline(relationship,lwd=3,col="black")  
### Gives the statistical results: the r-squared value and 
#the signficance of the explanatory variable
summary(relationship)    

#################################################################################################################################
### QUESTION: 
### WHAT IS THE R-SQUARED VALUE OF THE LINEAR MODEL?
### IS THERE A SIGNIFICANT RELATIONSHIP BETWEEN GENETIC DISTANCE 
#AND GEOGRAPHIC DISTANCE? 
##################################################################################################################################
