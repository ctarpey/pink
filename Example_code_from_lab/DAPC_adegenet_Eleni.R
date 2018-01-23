### Create DAPC plots for Populations
###   Written by Eleni 
###      
###   Eleni Petrou | January 2018
### ---------------------------------------

#This code requires the lists of paralogs from HDPLOT_forPinks.R
#The one SNP per tag genotypes should be filtered through SecondPinkFiltering.R already
#The haplotype file is from STACKS using the 31485 whitelist after initialPinkFiltering.R
#Load the necessary libaries

library(adegenet)
library(ggplot2)
library(hierfstat)

# Set your working directory
setwd("D:/sequencing_data/Herring_Coastwide_PopulationStructure/Genepop_files")

################################################### DAPC using known groups
# What is a Discriminant Analysis of Principal Components (DAPC)?
# DAPC is a multivariate approach that relies on finding synthetic variables built as linear combinations of alleles. 
# DAPC optimizes the variance between groups, while minimizing the variance within groups: it seeks synthetic variables, 
# the discriminant functions, which show differences between groups as best as possible while minimizing variation within clusters. 


# Let's try this shit out. 


###############################################################################################
# READ IN DATA


# First, read in your data as a genepop file. Specify how many characters code each allele with ncode. 

#all pops
#data_all_loci <-read.genepop("batch_1.gen", ncode = 2)

#BC/AK pops only
  #data_all_loci <-read.genepop("batch_1_BC_AK.gen", ncode = 2)

#BC/AK pops only- no late spawners
  #data_all_loci <-read.genepop("batch_1_BC_AK_nolatespawn.gen", ncode = 2)

#BC- primary spawners only
  data_all_loci <-read.genepop("batch_1_BC_primary.gen", ncode = 2)

#specify population names
# all pops
mypop_names <- c("Bern16", "Bute04", "Case07",  "ChPt14",   "ChPt16",
                 "ElBy15", "Ells15", "Gabr15" ,
                 "Harr14", "Knight",
                 "Kwak15", "Mass03",      
                 "Mass16", "Metla02",
                 "Newb14" , "Port14", "PtGb14",   
                 "QlBy14",  "RivIn01",    
                "Sali16",  "Sitka17", "Skid14",
                 "Skid99", "SmBy15", "Spch14",  
                 "SpCh15", "Squa14", "Venn99")

#BC/AK pops only
#mypop_names <- c("Bern16", "Bute04",
#                  "Ells15", "Gabr15" ,
#                "Harr14", "Knight",
#                 "Kwak15", "Mass03",      
#                  "Mass16", "Metla02",
#                  "Newb14" , "RivIn01",    
#                  "Sali16",  "Sitka17", "Skid14",
#                  "Skid99",  "Spch14",  
#                  "SpCh15",  "Venn99")

#BC/AK pops only- no late spawners
#mypop_names <- c("Bern16", "Bute04",
#                 "Ells15", "Gabr15" ,
#                 "Harr14", "Knight",
#                 "Kwak15", 
#                 "Newb14" , "RivIn01",    
#                 "Sali16",  "Sitka17",  "Spch14",  
#                 "SpCh15",  "Venn99")

#BC-primary only
mypop_names <- c("Ells15", "Gabr15" ,
                "Harr14", 
                 "Kwak15", 
                 "Newb14" ,     
                "Spch14",  
                 "SpCh15",  "Venn99")


# Check out the different levels in the data:
      names(data_all_loci)
  # To access on of these levels, type your genind object name + $ + level name
    # Check locus names
      data_all_loci$loc.fac
    # Check number of alleles per locus and the allele names
      data_all_loci$loc.n.all
      data_all_loci$all.names
    # Check number of populations and their names
      data_all_loci$pop

##################################################################################################
#  CONDUCT DAPC
    
# Find clusters of populations in the data. Start with 100 and look at output.
# Then specify the maximum number of PCAs to retain (by looking at the graph). There is no reason to discard PCAs
# other than computational time, so keep all of them.
# Finally, look at the graph output: Value of BIC versus number of clusters; the optimal number of clusters
# is the minimum BIC score (at the elbow of the graph). 
# DAPC function transforms the data using PCA and then performs a discriminant analysis on the retained principal components. 
#dapc_all <- dapc(data_all_loci,data_all_loci$pop)
myclusters <- find.clusters(data_all_loci)

dapc_all <- dapc(data_all_loci,data_all_loci$pop,n.pca=350,n.da=8) ##Retain all, then identify optimal number by optim.a.score
test_a_score <- optim.a.score(dapc_all)
dapc_all <- dapc(data_all_loci,data_all_loci$pop,n.pca=95,n.da=8) ##25 PC's is the optimal number
test_a_score <- optim.a.score(dapc_all)
###################################################################################################
# PLOTTING THE DAPC

#specify colors for the plot

pop_cols <- rainbow(length(mypop_names))


#2D plot. NOTE THAT LABEL ORDER SHOULD BE IN SAME ORDER AS GENEPOP FILE ORDER
scatter(dapc_all,1,2,scree.da=FALSE,cellipse=0,leg=FALSE,cstar = 0,
        label=mypop_names,
        posi.da="bottomleft",csub=2,col=pop_cols,cex=1.5,
        clabel=1,pch=15,solid=1)


legend(x=4.0,y=1.8,bty='n',
       legend= mypop_names,
       pch=15,col=pop_cols,cex=1.0)


#names(dapc_all)
#write.csv(dapc_all$var.contr,file="DAPC Variable Contributions.csv")

##################################################################################
# SUMMARY OF DAPC ANALYSIS

dapc_all  ####Gives summary of the DAPC - Assign.per.pop gives proportion of successful reassignments to original pops based on DF's
dapc_all$var  ### Proportion of variance conserved by the principal components (18.7%)
dapc_all$prior #### Numeric vector giving prior group probabilities
dapc_all$assign ## Posterior group assignment
dapc_all$posterior
dapc_all$eig[1]/sum(dapc_all$eig)  ### Variance explained by first discriminant function (75.1%)
dapc_all$eig[2]/sum(dapc_all$eig)  ### Variance explained by second discriminant function (17.5%)
dapc_all$eig[3]/sum(dapc_all$eig)

########################################################################################################
##### ASSESSING WHICH LOCI ARE IMPORTANT IN THE DATASET

# We can assess which alleles pull apart the DAPC clusters using the command loadingplot. 
# Variable contributions are stored in the var.contr slot of a dapc object

#Along Axis 1
contrib <- loadingplot(dapc_all$var.contr, axis = 1, lab.jitter = 1)

# Along Axis 2
contrib <- loadingplot(dapc_all$var.contr, axis = 2, lab.jitter = 1)

# Identify structural loci of the DAPC
# The function snpzip identifies the set of alleles which contribute most significantly to phenotypic structure
#test_snpzip<-snpzip(data_all_loci,dapc_all,loading.plot=TRUE,method="median") 
#Keep on pressing "return" to scroll through all the DF axes



########################################################################################################################
# PLOT THE ALLELE FREQUENCIES OF INTERESTING LOCI


# Gather the allele frequencies for loci that are important along DA 1
#DA Axis1
freq8468 <- tab(genind2genpop(data_all_loci[loc = c("L_8468")]), freq = TRUE)
freq30660 <- tab(genind2genpop(data_all_loci[loc = c("L_30660")]), freq = TRUE)
freq22402 <- tab(genind2genpop(data_all_loci[loc = c("L_22402")]), freq = TRUE)
freq37605 <- tab(genind2genpop(data_all_loci[loc = c("L_37605")]), freq = TRUE)
freq7529 <- tab(genind2genpop(data_all_loci[loc = c("L_7529")]), freq = TRUE)
freq24510 <- tab(genind2genpop(data_all_loci[loc = c("L_24510")]), freq = TRUE)
freq23188 <- tab(genind2genpop(data_all_loci[loc = c("L_23188")]), freq = TRUE)
freq6354 <- tab(genind2genpop(data_all_loci[loc = c("L_6354")]), freq = TRUE)

# Make a data frame to store this information

Locus_names <- c( rep("Locus_8468",8) ,rep("Locus_30660",8), rep("Locus_22402",8) , rep("Locus_37605",8) ,
                 rep("Locus_7529",8), rep("Locus_24510", 8), rep("Locus_23188", 8), rep("Locus_6354",8))


# Save your population names as a factor with explicitly specified levels, so that 
# ggplot respects the order of the populations and does not order them alphabetically. 
Pop_name <- factor(c("Port Orchard", "Squaxin", "Similk", "Port Gamble", 
              "Quilcene", "Elliot Bay", "Cherry Point 2014", "Cherry Point 2016"), 
              levels = c("Port Orchard", "Squaxin", "Similk", "Port Gamble", 
                          "Quilcene", "Elliot Bay", "Cherry Point 2014", "Cherry Point 2016"))


allele_freq <- c(freq8468[,1], freq30660[,1], freq22402[,1], freq37605[,1], 
                 freq7529[,1], freq24510[,1], freq23188[,2], freq6354[,1])

allelefreq_df = data.frame(Locus_names, Pop_name, allele_freq)



# Plot these using ggplot
# Set the bw theme (layout) for ggplot.
theme_set(theme_bw()) 

ggplot(data = allelefreq_df, aes(x = Pop_name, y = allele_freq, group = Locus_names)) +
  geom_line(aes(color = Locus_names), size = 1.2) +
  geom_point(aes(color = Locus_names), size = 1.5) +
  labs(x="Population", y="Frequency of major allele")+
  ggtitle("Loci with high loadings on DA 1") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
#######################################################################
# Gather the allele frequencies for loci that are important along DA 1
#DA Axis2
freq12303 <- tab(genind2genpop(data_all_loci[loc = c("L_12303")]), freq = TRUE)
freq20984 <- tab(genind2genpop(data_all_loci[loc = c("L_20984")]), freq = TRUE)
freq22785 <- tab(genind2genpop(data_all_loci[loc = c("L_22785")]), freq = TRUE)
freq18340 <- tab(genind2genpop(data_all_loci[loc = c("L_18340")]), freq = TRUE)
freq30828 <- tab(genind2genpop(data_all_loci[loc = c("L_30828")]), freq = TRUE)
freq32581 <- tab(genind2genpop(data_all_loci[loc = c("L_32581")]), freq = TRUE)
freq26497 <- tab(genind2genpop(data_all_loci[loc = c("L_26497")]), freq = TRUE)
freq14706 <- tab(genind2genpop(data_all_loci[loc = c("L_14706")]), freq = TRUE)

# Make a data frame to store this information

Locus_names2 <- c( rep("Locus_12303",8) ,rep("Locus_20984",8), rep("Locus_22785",8) , rep("Locus_18340",8), 
                   rep("Locus_30828",8), rep("Locus_32581",8), rep("Locus_26497",8), rep("Locus_14706",8))


allele_freq2 <- c(freq12303[,1], freq20984[,2], freq22785[,2], freq18340[,1], freq30828[,1], 
                  freq32581[,1], freq26497[,1], freq14706[,1] )

allelefreq_df2 = data.frame(Locus_names2, Pop_name, allele_freq2)

# Plot these using ggplot
# Set the bw theme (layout) for ggplot.
theme_set(theme_bw()) 

ggplot(data = allelefreq_df2, aes(x = Pop_name, y = allele_freq2, group = Locus_names2)) +
  geom_line(aes(color = Locus_names2), size = 1.2) +
  geom_point(aes(color = Locus_names2), size = 1.5) +
  labs(x="Population", y="Frequency of major allele")+
  ggtitle("Loci with high loadings on DA 2") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))


#########################################################################################################################



#others
freq29131 <- tab(genind2genpop(data_all_loci[loc = c("L_29131")]), freq = TRUE)
matplot(freq29131, pch = c("N", "N"), type = "b", 
        xlab = "population", ylab = "allele frequency", xaxt = "n", cex = 1, main = "Locus 29131")

freq18960 <- tab(genind2genpop(data_all_loci[loc = c("L_18960")]), freq = TRUE)
matplot(freq18960, pch = c("N", "N"), type = "b", 
        xlab = "population", ylab = "allele frequency", xaxt = "n", cex = 1, main = "Locus 18960")




