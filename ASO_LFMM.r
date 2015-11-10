### Performing LFMM analyis on population data
###    using the package LEA from http://www.bioconductor.org/packages/release/bioc/html/LEA.html
### Carolyn Tarpey | October 2015 
### ---------------------------------------

#set the working directory 
setwd("C:/Users/ctarpey/Desktop/LFMM_sm_runs/ASO")

library(LEA)

#names of the populations in order
names <- read.table("POPnames.txt", row.names = 1 )

#convert the edited ped file into the LFMM format
ped2lfmm("ASO_ped.ped","ASO_lfmm.lfmm" )

#convert the lfmm format into the genotype format
ASOgenotype <- lfmm2geno("ASO_lfmm.lfmm", "ASO_geno.geno")

#figuring out K, the number of ancestral populations
ASO.snmf = snmf(ASOgenotype, K = 1:4, entropy = TRUE, ploidy = 2, rep = 10, project = "new")


#$$$$$$$$$$$$$$$$$


#show and summarize the results 
show(ASO.snmf)
summary(ASO.snmf)

#Plotting the values of the cross-entropy criterion for each K; the value for which the function plateaus or increases is our estimate of K
plot(ASO.snmf)

#plot cross-entropy criterion of all runs of the project
plot(ASO.snmf, lwd= 5, col = "red", pch=1)

#to get the cross-entropy of each run for K = 4:8
ce_1 <- cross.entropy(ASO.snmf, K= 1)
ce_2 <- cross.entropy(ASO.snmf, K= 2)
ce_3 <- cross.entropy(ASO.snmf, K= 3)
ce_4 <- cross.entropy(ASO.snmf, K= 4)


#select the run with the lowest cross-entropy 
best_1 <- which.min(ce_1)
best_2 <- which.min(ce_2)
best_3 <- which.min(ce_3)
best_4 <- which.min(ce_4)


#visualize a bar plot of ancestry coefficients 
windows(height=7, width=7)

barplot(t(Q(ASO.snmf, K = 1, run = best_1)), col = 1)
barplot(t(Q(ASO.snmf, K = 2, run = best_2)), col = 1:2)
barplot(t(Q(ASO.snmf, K = 3, run = best_3)), col = 1:3)
barplot(t(Q(ASO.snmf, K = 4, run = best_4)), col = 1:4)

#genome scan for selection using env variables
#running it with the standard longitude as the first variable and the stand lat as the second 


ASO.lfmm_LL_k1_d1_r5 <- lfmm("ASO_lfmm.lfmm","ASO_lat_long.env", project = "continue", 
                             K = 1,  all = FALSE, missing.data = TRUE, CPU = 1, d = 1, 
                             iterations = 10000, burnin = 5000, repetitions = 5)

ASO.lfmm_LL_k2_d1_r5 <- lfmm("ASO_lfmm.lfmm","ASO_lat_long.env", project = "continue", 
                             K = 2,  all = FALSE, missing.data = TRUE, CPU = 1, d = 1, 
                             iterations = 10000, burnin = 5000, repetitions = 5)

ASO.lfmm_LL_k1_d2_r5 <- lfmm("ASO_lfmm.lfmm","ASO_lat_long.env", project = "continue", 
                             K = 1,  all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                             iterations = 10000, burnin = 5000, repetitions = 5)

ASO.lfmm_LL_k2_d2_r5 <- lfmm("ASO_lfmm.lfmm","ASO_lat_long.env", project = "continue", 
                             K = 2,  all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                             iterations = 10000, burnin = 5000, repetitions = 5)



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#genome scan for selection using env variables
#running it with the temperature and precipitation variables

########Repetition of 1



#all pc axis together- one and two. 

ASO.lfmm_PT_k1_all_r5 <- lfmm("ASO_lfmm.lfmm","ASO_PT.env", project = "continue", K = 1,  
                              all = TRUE, missing.data = TRUE, CPU = 1, 
                              iterations = 10000, burnin = 5000, repetitions = 5)

ASO.lfmm_PT_k2_all_r5 <- lfmm("ASO_lfmm.lfmm","ASO_PT.env", project = "continue", K = 2,  
                              all = TRUE, missing.data = TRUE, CPU = 1, 
                              iterations = 10000, burnin = 5000, repetitions = 5)

## pc axis one alone d= 1
ASO.lfmm_PT_k1_d1_r5 <- lfmm("ASO_lfmm.lfmm","ASO_PT.env", project = "continue", K = 1,  
                             all = FALSE, missing.data = TRUE, CPU = 1, d = 1, 
                             iterations = 10000, burnin = 5000, repetitions = 5)

ASO.lfmm_PT_k2_d1_r5 <- lfmm("ASO_lfmm.lfmm","ASO_PT.env", project = "continue", K = 2,  
                             all = FALSE, missing.data = TRUE, CPU = 1, d = 1, 
                             iterations = 10000, burnin = 5000, repetitions = 5)

## pc axis two alone d = 2

ASO.lfmm_PT_k1_d2_r5 <- lfmm("ASO_lfmm.lfmm","ASO_PT.env", project = "continue", K = 1,  
                             all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                             iterations = 10000, burnin = 5000, repetitions = 5)

ASO.lfmm_PT_k2_d2_r5 <- lfmm("ASO_lfmm.lfmm","ASO_PT.env", project = "continue", K = 2,  
                             all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                             iterations = 10000, burnin = 5000, repetitions = 5)

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  #from here use the Post_process_LFMM.r code to process the output of these runs. 
  