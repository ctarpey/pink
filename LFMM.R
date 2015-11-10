### Performing LFMM analyis on population data
###    using the package LEA from http://www.bioconductor.org/packages/release/bioc/html/LEA.html
### Carolyn Tarpey | September 2015 
### ---------------------------------------

#set the working directory 
#setwd("C:/Users/Carolyn/Documents/GitHub/pink/LFMM")

setwd("G:/Analysis/Pop_analysis/Populations_b3_may/LFMM")

#install.packages("biocLite", source = "https://bioconductor.org/biocLite.R")
#biocLite("LEA")

library(LEA)


#to load project, use:
project <- load.snmfProject("G:/Analysis/Pop_analysis/Populations_b3_may/LFMM/LFMM.snmfProject")

#names of the populations in order
names <- read.table("G:/Analysis/Pop_analysis/Populations_b3_may/LFMM/POPnames.txt", row.names = 1 )

#convert the edited ped file into the LFMM format
ped2lfmm("G:/Analysis/Pop_analysis/Populations_b3_may/LFMM/ConvertFiles/LFMM_out_ped.ped",
         "G:/Analysis/Pop_analysis/Populations_b3_may/LFMM/LFMM2.lfmm" )

#convert the lfmm format into the genotype format
genotype <- lfmm2geno("G:/Analysis/Pop_analysis/Populations_b3_may/LFMM/LFMM2.lfmm", 
                      "G:/Analysis/Pop_analysis/Populations_b3_may/LFMM/LFMM.geno")

#figuring out K, the number of ancestral populations
obj.snmf = snmf(genotype, K = 1:14, entropy = TRUE, ploidy = 2, rep = 10, project = "new")

#show and summarize the results 
show(obj.snmf)
summary(obj.snmf)

#Plotting the values of the cross-entropy criterion for each K; the value for which the function plateaus or increases is our estimate of K
plot(obj.snmf)

#plot cross-entropy criterion of all runs of the project
plot(obj.snmf, lwd= 5, col = "red", pch=1)

#to get the cross-entropy of each run for K = 4:8
ce_4 <- cross.entropy(obj.snmf, K= 4)
ce_5 <- cross.entropy(obj.snmf, K= 5)
ce_6 <- cross.entropy(obj.snmf, K= 6)
ce_7 <- cross.entropy(obj.snmf, K= 7)
ce_8 <- cross.entropy(obj.snmf, K= 8)

#select the run with the lowest cross-entropy 
best_4 <- which.min(ce_4)
best_5 <- which.min(ce_5)
best_6 <- which.min(ce_6)
best_7 <- which.min(ce_7)
best_8 <- which.min(ce_8)

#visualize a bar plot of ancestry coefficients 
windows(height=7, width=7)

barplot(t(Q(obj.snmf, K = 4, run = best_4)), col = 1:4)
barplot(t(Q(obj.snmf, K = 5, run = best_5)), col = 1:5)
barplot(t(Q(obj.snmf, K = 6, run = best_6)), col = 1:6)
barplot(t(Q(obj.snmf, K = 7, run = best_7)), col = 1:7)
barplot(t(Q(obj.snmf, K = 8, run = best_8)), col = 1:8)

#genome scan for selection using env variables
#running it with the standard longitude as the first variable and the stand lat as the second 


obj.lfmm_4_d1 <- lfmm("LFMM2.lfmm","Long_Lat.env", project = "new", 
     K = 4,  all = FALSE, missing.data = TRUE, CPU = 1, d = 1, 
     iterations = 10000, burnin = 5000, repetitions = 1)

obj.lfmm_5_d1 <-lfmm("LFMM2.lfmm","Long_Lat.env", project = "new", K = 5, 
     all = FALSE, missing.data = TRUE, CPU = 1, d = 1,
     iterations = 10000, burnin = 5000, repetitions = 1 )

obj.lfmm_6_d1 <-lfmm("LFMM2.lfmm","Long_Lat.env", project = "new",  K = 6, 
      all = FALSE, missing.data = TRUE, CPU = 1, d = 1,
     iterations = 10000, burnin = 5000, repetitions = 1 )

obj.lfmm_7_d1 <-lfmm("LFMM2.lfmm","Long_Lat.env", project = "new",  K = 7, 
     all = FALSE, missing.data = TRUE, CPU = 1, d = 1,
     iterations = 10000, burnin = 5000, repetitions = 1 )

obj.lfmm_8_d1 <-lfmm("LFMM2.lfmm","Long_Lat.env", project = "new", K = 8, 
      all = FALSE, missing.data = TRUE, CPU = 1, d = 1,
     iterations = 10000, burnin = 5000, repetitions = 1 )

obj.lfmm_4_d2 <- lfmm("LFMM2.lfmm","Long_Lat.env", project = "new",  K = 4, 
    all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
    iterations = 10000, burnin = 5000, repetitions = 1)

obj.lfmm_5_d2 <-lfmm("LFMM2.lfmm","Long_Lat.env", project = "new",  K = 5, 
      all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
     iterations = 10000, burnin = 5000, repetitions = 1)

obj.lfmm_6_d2 <-lfmm("LFMM2.lfmm","Long_Lat.env", project = "new",  K = 6, 
      all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
     iterations = 10000, burnin = 5000, repetitions = 1 )

obj.lfmm_7_d2 <-lfmm("LFMM2.lfmm","Long_Lat.env", project = "new",  K = 7, 
      all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
     iterations = 10000, burnin = 5000, repetitions = 1 )

obj.lfmm_8_d2 <-lfmm("LFMM2.lfmm","Long_Lat.env", project = "new",  K = 8, 
     all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
     iterations = 10000, burnin = 5000, repetitions = 1 )


#genome scan for selection using env variables
#running it with the standard longitude as the first variable and the stand lat as the second 
#repition of 5

obj.lfmm_4_d1_r5 <- lfmm("LFMM2.lfmm","Long_Lat.env", project = "new", 
                         K = 4,  all = FALSE, missing.data = TRUE, CPU = 1, d = 1, 
                         iterations = 10000, burnin = 5000, repetitions = 5)

obj.lfmm_5_d1_r5 <-lfmm("LFMM2.lfmm","Long_Lat.env", project = "new", K = 5, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 1,
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_6_d1_r5 <-lfmm("LFMM2.lfmm","Long_Lat.env", project = "new",  K = 6, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 1,
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_7_d1_r5 <-lfmm("LFMM2.lfmm","Long_Lat.env", project = "new",  K = 7, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 1,
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_8_d1_r5 <-lfmm("LFMM2.lfmm","Long_Lat.env", project = "new", K = 8, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 1,
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_4_d2_r5 <- lfmm("LFMM2.lfmm","Long_Lat.env", project = "new",  K = 4, 
                         all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                         iterations = 10000, burnin = 5000, repetitions = 5)

obj.lfmm_5_d2_r5 <-lfmm("LFMM2.lfmm","Long_Lat.env", project = "new",  K = 5, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_6_d2_r5 <-lfmm("LFMM2.lfmm","Long_Lat.env", project = "new",  K = 6, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_7_d2_r5 <-lfmm("LFMM2.lfmm","Long_Lat.env", project = "new",  K = 7, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_8_d2_r5 <-lfmm("LFMM2.lfmm","Long_Lat.env", project = "new",  K = 8, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                        iterations = 10000, burnin = 5000, repetitions = 5 )

#genome scan for selection using env variables 
#running it with the bioclims 

obj.lfmm_4_all_temp <- lfmm("LFMM2.lfmm","temp.env", project = "new", 
                         K = 4,  all = TRUE, missing.data = TRUE, CPU = 2,  
                         iterations = 10000, burnin = 5000, repetitions = 5)

obj.lfmm_5_all_temp <-lfmm("LFMM2.lfmm","temp.env", project = "new", K = 5, 
                        all = TRUE, missing.data = TRUE, CPU = 2, 
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_6_all_temp <-lfmm("LFMM2.lfmm","temp.env", project = "new",  K = 6, 
                        all = TRUE, missing.data = TRUE, CPU = 2,
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_7_all_temp <-lfmm("LFMM2.lfmm","temp.env", project = "new",  K = 7, 
                        all = TRUE, missing.data = TRUE, CPU = 2,
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_8_all_temp<-lfmm("LFMM2.lfmm","temp.env", project = "new", K = 8, 
                        all = TRUE, missing.data = TRUE, CPU = 2, 
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_4_d1_temp <- lfmm("LFMM2.lfmm","temp.env", project = "new",  K = 4, 
                         all = FALSE, missing.data = TRUE, CPU = 2, d = 1, 
                         iterations = 10000, burnin = 5000, repetitions = 5)

obj.lfmm_5_d1_temp <-lfmm("LFMM2.lfmm","temp.env", project = "new",  K = 5, 
                        all = FALSE, missing.data = TRUE, CPU = 2, d = 1, 
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_6_d1_temp <-lfmm("LFMM2.lfmm","temp.env", project = "new",  K = 6, 
                        all = FALSE, missing.data = TRUE, CPU = 2, d = 1, 
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_7_d1_temp <-lfmm("LFMM2.lfmm","temp.env", project = "new",  K = 7, 
                        all = FALSE, missing.data = TRUE, CPU = 2, d = 1, 
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_8_d1_temp <-lfmm("LFMM2.lfmm","temp.env", project = "new",  K = 8, 
                        all = FALSE, missing.data = TRUE, CPU = 2, d = 1, 
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_4_d2_temp <- lfmm("LFMM2.lfmm","temp.env", project = "new",  K = 4, 
                           all = FALSE, missing.data = TRUE, CPU = 2, d = 2, 
                           iterations = 10000, burnin = 5000, repetitions = 5)

obj.lfmm_5_d2_temp <-lfmm("LFMM2.lfmm","temp.env", project = "new",  K = 5, 
                          all = FALSE, missing.data = TRUE, CPU = 2, d = 2, 
                          iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_6_d2_temp <-lfmm("LFMM2.lfmm","temp.env", project = "new",  K = 6, 
                          all = FALSE, missing.data = TRUE, CPU = 2, d = 2, 
                          iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_7_d2_temp <-lfmm("LFMM2.lfmm","temp.env", project = "new",  K = 7, 
                          all = FALSE, missing.data = TRUE, CPU = 2, d = 2, 
                          iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_8_d2_temp <-lfmm("LFMM2.lfmm","temp.env", project = "new",  K = 8, 
                          all = FALSE, missing.data = TRUE, CPU = 2, d = 2, 
                          iterations = 10000, burnin = 5000, repetitions = 5 )
