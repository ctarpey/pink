### Performing LFMM analyis on population data
###    using the package LEA from http://www.bioconductor.org/packages/release/bioc/html/LEA.html
### Carolyn Tarpey | September 2015 
### ---------------------------------------

#set the working directory 
#setwd("C:/Users/Carolyn/Documents/GitHub/pink/LFMM")
#setwd("G:/Analysis/Pop_analysis/Populations_b3_may/LFMM")

setwd("U:/LFMM")
setwd("C:/Users/ctarpey/Desktop/copy_of_Udrive/LFMM_complete")

#install.packages("biocLite", source = "https://bioconductor.org/biocLite.R")

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("LEA")
library(LEA)


#to load project, use:
#project <- load.snmfProject("G:/Analysis/Pop_analysis/Populations_b3_may/LFMM/LFMM.snmfProject")

#names of the populations in order
names <- read.table("POPnames.txt", row.names = 1 )

#convert the edited ped file into the LFMM format
ped2lfmm("LFMM_out_ped.ped", "LFMM_in.lfmm" )

#convert the lfmm format into the genotype format
genotype <- lfmm2geno("LFMM_in.lfmm", "LFMM_in.geno")

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#genome scan for selection using env variables
#running it with the standard longitude as the first variable and the stand lat as the second 

#project <- "SecondTest"

obj.lfmm_4_LL_d1 <- lfmm("LFMM_in.lfmm","Long_Lat_new.env", project = "continue", K = 4,  
      all = FALSE, missing.data = TRUE, CPU = 1, d = 1, 
     iterations = 10000, burnin = 5000, repetitions = 1)

obj.lfmm_5_LL_d1 <-lfmm("LFMM_in.lfmm","Long_Lat_new.env", project = "continue", K = 5, 
     all = FALSE, missing.data = TRUE, CPU = 1, d = 1,
     iterations = 10000, burnin = 5000, repetitions = 1 )

obj.lfmm_6_LL_d1 <-lfmm("LFMM_in.lfmm","Long_Lat_new.env", project = "continue",  K = 6, 
      all = FALSE, missing.data = TRUE, CPU = 1, d = 1,
     iterations = 10000, burnin = 5000, repetitions = 1 )

obj.lfmm_4_LL_d2 <- lfmm("LFMM_in.lfmm","Long_Lat_new.env", project = "continue",  K = 4, 
    all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
    iterations = 10000, burnin = 5000, repetitions = 1)

obj.lfmm_5_LL_d2 <-lfmm("LFMM_in.lfmm","Long_Lat_new.env", project = "continue",  K = 5, 
      all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
     iterations = 10000, burnin = 5000, repetitions = 1)

obj.lfmm_6_LL_d2 <-lfmm("LFMM_in.lfmm","Long_Lat_new.env", project = "continue",  K = 6, 
      all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
     iterations = 10000, burnin = 5000, repetitions = 1 )

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#genome scan for selection using env variables
#running it with the temperature and precipitation variables

########Repetition of 1



#all pc axis together- one and two. 

obj.lfmm_4_PT_all <- lfmm("LFMM_in.lfmm","precip_temp.env", project = "continue", K = 4,  
                         all = TRUE, missing.data = TRUE, CPU = 1, 
                         iterations = 10000, burnin = 5000, repetitions = 1)

obj.lfmm_5_PT_all <-lfmm("LFMM_in.lfmm","precip_temp.env", project = "continue", K = 5, 
                        all = TRUE, missing.data = TRUE, CPU = 1,
                        iterations = 10000, burnin = 5000, repetitions = 1 )

obj.lfmm_6_PT_all <-lfmm("LFMM_in.lfmm","precip_temp.env", project = "continue",  K = 6, 
                        all = TRUE, missing.data = TRUE, CPU = 1,
                        iterations = 10000, burnin = 5000, repetitions = 1 )

## pc axis one alone d= 1

obj.lfmm_4_PT_d1 <-lfmm("LFMM_in.lfmm","precip_temp.env", project = "continue",  K = 4, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 1,
                        iterations = 10000, burnin = 5000, repetitions = 1 )

obj.lfmm_5_PT_d1 <-lfmm("LFMM_in.lfmm","precip_temp.env", project = "continue", K = 5, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 1,
                        iterations = 10000, burnin = 5000, repetitions = 1 )

obj.lfmm_6_PT_d1 <- lfmm("LFMM_in.lfmm","precip_temp.env", project = "continue",  K = 6, 
                         all = FALSE, missing.data = TRUE, CPU = 1, d = 1, 
                         iterations = 10000, burnin = 5000, repetitions = 1)

## pc axis two alone d = 2

obj.lfmm_4_PT_d2 <-lfmm("LFMM_in.lfmm","precip_temp.env", project = "continue",  K = 4, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                        iterations = 10000, burnin = 5000, repetitions = 1 )

obj.lfmm_5_PT_d2 <-lfmm("LFMM_in.lfmm","precip_temp.env", project = "continue",  K = 5, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                        iterations = 10000, burnin = 5000, repetitions = 1 )

obj.lfmm_6_PT_d2 <-lfmm("LFMM_in.lfmm","precip_temp.env", project = "continue",  K = 6, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                        iterations = 10000, burnin = 5000, repetitions = 1 )


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
#genome scan for selection using env variables
#running it with the standard longitude as the first variable and the stand lat as the second 
#repetition of 5


obj.lfmm_4_LL_d1_5 <- lfmm("LFMM_in.lfmm","Long_Lat_5_new.env", project = "continue", K = 4,  
                         all = FALSE, missing.data = TRUE, CPU = 1, d = 1, 
                         iterations = 10000, burnin = 5000, repetitions = 5)

obj.lfmm_5_LL_d1_5 <-lfmm("LFMM_in.lfmm","Long_Lat_5_new.env", project = "continue", K = 5, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 1,
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_6_LL_d1_5 <-lfmm("LFMM_in.lfmm","Long_Lat_5_new.env", project = "continue",  K = 6, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 1,
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_4_LL_d2_5 <- lfmm("LFMM_in.lfmm","Long_Lat_5_new.env", project = "continue",  K = 4, 
                         all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                         iterations = 10000, burnin = 5000, repetitions = 5)

obj.lfmm_5_LL_d2_5 <-lfmm("LFMM_in.lfmm","Long_Lat_5_new.env", project = "continue",  K = 5, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                        iterations = 10000, burnin = 5000, repetitions = 5)

obj.lfmm_6_LL_d2_5 <-lfmm("LFMM_in.lfmm","Long_Lat_5_new.env", project = "continue",  K = 6, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                        iterations = 10000, burnin = 5000, repetitions = 5 )  
  
  
  

########Repetition of 5


#all pc axis together- one and two. 

obj.lfmm_4_PT_all_5 <- lfmm("LFMM_in.lfmm","precip_temp_all.env", project = "continue", K = 4,  
                          all = TRUE, missing.data = TRUE, CPU = 1, 
                          iterations = 10000, burnin = 5000, repetitions = 5)

obj.lfmm_5_PT_all_5 <-lfmm("LFMM_in.lfmm","precip_temp_all.env", project = "continue", K = 5, 
                         all = TRUE, missing.data = TRUE, CPU = 1,
                         iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_6_PT_all_5 <-lfmm("LFMM_in.lfmm","precip_temp_all.env", project = "continue",  K = 6, 
                         all = TRUE, missing.data = TRUE, CPU = 1,
                         iterations = 10000, burnin = 5000, repetitions = 5 )

## pc axis one alone d= 1

obj.lfmm_4_PT_d1_5 <-lfmm("LFMM_in.lfmm","precip_temp_all.env", project = "continue",  K = 4, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 1,
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_5_PT_d1_5 <-lfmm("LFMM_in.lfmm","precip_temp_all.env", project = "continue", K = 5, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 1,
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_6_PT_d1_5 <- lfmm("LFMM_in.lfmm","precip_temp_all.env", project = "continue",  K = 6, 
                         all = FALSE, missing.data = TRUE, CPU = 1, d = 1, 
                         iterations = 10000, burnin = 5000, repetitions = 5)

## pc axis two alone d = 2

obj.lfmm_4_PT_d2_5 <-lfmm("LFMM_in.lfmm","precip_temp_all.env", project = "continue",  K = 4, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_5_PT_d2_5 <-lfmm("LFMM_in.lfmm","precip_temp_all.env", project = "continue",  K = 5, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                        iterations = 10000, burnin = 5000, repetitions = 5 )

obj.lfmm_6_PT_d2_5 <-lfmm("LFMM_in.lfmm","precip_temp_all.env", project = "continue",  K = 6, 
                        all = FALSE, missing.data = TRUE, CPU = 1, d = 2, 
                        iterations = 10000, burnin = 5000, repetitions = 5 )

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  #from here use the Post_process_LFMM.r code to process the output of these runs. 
