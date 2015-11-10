### Performing LFMM analyis on population data
###    using the package LEA from http://www.bioconductor.org/packages/release/bioc/html/LEA.html
### Carolyn Tarpey | October 2015 
### ---------------------------------------


#source("https://bioconductor.org/biocLite.R")
#biocLite()
#biocLite("LEA")

library(LEA)


#The lfmm function returns a project object containing all lfmm runs. 
#When performing additional runs, the function enables the project to be included as a parameter 
#to add more runs. Performing several runs for various values of the number of latent factors (K) 
#is recommended. 

#these are the steps listed in the manual for post-processing the LFMM runs: 


#setwd("U:/LFMM")
#setwd("G:/Analysis/Pop_analysis/Populations_b3_may/LFMM/LFMM_complete")
setwd("C:/Users/ctarpey/Desktop/LFMM_sm_runs/NAO")



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ LATITUDE and LONGITUDE


#to load project, use:
project_LL <- load.lfmmProject("NAO_lfmm_NAO_long_lat_PT.lfmmProject")

# show the project
show(project_LL)

# summary of the project

summary(project_LL)


# get the zscores of each run for K 
zs_3_LL = z.scores(project_LL, d = 1, K = 3)
zs_2_LL = z.scores(project_LL, d = 1, K = 2)

zs_3_LL_2 = z.scores(project_LL, d = 2, K = 3)
zs_2_LL_2 = z.scores(project_LL, d = 2, K = 2)


# Combine the z-scores using the Stouffer method
zs_3.stouffer_LL = apply(zs_3_LL, MARGIN = 1, median)
zs_2.stouffer_LL = apply(zs_2_LL, MARGIN = 1, median)

zs_3.stouffer_LL_2 = apply(zs_3_LL, MARGIN = 1, median)
zs_2.stouffer_LL_2 = apply(zs_2_LL, MARGIN = 1, median)

# calculate the inflation factor
lambda_3_LL = median(zs_3.stouffer_LL^2)/.456
lambda_2_LL = median(zs_2.stouffer_LL^2)/.456

lambda_3_LL_2 = median(zs_3.stouffer_LL_2^2)/.456
lambda_2_LL_2 = median(zs_2.stouffer_LL_2^2)/.456

# calculate adjusted p-values
cp_3.values_LL = pchisq(zs_3.stouffer_LL^2/lambda_3_LL, df = 1, lower = FALSE)
cp_2.values_LL = pchisq(zs_2.stouffer_LL^2/lambda_2_LL, df = 1, lower = FALSE)


cp_3.values_LL_2 = pchisq(zs_3.stouffer_LL_2^2/lambda_3_LL_2, df = 1, lower = FALSE)
cp_2.values_LL_2 = pchisq(zs_2.stouffer_LL_2^2/lambda_2_LL_2, df = 1, lower = FALSE)



for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 3 expected FDR:", alpha))
  L = length(cp_3.values_LL)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_3.values_LL) < alpha * (1:L) / L)
  candidates = order(cp_3.values_LL)[w]
  
  # estimated FDR and True Positif
  estimated.FDR = length(which(candidates <= 350))/length(candidates)
  estimated.TP = length(which(candidates > 350))/50
  print(paste("K = 3 FDR:", estimated.FDR, "True Positive:", estimated.TP))
}


for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 2 expected FDR:", alpha))
  L = length(cp_2.values_LL)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_2.values_LL) < alpha * (1:L) / L)
  candidates = order(cp_2.values_LL)[w]
  
  # estimated FDR and True Positif
  estimated.FDR = length(which(candidates <= 350))/length(candidates)
  estimated.TP = length(which(candidates > 350))/50
  print(paste("K = 2 FDR:", estimated.FDR, "True Positive:", estimated.TP))
}


for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 3 expected FDR:", alpha))
  L = length(cp_3.values_LL_2)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_3.values_LL_2) < alpha * (1:L) / L)
  candidates_LL = order(cp_3.values_LL_2)[w]
  
  # estimated FDR and True Positif
  estimated.FDR_LL = length(which(candidates_LL <= 350))/length(candidates_LL)
  estimated.TP_LL = length(which(candidates_LL > 350))/50
  print(paste("k = 3 FDR:", estimated.FDR_LL, "True Positive:", estimated.TP_LL))
}



for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 2 expected FDR:", alpha))
  L = length(cp_2.values_LL_2)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_2.values_LL_2) < alpha * (1:L) / L)
  candidates_LL = order(cp_2.values_LL_2)[w]
  
  # estimated FDR and True Positif
  estimated.FDR_LL = length(which(candidates_LL <= 350))/length(candidates_LL)
  estimated.TP_LL = length(which(candidates_LL > 350))/50
  print(paste("k = 2 FDR:", estimated.FDR_LL, "True Positive:", estimated.TP_LL))
}




#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ PRECIPITATION and TEMPERATURE

#to load project, use:
project_PT <- load.lfmmProject("NAO_lfmm_NAO_PT.lfmmProject")

# show the project
show(project_PT)

# summary of the project
summary(project_PT)

# get the zscores of each run for K 
zs_3_PT = z.scores(project_PT, d = 1, K = 3)
zs_2_PT = z.scores(project_PT, d = 1, K = 2)

zs_3_PT_2 = z.scores(project_PT, d = 2, K = 3)
zs_2_PT_2 = z.scores(project_PT, d = 2, K = 2)


# Combine the z-scores using the Stouffer method
zs_3.stouffer_PT = apply(zs_3_PT, MARGIN = 1, median)
zs_2.stouffer_PT = apply(zs_2_PT, MARGIN = 1, median)

zs_3.stouffer_PT_2 = apply(zs_3_PT, MARGIN = 1, median)
zs_2.stouffer_PT_2 = apply(zs_2_PT, MARGIN = 1, median)

# calculate the inflation factor
lambda_3_PT = median(zs_3.stouffer_PT^2)/.456
lambda_2_PT = median(zs_2.stouffer_PT^2)/.456

lambda_3_PT_2 = median(zs_3.stouffer_PT_2^2)/.456
lambda_2_PT_2 = median(zs_2.stouffer_PT_2^2)/.456

# calculate adjusted p-values
cp_3.values_PT = pchisq(zs_3.stouffer_PT^2/lambda_3_PT, df = 1, lower = FALSE)
cp_2.values_PT = pchisq(zs_2.stouffer_PT^2/lambda_2_PT, df = 1, lower = FALSE)

cp_3.values_PT_2 = pchisq(zs_3.stouffer_PT_2^2/lambda_3_PT_2, df = 1, lower = FALSE)
cp_2.values_PT_2 = pchisq(zs_2.stouffer_PT_2^2/lambda_2_PT_2, df = 1, lower = FALSE)


for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 3 expected FDR:", alpha))
  L = length(cp_3.values_PT)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_3.values_PT) < alpha * (1:L) / L)
  candidates = order(cp_3.values_PT)[w]
  
  # estimated FDR and True Positif
  estimated.FDR = length(which(candidates <= 350))/length(candidates)
  estimated.TP = length(which(candidates > 350))/50
  print(paste("k = 3 FDR:", estimated.FDR, "True Positive:", estimated.TP))
}


for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 2 expected FDR:", alpha))
  L = length(cp_2.values_PT)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_2.values_PT) < alpha * (1:L) / L)
  candidates = order(cp_2.values_PT)[w]
  
  # estimated FDR and True Positif
  estimated.FDR = length(which(candidates <= 350))/length(candidates)
  estimated.TP = length(which(candidates > 350))/50
  print(paste("K = 2 FDR:", estimated.FDR, "True Positive:", estimated.TP))
}

for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 3 expected FDR:", alpha))
  L = length(cp_3.values_PT_2)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_3.values_PT_2) < alpha * (1:L) / L)
  candidates = order(cp_3.values_PT_2)[w]
  
  # estimated FDR and True Positif
  estimated.FDR = length(which(candidates <= 350))/length(candidates)
  estimated.TP = length(which(candidates > 350))/50
  print(paste("k = 3 FDR:", estimated.FDR, "True Positive:", estimated.TP))
}


for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 2 expected FDR:", alpha))
  L = length(cp_2.values_PT_2)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_2.values_PT_2) < alpha * (1:L) / L)
  candidates = order(cp_2.values_PT_2)[w]
  
  # estimated FDR and True Positif
  estimated.FDR = length(which(candidates <= 350))/length(candidates)
  estimated.TP = length(which(candidates > 350))/50
  print(paste("K = 2 FDR:", estimated.FDR, "True Positive:", estimated.TP))
}

