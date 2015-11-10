### Performing LFMM analyis on population data
###    using the package LEA from http://www.bioconductor.org/packages/release/bioc/html/LEA.html
### Carolyn Tarpey | October 2015 
### ---------------------------------------


source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("LEA")
library(LEA)




#The lfmm function returns a project object containing all lfmm runs. 
#When performing additional runs, the function enables the project to be included as a parameter 
#to add more runs. Performing several runs for various values of the number of latent factors (K) 
#is recommended. 

#these are the steps listed in the manual for post-processing the LFMM runs: 


#setwd("U:/LFMM")
#setwd("G:/Analysis/Pop_analysis/Populations_b3_may/LFMM/LFMM_complete")
setwd("C:/Users/ctarpey/Desktop/copy_of_Udrive/LFMM_complete")



#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ LATITUDE and LONGITUDE


#to load project, use:
project_LL <- load.lfmmProject("LFMM_in_Long_Lat_5_new.lfmmProject")

# show the project
show(project_LL)

# summary of the project
summary(project_LL)


# get the zscores of each run for K 
zs_4_LL = z.scores(project_LL, d = 1, K = 4)
zs_5_LL = z.scores(project_LL, d = 1, K = 5)
zs_6_LL = z.scores(project_LL, d = 1 ,K = 6)

zs_4_LL_2 = z.scores(project_LL, d = 2, K = 4)
zs_5_LL_2 = z.scores(project_LL, d = 2, K = 5)
zs_6_LL_2 = z.scores(project_LL, d = 2 ,K = 6)


# Combine the z-scores using the Stouffer method
zs_4.stouffer_LL = apply(zs_4_LL, MARGIN = 1, median)
zs_5.stouffer_LL = apply(zs_5_LL, MARGIN = 1, median)
zs_6.stouffer_LL = apply(zs_6_LL, MARGIN = 1, median)

zs_4.stouffer_LL_2 = apply(zs_4_LL, MARGIN = 1, median)
zs_5.stouffer_LL_2 = apply(zs_5_LL, MARGIN = 1, median)
zs_6.stouffer_LL_2 = apply(zs_6_LL, MARGIN = 1, median)

# calculate the inflation factor
lambda_4_LL = median(zs_4.stouffer_LL^2)/.456
lambda_5_LL = median(zs_5.stouffer_LL^2)/.456
lambda_6_LL = median(zs_6.stouffer_LL^2)/.456

lambda_4_LL_2 = median(zs_4.stouffer_LL_2^2)/.456
lambda_5_LL_2 = median(zs_5.stouffer_LL_2^2)/.456
lambda_6_LL_2 = median(zs_6.stouffer_LL_2^2)/.456

# calculate adjusted p-values
cp_4.values_LL = pchisq(zs_4.stouffer_LL^2/lambda_4_LL, df = 1, lower = FALSE)
cp_5.values_LL = pchisq(zs_5.stouffer_LL^2/lambda_5_LL, df = 1, lower = FALSE)
cp_6.values_LL = pchisq(zs_6.stouffer_LL^2/lambda_6_LL, df = 1, lower = FALSE)

cp_4.values_LL_2 = pchisq(zs_4.stouffer_LL_2^2/lambda_4_LL_2, df = 1, lower = FALSE)
cp_5.values_LL_2 = pchisq(zs_5.stouffer_LL_2^2/lambda_5_LL_2, df = 1, lower = FALSE)
cp_6.values_LL_2 = pchisq(zs_6.stouffer_LL_2^2/lambda_6_LL_2, df = 1, lower = FALSE)


for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 4 expected FDR:", alpha))
  L = length(cp_4.values_LL)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_4.values_LL) < alpha * (1:L) / L)
  candidates_LL = order(cp_4.values_LL)[w]
  
  # estimated FDR and True Positif
  estimated.FDR_LL = length(which(candidates_LL <= 350))/length(candidates_LL)
  estimated.TP_LL = length(which(candidates_LL > 350))/50
  print(paste("k = 4 FDR:", estimated.FDR_LL, "True Positive:", estimated.TP_LL))
}


for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 5 expected FDR:", alpha))
  L = length(cp_5.values_LL)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_5.values_LL) < alpha * (1:L) / L)
  candidates_LL = order(cp_5.values_LL)[w]
  
  # estimated FDR and True Positif
  estimated.FDR_LL = length(which(candidates_LL <= 350))/length(candidates_LL)
  estimated.TP_LL = length(which(candidates_LL > 350))/50
  print(paste("K = 5 FDR:", estimated.FDR_LL, "True Positive:", estimated.TP_LL))
}

for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 6 expected FDR:", alpha))
  L = length(cp_6.values_LL)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_6.values_LL) < alpha * (1:L) / L)
  candidates = order(cp_6.values_LL)[w]
  
  # estimated FDR and True Positif
  estimated.FDR_LL = length(which(candidates <= 350))/length(candidates)
  estimated.TP_LL = length(which(candidates > 350))/50
  print(paste("K = 6 FDR:", estimated.FDR_LL, "True Positive:", estimated.TP_LL))
}


for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 4 expected FDR:", alpha))
  L = length(cp_4.values_LL_2)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_4.values_LL_2) < alpha * (1:L) / L)
  candidates_LL = order(cp_4.values_LL_2)[w]
  
  # estimated FDR and True Positif
  estimated.FDR_LL = length(which(candidates_LL <= 350))/length(candidates_LL)
  estimated.TP_LL = length(which(candidates_LL > 350))/50
  print(paste("k = 4 FDR:", estimated.FDR_LL, "True Positive:", estimated.TP_LL))
}


for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 5 expected FDR:", alpha))
  L = length(cp_5.values_LL_2)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_5.values_LL_2) < alpha * (1:L) / L)
  candidates_LL = order(cp_5.values_LL_2)[w]
  
  # estimated FDR and True Positif
  estimated.FDR_LL = length(which(candidates_LL <= 350))/length(candidates_LL)
  estimated.TP_LL = length(which(candidates_LL > 350))/50
  print(paste("K = 5 FDR:", estimated.FDR_LL, "True Positive:", estimated.TP_LL))
}

for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 6 expected FDR:", alpha))
  L = length(cp_6.values_LL_2)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_6.values_LL_2) < alpha * (1:L) / L)
  candidates = order(cp_6.values_LL_2)[w]
  
  # estimated FDR and True Positif
  estimated.FDR_LL = length(which(candidates <= 350))/length(candidates)
  estimated.TP_LL = length(which(candidates > 350))/50
  print(paste("K = 6 FDR:", estimated.FDR_LL, "True Positive:", estimated.TP_LL))
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ PRECIPITATION and TEMPERATURE

#to load project, use:
project_PT <- load.lfmmProject("LFMM_in_precip_temp_all.lfmmProject")

# show the project
show(project_PT)

# summary of the project
summary(project_PT)

# get the zscores of each run for K 
zs_4_PT = z.scores(project_PT, d = 1, K = 4)
zs_5_PT = z.scores(project_PT, d = 1, K = 5)
zs_6_PT = z.scores(project_PT, d = 1 ,K = 6)

zs_4_PT_2 = z.scores(project_PT, d = 2, K = 4)
zs_5_PT_2 = z.scores(project_PT, d = 2, K = 5)
zs_6_PT_2 = z.scores(project_PT, d = 2, K = 6)


# Combine the z-scores using the Stouffer method
zs_4.stouffer_PT = apply(zs_4_PT, MARGIN = 1, median)
zs_5.stouffer_PT = apply(zs_5_PT, MARGIN = 1, median)
zs_6.stouffer_PT = apply(zs_6_PT, MARGIN = 1, median)

zs_4.stouffer_PT_2 = apply(zs_4_PT, MARGIN = 1, median)
zs_5.stouffer_PT_2 = apply(zs_5_PT, MARGIN = 1, median)
zs_6.stouffer_PT_2 = apply(zs_6_PT, MARGIN = 1, median)

# calculate the inflation factor
lambda_4_PT = median(zs_4.stouffer_PT^2)/.456
lambda_5_PT = median(zs_5.stouffer_PT^2)/.456
lambda_6_PT = median(zs_6.stouffer_PT^2)/.456

lambda_4_PT_2 = median(zs_4.stouffer_PT_2^2)/.456
lambda_5_PT_2 = median(zs_5.stouffer_PT_2^2)/.456
lambda_6_PT_2 = median(zs_6.stouffer_PT_2^2)/.456

# calculate adjusted p-values
cp_4.values_PT = pchisq(zs_4.stouffer_PT^2/lambda_4_PT, df = 1, lower = FALSE)
cp_5.values_PT = pchisq(zs_5.stouffer_PT^2/lambda_5_PT, df = 1, lower = FALSE)
cp_6.values_PT = pchisq(zs_6.stouffer_PT^2/lambda_6_PT, df = 1, lower = FALSE)

cp_4.values_PT_2 = pchisq(zs_4.stouffer_PT_2^2/lambda_4_PT_2, df = 1, lower = FALSE)
cp_5.values_PT_2 = pchisq(zs_5.stouffer_PT_2^2/lambda_5_PT_2, df = 1, lower = FALSE)
cp_6.values_PT_2 = pchisq(zs_6.stouffer_PT_2^2/lambda_6_PT_2, df = 1, lower = FALSE)


for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 4 expected FDR:", alpha))
  L = length(cp_4.values_PT)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_4.values_PT) < alpha * (1:L) / L)
  candidates = order(cp_4.values_PT)[w]
  
  # estimated FDR and True Positif
  estimated.FDR = length(which(candidates <= 350))/length(candidates)
  estimated.TP = length(which(candidates > 350))/50
  print(paste("k = 4 FDR:", estimated.FDR, "True Positive:", estimated.TP))
}


for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 5 expected FDR:", alpha))
  L = length(cp_5.values_PT)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_5.values_PT) < alpha * (1:L) / L)
  candidates = order(cp_5.values_PT)[w]
  
  # estimated FDR and True Positif
  estimated.FDR = length(which(candidates <= 350))/length(candidates)
  estimated.TP = length(which(candidates > 350))/50
  print(paste("K = 5 FDR:", estimated.FDR, "True Positive:", estimated.TP))
}

for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 6 expected FDR:", alpha))
  L = length(cp_6.values_PT)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_6.values_PT) < alpha * (1:L) / L)
  candidates = order(cp_6.values_PT)[w]
  
  # estimated FDR and True Positif
  estimated.FDR = length(which(candidates <= 350))/length(candidates)
  estimated.TP = length(which(candidates > 350))/50
  print(paste("K = 6 FDR:", estimated.FDR, "True Positive:", estimated.TP))
}


for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 4 expected FDR:", alpha))
  L = length(cp_4.values_PT_2)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_4.values_PT_2) < alpha * (1:L) / L)
  candidates = order(cp_4.values_PT_2)[w]
  
  # estimated FDR and True Positif
  estimated.FDR = length(which(candidates <= 350))/length(candidates)
  estimated.TP = length(which(candidates > 350))/50
  print(paste("k = 4 FDR:", estimated.FDR, "True Positive:", estimated.TP))
}


for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 5 expected FDR:", alpha))
  L = length(cp_5.values_PT_2)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_5.values_PT_2) < alpha * (1:L) / L)
  candidates = order(cp_5.values_PT_2)[w]
  
  # estimated FDR and True Positif
  estimated.FDR = length(which(candidates <= 350))/length(candidates)
  estimated.TP = length(which(candidates > 350))/50
  print(paste("K = 5 FDR:", estimated.FDR, "True Positive:", estimated.TP))
}

for (alpha in c(.05,.1,.15,.2)) {
  # expected FDR
  print(paste("K = 6 expected FDR:", alpha))
  L = length(cp_6.values_PT_2)
  # return a list of candidates with an expected FDR of alpha.
  w = which(sort(cp_6.values_PT_2) < alpha * (1:L) / L)
  candidates = order(cp_6.values_PT_2)[w]
  
  # estimated FDR and True Positif
  estimated.FDR = length(which(candidates <= 350))/length(candidates)
  estimated.TP = length(which(candidates > 350))/50
  print(paste("K = 6 FDR:", estimated.FDR, "True Positive:", estimated.TP))
}

