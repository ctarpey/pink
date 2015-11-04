### Check whether the MCMC has reached convergence
###    using the package coda, and code from: https://github.com/devillemereuil/bayescenv/wiki/5.-How-to-calibrate-the-MCMC
### Carolyn Tarpey | October 2015 
### ---------------------------------------


install.packages("coda")

library(coda)


#@@@@ these are the orginal run of the bayescenv that are all the defaults 


#set the working directory 
setwd("G:/Analysis/Pop_analysis/Populations_b3_may/BayeScEnv")

### 16681_80_bayescan.sel

chain_16681 <-read.table("16681_80_bayescan.sel",header=TRUE)

chain_16681 <-mcmc(chain_16681, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_16681)             #To test for convergence
effectiveSize(chain_16681)           #To compute effective sample size
autocorr.diag(chain_16681)           #To look for auto-correlation
plot(chain_16681)                    #To plot the "trace" and the posterior distribution

x <- pairs(data.frame(chain_16681)) # to test for correlation


### NA_even_bayescan.sel

chain_NAE <-read.table("NA_even_bayescan.sel",header=TRUE)
chain_NAE <-mcmc(chain_NAE, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_NAE)             #To test for convergence
effectiveSize(chain_NAE)           #To compute effective sample size
autocorr.diag(chain_NAE)           #To look for auto-correlation
plot(chain_NAE)                    #To plot the "trace" and the posterior distribution

x_NAE <- pairs(data.frame(chain_NAE)) # to test for correlation
 

### NA_odd_bayescan.sel

chain_NAO <-read.table("NA_odd_bayescan.sel",header=TRUE)
chain_NAO <-mcmc(chain_NAO, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_NAO)             #To test for convergence
effectiveSize(chain_NAO)           #To compute effective sample size
autocorr.diag(chain_NAO)           #To look for auto-correlation
plot(chain_NAO)                    #To plot the "trace" and the posterior distribution

x_NAO <- pairs(data.frame(chain_NAO)) # to test for correlation


### asia_even_bayescan.sel

chain_ASE <-read.table("asia_even_bayescan.sel",header=TRUE)
chain_ASE <-mcmc(chain_ASE ,thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_ASE)             #To test for convergence
effectiveSize(chain_ASE)           #To compute effective sample size
autocorr.diag(chain_ASE)           #To look for auto-correlation
plot(chain_ASE)                    #To plot the "trace" and the posterior distribution

x_ASE <- pairs(data.frame(chain_ASE)) # to test for correlation

### asia_odd_bayescan.sel

chain_ASO <-read.table("asia_odd_bayescan.sel",header=TRUE)
chain_ASO <-mcmc(chain_ASO ,thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_ASO)             #To test for convergence
effectiveSize(chain_ASO)           #To compute effective sample size
autocorr.diag(chain_ASO)           #To look for auto-correlation
plot(chain_ASO)                    #To plot the "trace" and the posterior distribution

x_ASO <- pairs(data.frame(chain_ASO)) # to test for correlation




#@@ These runs are from the DOUBLE set that i ran taht have an extra 0 on the end of the # of iterations


#set the working directory 
setwd("G:/Analysis/Pop_analysis/Populations_b3_may/BayeScEnv/Bayescenv_output")

### 16681_80_bayescan.sel

chain_16681_d <-read.table("16681_80_bayescan.sel",header=TRUE)

chain_16681_d <-mcmc(chain_16681_d, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_16681_d)             #To test for convergence
effectiveSize(chain_16681_d)           #To compute effective sample size
autocorr.diag(chain_16681_d)           #To look for auto-correlation
plot(chain_16681_d)                    #To plot the "trace" and the posterior distribution

x_d <- pairs(data.frame(chain_16681_d)) # to test for correlation


### NA_even_bayescan.sel

chain_NAE_d <-read.table("NA_even_bayescan.sel",header=TRUE)
chain_NAE_d <-mcmc(chain_NAE_d, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_NAE_d)             #To test for convergence
effectiveSize(chain_NAE_d)           #To compute effective sample size
autocorr.diag(chain_NAE_d)           #To look for auto-correlation
plot(chain_NAE_d)                    #To plot the "trace" and the posterior distribution

x_NAE_d <- pairs(data.frame(chain_NAE_d)) # to test for correlation


### NA_odd_bayescan.sel

chain_NAO_d <-read.table("NA_odd_bayescan.sel",header=TRUE)
chain_NAO_d <-mcmc(chain_NAO_d, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_NAO_d)             #To test for convergence
effectiveSize(chain_NAO_d)           #To compute effective sample size
autocorr.diag(chain_NAO_d)           #To look for auto-correlation
plot(chain_NAO_d)                    #To plot the "trace" and the posterior distribution

x_NAO_d <- pairs(data.frame(chain_NAO_d)) # to test for correlation


### asia_even_bayescan.sel

chain_ASE_d <-read.table("asia_even_bayescan.sel",header=TRUE)
chain_ASE_d <-mcmc(chain_ASE_d ,thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_ASE_d)             #To test for convergence
effectiveSize(chain_ASE_d)           #To compute effective sample size
autocorr.diag(chain_ASE_d)           #To look for auto-correlation
plot(chain_ASE_d)                    #To plot the "trace" and the posterior distribution

x_ASE_d <- pairs(data.frame(chain_ASE_d)) # to test for correlation

### asia_odd_bayescan.sel

chain_ASO_d <-read.table("asia_odd_bayescan.sel",header=TRUE)
chain_ASO_d <-mcmc(chain_ASO_d ,thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_ASO_d)             #To test for convergence
effectiveSize(chain_ASO_d)           #To compute effective sample size
autocorr.diag(chain_ASO_d)           #To look for auto-correlation
plot(chain_ASO_d)                    #To plot the "trace" and the posterior distribution

x_ASO_d <- pairs(data.frame(chain_ASO_d)) # to test for correlation

