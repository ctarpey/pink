### Check whether the MCMC has reached convergence
###    using the package coda, and code from: https://github.com/devillemereuil/bayescenv/wiki/5.-How-to-calibrate-the-MCMC
### Carolyn Tarpey | October 2015 
### ---------------------------------------


install.packages("coda")

library(coda)


#@@@ These are the default settings, and the postitively standardized environmental variables. 


#set the working directory 
setwd("G:/Analysis/Pop_analysis/Populations_b3_may/BayeScEnv/ALL_standPos")

### 16681_bayescenv_LA.sel

chain_16681_LA <-read.table("16681_bayescenv_LA.sel",header=TRUE)

chain_16681_LA <-mcmc(chain_16681_LA, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_16681_LA)             #To test for convergence
effectiveSize(chain_16681_LA)           #To compute effective sample size
autocorr.diag(chain_16681_LA)           #To look for auto-correlation
pdf("Plots/16681_bayescenv_LA.pdf", width=8, height=8)

plot(chain_16681_LA)                    #To plot the "trace" and the posterior distribution




x_LA <- pairs(data.frame(chain_16681_LA)) # to test for correlation

dev.off()


### 16681_bayescenv_LO.sel

chain_16681_LO <-read.table("16681_bayescenv_LO.sel",header=TRUE)

chain_16681_LO <-mcmc(chain_16681_LO, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_16681_LO)             #To test for convergence
effectiveSize(chain_16681_LO)           #To compute effective sample size
autocorr.diag(chain_16681_LO)           #To look for auto-correlation
pdf("Plots/16681_bayescenv_LO.pdf", width=8, height=8)

plot(chain_16681_LO)                    #To plot the "trace" and the posterior distribution




x_LO <- pairs(data.frame(chain_16681_LO)) # to test for correlation

dev.off()

  
### 16681_bayescenv_PT1.sel

chain_16681_PT1 <-read.table("16681_bayescenv_PT1.sel",header=TRUE)

chain_16681_PT1 <-mcmc(chain_16681_PT1, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_16681_PT1)             #To test for convergence
effectiveSize(chain_16681_PT1)           #To compute effective sample size
autocorr.diag(chain_16681_PT1)           #To look for auto-correlation
pdf("Plots/16681_bayescenv_PT1.pdf", width=8, height=8)

plot(chain_16681_PT1)                    #To plot the "trace" and the posterior distribution





x_PT1 <- pairs(data.frame(chain_16681_PT1)) # to test for correlation

dev.off()

### 16681_bayescenv_PT2.sel

chain_16681_PT2 <-read.table("16681_bayescenv_PT2.sel",header=TRUE)

chain_16681_PT2 <-mcmc(chain_16681_PT2, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_16681_PT2)             #To test for convergence
effectiveSize(chain_16681_PT2)           #To compute effective sample size
autocorr.diag(chain_16681_PT2)           #To look for auto-correlation
pdf("Plots/16681_bayescenv_PT2.pdf", width=8, height=8)

plot(chain_16681_PT2)                    #To plot the "trace" and the posterior distribution





x_PT2 <- pairs(data.frame(chain_16681_PT2)) # to test for correlation

dev.off()
#_____________________________________________________________________________________________

### NAE_bayescenv_LO.sel

chain_NAE_LO <-read.table("NAE_bayescenv_LO.sel",header=TRUE)
chain_NAE_LO <-mcmc(chain_NAE_LO, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_NAE_LO)             #To test for convergence
effectiveSize(chain_NAE_LO)           #To compute effective sample size
autocorr.diag(chain_NAE_LO)           #To look for auto-correlation
pdf("Plots/NAE_bayescenv_LO.pdf", width=8, height=8)

plot(chain_NAE_LO)                    #To plot the "trace" and the posterior distribution

x_NAE_LO <- pairs(data.frame(chain_NAE_LO)) # to test for correlation

dev.off()
### NAE_bayescenv_LA.sel

chain_NAE_LA <-read.table("NAE_bayescenv_LA.sel",header=TRUE)
chain_NAE_LA <-mcmc(chain_NAE_LA, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_NAE_LA)             #To test for convergence
effectiveSize(chain_NAE_LA)           #To compute effective sample size
autocorr.diag(chain_NAE_LA)           #To look for auto-correlation
pdf("Plots/NAE_bayescenv_LA.pdf", width=8, height=8)

plot(chain_NAE_LA)                    #To plot the "trace" and the posterior distribution

x_NAE_LA <- pairs(data.frame(chain_NAE_LA)) # to test for correlation

dev.off()
### NAE_bayescenv_PT1.sel

chain_NAE_PT1 <-read.table("NAE_bayescenv_PT1.sel",header=TRUE)
chain_NAE_PT1 <-mcmc(chain_NAE_PT1, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_NAE_PT1)             #To test for convergence
effectiveSize(chain_NAE_PT1)           #To compute effective sample size
autocorr.diag(chain_NAE_PT1)           #To look for auto-correlation
pdf("Plots/NAE_bayescenv_PT1.pdf", width=8, height=8)

plot(chain_NAE_PT1)                    #To plot the "trace" and the posterior distribution

x_NAE_PT1 <- pairs(data.frame(chain_NAE_PT1)) # to test for correlation

dev.off()
### NAE_bayescenv_PT2.sel

chain_NAE_PT2 <-read.table("NAE_bayescenv_PT2.sel",header=TRUE)
chain_NAE_PT2 <-mcmc(chain_NAE_PT2, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_NAE_PT2)             #To test for convergence
effectiveSize(chain_NAE_PT2)           #To compute effective sample size
autocorr.diag(chain_NAE_PT2)           #To look for auto-correlation
pdf("Plots/NAE_bayescenv_PT2.pdf", width=8, height=8)

plot(chain_NAE_PT2)                    #To plot the "trace" and the posterior distribution


x_NAE_PT2 <- pairs(data.frame(chain_NAE_PT2)) # to test for correlation

dev.off()
#_____________________________________________________________________________________________

### NAO_bayescenv_LO.sel

chain_NAO_LO <-read.table("NAO_bayescenv_LO.sel",header=TRUE)
chain_NAO_LO <-mcmc(chain_NAO_LO, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_NAO_LO)             #To test for convergence
effectiveSize(chain_NAO_LO)           #To compute effective sample size
autocorr.diag(chain_NAO_LO)           #To look for auto-correlation
pdf("Plots/NAO_bayescenv_LO.pdf", width=8, height=8)

plot(chain_NAO_LO)                    #To plot the "trace" and the posterior distribution

x_NAO_LO <- pairs(data.frame(chain_NAO_LO)) # to test for correlation

dev.off()
### NAO_bayescenv_LA.sel

chain_NAO_LA <-read.table("NAO_bayescenv_LA.sel",header=TRUE)
chain_NAO_LA <-mcmc(chain_NAO_LA, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_NAO_LA)             #To test for convergence
effectiveSize(chain_NAO_LA)           #To compute effective sample size
autocorr.diag(chain_NAO_LA)           #To look for auto-correlation
pdf("Plots/NAO_bayescenv_LA.pdf", width=8, height=8)

plot(chain_NAO_LA)                    #To plot the "trace" and the posterior distribution

x_NAO_LA <- pairs(data.frame(chain_NAO_LA)) # to test for correlation

dev.off()
### NAO_bayescenv_PT1.sel

chain_NAO_PT1 <-read.table("NAO_bayescenv_PT1.sel",header=TRUE)
chain_NAO_PT1 <-mcmc(chain_NAO, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_NAO_PT1)             #To test for convergence
effectiveSize(chain_NAO_PT1)           #To compute effective sample size
autocorr.diag(chain_NAO_PT1)           #To look for auto-correlation
pdf("Plots/NAO_bayescenv_PT1.pdf", width=8, height=8)

plot(chain_NAO_PT1)                    #To plot the "trace" and the posterior distribution

x_NAO_PT1 <- pairs(data.frame(chain_NAO_PT1)) # to test for correlation

dev.off()
### NAO_bayescenv_PT2.sel

chain_NAO_PT2 <-read.table("NAO_bayescenv_PT2.sel",header=TRUE)
chain_NAO_PT2 <-mcmc(chain_NAO_PT2, thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_NAO_PT2)             #To test for convergence
effectiveSize(chain_NAO_PT2)           #To compute effective sample size
autocorr.diag(chain_NAO_PT2)           #To look for auto-correlation
pdf("Plots/NAO_bayescenv_PT2.pdf", width=8, height=8)

plot(chain_NAO_PT2)                    #To plot the "trace" and the posterior distribution

x_NAO_PT2 <- pairs(data.frame(chain_NAO_PT2)) # to test for correlation

dev.off()
#_____________________________________________________________________________________________

### ASE_bayescenv_LA.sel

chain_ASE_LA <-read.table("ASE_bayescenv_LA.sel",header=TRUE)
chain_ASE_LA <-mcmc(chain_ASE_LA,thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_ASE_LA)             #To test for convergence
effectiveSize(chain_ASE_LA)           #To compute effective sample size
autocorr.diag(chain_ASE_LA)           #To look for auto-correlation
pdf("Plots/ASE_bayescenv_LA.pdf", width=8, height=8)

plot(chain_ASE_LA)                    #To plot the "trace" and the posterior distribution

x_ASE_LA <- pairs(data.frame(chain_ASE_LA)) # to test for correlation

dev.off()
### ASE_bayescenv_LO.sel

chain_ASE_LO <-read.table("ASE_bayescenv_LO.sel",header=TRUE)
chain_ASE_LO <-mcmc(chain_ASE_LO ,thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_ASE_LO)             #To test for convergence
effectiveSize(chain_ASE_LO)           #To compute effective sample size
autocorr.diag(chain_ASE_LO)           #To look for auto-correlation
pdf("Plots/ASE_bayescenv_LO.pdf", width=8, height=8)

plot(chain_ASE_LO)                    #To plot the "trace" and the posterior distribution

x_ASE_LO <- pairs(data.frame(chain_ASE_LO)) # to test for correlation

dev.off()
### ASE_bayescenv_PT1.sel

chain_ASE_PT1 <-read.table("ASE_bayescenv_PT1.sel",header=TRUE)
chain_ASE_PT1 <-mcmc(chain_ASE ,thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_ASE_PT1)             #To test for convergence
effectiveSize(chain_ASE_PT1)           #To compute effective sample size
autocorr.diag(chain_ASE_PT1)           #To look for auto-correlation
pdf("Plots/ASE_bayescenv_PT1.pdf", width=8, height=8)

plot(chain_ASE_PT1)                    #To plot the "trace" and the posterior distribution

x_ASE_PT1 <- pairs(data.frame(chain_ASE_PT1)) # to test for correlation

dev.off()
### ASE_bayescenv_PT2.sel

chain_ASE_PT2 <-read.table("ASE_bayescenv_PT2.sel",header=TRUE)
chain_ASE_PT2 <-mcmc(chain_ASE_PT2 ,thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_ASE_PT2)             #To test for convergence
effectiveSize(chain_ASE_PT2)           #To compute effective sample size
autocorr.diag(chain_ASE_PT2)           #To look for auto-correlation
pdf("Plots/ASE_bayescenv_PT2.pdf", width=8, height=8)

plot(chain_ASE_PT2)                    #To plot the "trace" and the posterior distribution

x_ASE_PT2 <- pairs(data.frame(chain_ASE_PT2)) # to test for correlation

dev.off()
#_____________________________________________________________________________________________


### ASO_bayescenv_LA.sel

chain_ASO_LA <-read.table("ASO_bayescenv_LA.sel",header=TRUE)
chain_ASO_LA <-mcmc(chain_ASO_LA ,thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_ASO_LA)             #To test for convergence
effectiveSize(chain_ASO_LA)           #To compute effective sample size
autocorr.diag(chain_ASO_LA)           #To look for auto-correlation
pdf("Plots/ASO_bayescenv_LA.pdf", width=8, height=8)

plot(chain_ASO_LA)                    #To plot the "trace" and the posterior distribution

x_ASO_LA <- pairs(data.frame(chain_ASO_LA)) # to test for correlation

dev.off()
### ASO_bayescenv_LO.sel

chain_ASO_LO <-read.table("ASO_bayescenv_LO.sel",header=TRUE)
chain_ASO_LO <-mcmc(chain_ASO_LO ,thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_ASO_LO)             #To test for convergence
effectiveSize(chain_ASO_LO)           #To compute effective sample size
autocorr.diag(chain_ASO_LO)           #To look for auto-correlation
pdf("Plots/ASO_bayescenv_LO.pdf", width=8, height=8)

plot(chain_ASO_LO)                    #To plot the "trace" and the posterior distribution

x_ASO_LO <- pairs(data.frame(chain_ASO_LO)) # to test for correlation

dev.off()
### ASO_bayescenv_PT1.sel

chain_ASO_PT1 <-read.table("ASO_bayescenv_PT1.sel",header=TRUE)
chain_ASO_PT1 <-mcmc(chain_ASO_PT1 ,thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_ASO_PT1)             #To test for convergence
effectiveSize(chain_ASO_PT1)           #To compute effective sample size
autocorr.diag(chain_ASO_PT1)           #To look for auto-correlation
pdf("Plots/ASO_bayescenv_PT1.pdf", width=8, height=8)

plot(chain_ASO_PT1)                    #To plot the "trace" and the posterior distribution

x_ASO_PT1 <- pairs(data.frame(chain_ASO_PT1)) # to test for correlation

dev.off()
### ASO_bayescenv_PT2.sel

chain_ASO_PT2 <-read.table("ASO_bayescenv_PT2.sel",header=TRUE)
chain_ASO_PT2 <-mcmc(chain_ASO_PT2 ,thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain_ASO_PT2)             #To test for convergence
effectiveSize(chain_ASO_PT2)           #To compute effective sample size
autocorr.diag(chain_ASO_PT2)           #To look for auto-correlation
pdf("Plots/ASO_bayescenv_PT2.pdf", width=8, height=8)

plot(chain_ASO_PT2)                    #To plot the "trace" and the posterior distribution

x_ASO_PT2 <- pairs(data.frame(chain_ASO_PT2)) # to test for correlation

dev.off()
#_____________________________________________________________________________________________







































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

