### Trace the MCMC of BayeScEnv
###   code taken from https://github.com/devillemereuil/bayescenv/wiki/5.-How-to-calibrate-the-MCMC
### Carolyn Tarpey | September 2015 
### ---------------------------------------



library(coda)
chain<-read.table("out.sel",header=TRUE)
chain<-mcmc(chain,thin=10)     #Adapt thin to its actual value (10 is the default)
heidel.diag(chain)             #To test for convergence
effectiveSize(chain)           #To compute effective sample size
autocorr.diag(chain)           #To look for auto-correlation
plot(chain)                    #To plot the "trace" and the posterior distribution

##These small tests will help you ensure convergence is effective and auto-correlation not too much of a problem. 
##Auto-correlation is acceptable if the effective sample size is large enough (at least above 1000!).

##We recommend effective sizes of at least the Fst's to be reported in Supplementary Material when using the method.