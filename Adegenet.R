### Run Adegenet on the population data
###   for summary stats, AMOVA, dapc
### Carolyn Tarpey | September 2015 
### ---------------------------------------


install.packages("adegenet", dep=TRUE)
install.packages("poppr", dep =TRUE)


library(adegenet)
library(ape)
library(pegas)
library(seqinr)
library(ggplot2)
library(poppr)
library(hierfstat)



##this one works, but the data has not been filtered for individuals
pops_16681 <- read.genepop("G:/Analysis/Pop_analysis/Populations_b3_may/Adegenet/16681_nopop_no80.gen")

#This one works too and is the final data set
pops_16681_edit <- read.genepop("G:/Analysis/Pop_analysis/Populations_b3_may/Adegenet/16681_80_edit.gen")

pops_16681_edit@pop

names(pops_16681_edit)

toto <- summary(pops_16681_edit)

edit(pops_16681_edit)
names(toto)

###plots that arent useful for these data.
#par(mfrow=c(2,2))

#plot(toto$  pop.eff, toto$  pop.nall, xlab=  "Colonies sample size",  ylab=  "Number of alleles", main=  "Alleles numbers and sample sizes",type=  "n")
#text(toto$  pop.eff, toto$  pop.nall, lab=  names(toto$  pop.eff))
#barplot(toto$  loc.n.all,ylab=  "Number of alleles",main=  "Number of alleles per locus")
#barplot(toto$  Hexp-  toto$  Hobs,main= "Heterozygosity: expected-observed",ylab=  "Hexp - Hobs")

## Q: Is the Mean observed H significantly lower than the mean expected H? 

bartlett.test(list(toto$Hexp, toto$Hobs))
t.test(toto$Hexp,toto$Hobs, pair = T, var.equal = TRUE, alter = "greater")

# this doesnt work
#Fst(as.loci(pops_16681_edit))
#Fst <- Fst(as.loci(pops_16681_edit))


#This is taking forever to run 
pairwise.fst(pops_16681_edit, pop= NULL, res.type = c("dist", "matrix"))

#im not really sure what this is, a discriminant analysis of principal components 
dapc.All<- dapc (pops_16681_edit, pops_16681_edit$pop)
??dapc
vignette("adegenet-dapc")


#####################################
###The rest of this code is derived from here and google searches on the internet
##http://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html



scatter(dapc.All, cell = 0, posi.da = "topleft", 
        scree.pca= TRUE, posi.pca= "bottomleft", ratio.da = .2, cex.lab = .5,
        ratio.pca= .2, pch = 10, cstar = 1,  axesell = TRUE, mstree = TRUE, lwd = 2, lty = 2)

names(dapc.All)

set.seed(14)

# change graphical parameter to subsequently overlay the labels withoutdrawing a new plot
par(new=TRUE)
# make a data frame of the dapc coordinates used in scatter
df <- data.frame(x = dapc.All$ind.coord[,1], y = dapc.All$ind.coord[,2])
# identify/ create a vector of names for the individuals in your plot
noms <- paste("ind", c(1:400), sep=".")
# use the text function to add labels to the positions given by the coordinates you used in plot
s.label(dfxy = df, xax=1, yax=2, label=noms,
        clabel=0.7, # change the size of the labels
        boxes=TRUE, # if points are spaced wide enough, can use TRUE to add boxes around the labels
        grid=FALSE, addaxes=FALSE) # do not draw lines or axes in addition to the labels



names(pops_16681_edit)
pops_16681_edit$loc.fac


contrib <- loadingplot(dapc.All$var.contr, axis= 2, thresh= .001,  lab.jitter= 1)

temp <- seploc(pops_16681_edit)

names(temp) #lists all the loci names

dapc.All$var.contr

dapc.df = print.dapc(dapc.All)

write.table(dapc.All, file= "G:/Analysis/Pop_analysis/Populations_b3_may/Adegenet/16681_R_dapc.txt", sep = '\t', eol = "\n", row.names = TRUE, col.names= TRUE)
