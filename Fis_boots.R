### Calculating Fis in HierFstat
###    want bootstrap confidence intervals per pop
### Carolyn Tarpey | June 2016
### ---------------------------------------


library(ape)
library(ggplot2)
library(adegenet)
library(hierfstat)
library(poppr)


#library(missingno)
#install.packages("poppr")
#install.packages("missingno")
#citation("hierfstat")


setwd('G:/Analysis/Pop_analysis/Populations_b3_may/hierFstat')

#convert the genepop file to genind file type
popdata <- read.genepop("batch_3_16681_pop_80_one_line_genepop.gen")
is.genind(popdata)

popsum<- summary(popdata)
names(popsum)


#calculating the basic statistics for the pops
b_stats<- basic.stats(popdata)

#figure out how to get them out of R
names(b_stats)
summary(b_stats)

b_stats$perloc

#write the basic stats  to a file
write.table(b_stats$Fis, file = "basic_statistics_Fis.txt")


data_AmE<- read.genepop("batch3_AMUR_e_GENEPOP.gen")
data_AmO
data_KaO
data_KaE
data_PwO
data_PwE
data_HoE
data_HoO
data_NoO
data_NoE
data_PuO
data_PuE
data_Mao
data_MaE




#calculate the Fis, with CI or some sig per pop 

#calculate the basic statistics for the population
### this is where there is an error in the NAs



b_stats_AmE <- basic.stats(data_AmE)                      
b_stats_AmO <- basic.stats(data_AmO_m) 
b_stats_KaO <- basic.stats(data_KaO_m) 
b_stats_KaE <- basic.stats(data_KaE_m)
b_stats_PwO <- basic.stats(data_PwO_m) 
b_stats_PwE <- basic.stats(data_PwE_m) 
b_stats_HoE <- basic.stats(data_HoE_m) 
b_stats_HoO <- basic.stats(data_HoO_m) 
b_stats_NoO <- basic.stats(data_NoO_m) 
b_stats_NoE <- basic.stats(data_NoE_m)
b_stats_PuO <- basic.stats(data_PuO_m)
b_stats_PuE <- basic.stats(data_PuE_m) 
b_stats_Mao <- basic.stats(data_Mao_m) 
b_stats_MaE <- basic.stats(data_MaE_m)

#####################################################

# The confidence intervals for the data

boots <- boot.ppfis(popdata)

#######################################################








write.table(a_rich, file = "allelic_richness2_16681.txt")


even_rich <- c(Amur_even_rich,Kamchatka_even_rich,PrinceWilliamSound_even_rich,Hokkaido_even_rich,
               NortonSound_even_rich,PugetSound_even_rich,Magadan_even_rich)
odd_rich <- c(Amur_odd_rich,Kamchatka_odd_rich,PrinceWilliamSound_odd_rich,Hokkaido_odd_rich,
              NortonSound_odd_rich,PugetSound_odd_rich,Magadan_odd_rich)

median(even_rich)
median(odd_rich)

mean(even_rich)
mean(odd_rich)

wilcox.test(even_rich, odd_rich, alternative = "two.sided", exact = TRUE )
wilcox.test(even_rich, odd_rich, alternative = "g", exact = TRUE )
wilcox.test(even_rich, odd_rich, alternative = "l", exact = TRUE )
wilcox.test(odd_rich, even_rich, alternative = "g", exact = TRUE )

Ber_even <- c(Amur_even_rich,Kamchatka_even_rich,Hokkaido_even_rich,
              NortonSound_even_rich,Magadan_even_rich)
N_even <- c(PrinceWilliamSound_even_rich,PugetSound_even_rich)

wilcox.test(Ber_even, N_even, alternative = "two.sided", exact = TRUE )
wilcox.test(Ber_even, N_even, alternative = "g", exact = TRUE )
wilcox.test(Ber_even, N_even, alternative = "l", exact = TRUE )


#######################################################################################
#######################################################################################
#######################################################################################
#######################################################################################

### Break up the entire data set into genind objects for each population
pop_labels <- c(rep("AmE",32),rep("AmO",32),
                rep("KaO",31),rep("KaE",31),rep("PwO",24),rep("PwE",24),
                rep("HoE",32),rep("HoO",32),rep("NoO",24),rep("NoE",18),rep("PuO",24),rep("PuE",24),rep("Mao",28),rep("MaE",27))  
## Creates a vector containing the population assignments of each individual
#numbers in code are numbers of individuals in each population


## Creates a list of genind objects for each population
pops_separated <- seppop(popdata,pop=pop_labels)
names(pops_separated)


#Creates a genind object comprising only the AD individuals
data_AmE <-pops_separated$AmE 

#Verify that the genind object has the correct number of individuals and loci
data_AmE                    

###Repeat for all populations
data_AmO <-pops_separated$AmO  
data_KaO <-pops_separated$KaO
data_KaE <-pops_separated$KaE
data_PwO <-pops_separated$PwO
data_PwE <-pops_separated$PwE
data_HoE <-pops_separated$HoE
data_HoO <-pops_separated$HoO
data_NoO <-pops_separated$NoO
data_NoE <-pops_separated$NoE
data_PuO <-pops_separated$PuO
data_PuE <-pops_separated$PuE
data_Mao <-pops_separated$Mao
data_MaE <-pops_separated$MaE


is.genind(data_AmE)

#remove the loci that are not polymorphic from each population
#That doesnt work. 

data_AmE_m <- missingno(data_AmE, type= "loci", cutoff= 0.10, quiet= FALSE, freq = FALSE) 
data_AmO_m <- missingno(data_AmO, type= "loci", cutoff= 0.10, quiet= FALSE, freq = FALSE) 
data_KaO_m <- missingno(data_KaO, type= "loci", cutoff= 0.10, quiet= FALSE, freq = FALSE) 
data_KaE_m <- missingno(data_KaE, type= "loci", cutoff= 0.10, quiet= FALSE, freq = FALSE) 
data_PwO_m <- missingno(data_PwO, type= "loci", cutoff= 0.10, quiet= FALSE, freq = FALSE) 
data_PwE_m <- missingno(data_PwE, type= "loci", cutoff= 0.10, quiet= FALSE, freq = FALSE) 
data_HoE_m <- missingno(data_HoE, type= "loci", cutoff= 0.10, quiet= FALSE, freq = FALSE) 
data_HoO_m <- missingno(data_HoO, type= "loci", cutoff= 0.10, quiet= FALSE, freq = FALSE) 
data_NoO_m <- missingno(data_NoO, type= "loci", cutoff= 0.10, quiet= FALSE, freq = FALSE) 
data_NoE_m <- missingno(data_NoE, type= "loci", cutoff= 0.10, quiet= FALSE, freq = FALSE) 
data_PuO_m <- missingno(data_PuO, type= "loci", cutoff= 0.10, quiet= FALSE, freq = FALSE) 
data_PuE_m <- missingno(data_PuE, type= "loci", cutoff= 0.10, quiet= FALSE, freq = FALSE) 
data_Mao_m <- missingno(data_Mao, type= "loci", cutoff= 0.10, quiet= FALSE, freq = FALSE) 
data_MaE_m <- missingno(data_MaE, type= "loci", cutoff= 0.10, quiet= FALSE, freq = FALSE) 





