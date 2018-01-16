### Calculating Allelic richness in HierFstat
###    
### Carolyn Tarpey | April 2016
### ---------------------------------------

#install.packages("adegenet")
#install.packages("hierfstat")

library(ape)
library(ggplot2)
library(adegenet)
library(hierfstat)

citation("hierfstat")


setwd('G:/Analysis/Pop_analysis/Populations_b3_may/hierFstat')

#convert the genepop file to genind file type
popdata <- read.genepop("batch_3_16681_pop_80_one_line_genepop.gen")
is.genind(popdata)

popsum<- summary(popdata)

#Do a bartlett test H0: Hexp = Hobs, but this is per locus, not pop and IDK what it means 

#bartlett.test(list(popsum$Hexp, popsum$Hobs))  


#calculating the basic statistics for the pops
b_stats<- basic.stats(popdata)

#figure out how to get them out of R
names(b_stats)
summary(b_stats)

#boxplot of Ho, Hs and Ht
boxplot(basic.stats(popdata)$perloc[,1:3])

#write the basic stats  to a file
write.table(b_stats$Ho, file = "basic_statistics_Ho.txt")
write.table(b_stats$Hs, file = "basic_statistics_Hs.txt")
write.table(b_stats$Fis, file = "basic_statistics_Fis.txt")
write.table(b_stats$perloc, file = "basic_statistics_perloc.txt")
write.table(b_stats$Ho, file = "basic_statistics_Ho.txt")
write.table(b_stats$overall, file = "basic_statistics_overall.txt") 

#calculate the rarefiled alleleic richness per pop 

a_rich<- allelic.richness(popdata)

write.table(a_rich, file = "allelic_richness2_16681.txt")

names(a_rich)

# Amur_even_rich <- mean(a_rich[,2])
# Amur_odd_rich <- mean(a_rich[,3])
# Kamchatka_odd_rich <-mean(a_rich[,4])
# Kamchatka_even_rich <- mean(a_rich[,5])
# PrinceWilliamSound_odd_rich <- mean(a_rich[,6])
# PrinceWilliamSound_even_rich <- mean(a_rich[,7])
# Hokkaido_even_rich <- mean(a_rich[,8])
# Hokkaido_odd_rich <- mean(a_rich[,9])
# NortonSound_odd_rich <- mean(a_rich[,10])
# NortonSound_even_rich <- mean(a_rich[,11])
# PugetSound_odd_rich <- mean(a_rich[,12])
# PugetSound_even_rich <- mean(a_rich[,13])
# Magadan_odd_rich <- mean(a_rich[,14])
# Magadan_even_rich <- mean(a_rich[,15])

Amur_even_rich	<- 1.625065549352910
Amur_odd_rich	<- 1.676718892367610
Kamchatka_odd_rich	<- 1.679049857851170
Kamchatka_even_rich	<- 1.618765778866160
PrinceWilliamSound_odd_rich	<- 1.623227804629000
PrinceWilliamSound_even_rich	<- 1.593153089744230
Hokkaido_even_rich	<- 1.595301561328270
Hokkaido_odd_rich	<- 1.645147883505370
NortonSound_odd_rich	<- 1.678231540796710
NortonSound_even_rich	<- 1.616838333315960
PugetSound_odd_rich	<- 1.555372820616390
PugetSound_even_rich	<- 1.557570128447880
Magadan_odd_rich	<- 1.676561458999230
Magadan_even_rich	<- 1.623081224903950



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





