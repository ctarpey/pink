### Testing differences in Heterozygosity
###   Between lineages 
###    Basic population indices, from genepop file, calc in excel
### Carolyn Tarpey | May 2016
### ---------------------------------------

setwd('G:/Analysis/Pop_analysis/Populations_b3_may/Genepop/observedHet')

##########################################

gen_het<- read.table("genhet.txt", header = TRUE)


head(gen_het)
Amur_even_het <- mean(gen_het[gen_het$pop=="Amur_even",6])
Amur_odd_het <- mean(gen_het[gen_het$pop=="Amur_odd",6])
Kamchatka_odd_het <- mean(gen_het[gen_het$pop=="Kamchatka_odd",6])
Kamchatka_even_het <- mean(gen_het[gen_het$pop=="Kamchatka_even",6])
PrinceWilliamSound_odd_het <- mean(gen_het[gen_het$pop=="PrinceWilliamSound_odd",6])
PrinceWilliamSound_even_het <- mean(gen_het[gen_het$pop=="PrinceWilliamSound_even",6])
Hokkaido_even_het <- mean(gen_het[gen_het$pop=="Hokkaido_even",6])
Hokkaido_odd_het <- mean(gen_het[gen_het$pop=="Hokkaido_odd",6])
NortonSound_odd_het <- mean(gen_het[gen_het$pop=="NortonSound_odd",6])
NortonSound_even_het <- mean(gen_het[gen_het$pop=="NortonSound_even",6])
PugetSound_odd_het <- mean(gen_het[gen_het$pop=="PugetSound_odd",6])
PugetSound_even_het <- mean(gen_het[gen_het$pop=="PugetSound_even",6])
Magadan_odd_het <- mean(gen_het[gen_het$pop=="Magadan_odd",6])
Magadan_even_het <- mean(gen_het[gen_het$pop=="Magadan_even",6])

even <- c(Amur_even_het,Kamchatka_even_het,PrinceWilliamSound_even_het,Hokkaido_even_het,
          NortonSound_even_het,PugetSound_even_het,Magadan_even_het)
odd <- c(Amur_odd_het,Kamchatka_odd_het,PrinceWilliamSound_odd_het,Hokkaido_odd_het,
         NortonSound_odd_het,PugetSound_odd_het,Magadan_odd_het)

median(even)
median(odd)

mean(even)
mean(odd)

wilcox.test(even, odd, alternative = "two.sided", exact = TRUE )
wilcox.test(even, odd, alternative = "g", exact = TRUE )
wilcox.test(even, odd, alternative = "l", exact = TRUE )
wilcox.test(odd, even, alternative = "g", exact = TRUE )






###########################################################
By_locus <- read.table("by_locus.txt", header = TRUE)

head(By_locus)



Amur_even_het <- mean(By_locus[By_locus$Pop=="1",4])
Amur_odd_het <- mean(By_locus[By_locus$Pop=="2",4])
Kamchatka_odd_het <- mean(By_locus[By_locus$Pop=="3",4])
Kamchatka_even_het <- mean(By_locus[By_locus$Pop=="4",4])
PrinceWilliamSound_odd_het <- mean(By_locus[By_locus$Pop=="5",4])
PrinceWilliamSound_even_het <- mean(By_locus[By_locus$Pop=="6",4])
Hokkaido_even_het <- mean(By_locus[By_locus$Pop=="7",4])
Hokkaido_odd_het <- mean(By_locus[By_locus$Pop=="8",4])
NortonSound_odd_het <- mean(By_locus[By_locus$Pop=="9",4])
NortonSound_even_het <- mean(By_locus[By_locus$Pop=="10",4])
PugetSound_odd_het <- mean(By_locus[By_locus$Pop=="11",4])
PugetSound_even_het <- mean(By_locus[By_locus$Pop=="12",4])
Magadan_odd_het <- mean(By_locus[By_locus$Pop=="13",4])
Magadan_even_het <- mean(By_locus[By_locus$Pop=="14",4])


Amur_even_het <- mean(By_locus[By_locus$Pop=="1",4])
Amur_odd_het <- mean(By_locus[By_locus$Pop=="2",4])
Kamchatka_odd_het <- mean(By_locus[By_locus$Pop=="3",4])
Kamchatka_even_het <- mean(By_locus[By_locus$Pop=="4",4])
PrinceWilliamSound_odd_het <- mean(By_locus[By_locus$Pop=="5",4])
PrinceWilliamSound_even_het <- mean(By_locus[By_locus$Pop=="6",4])
Hokkaido_even_het <- mean(By_locus[By_locus$Pop=="7",4])
Hokkaido_odd_het <- mean(By_locus[By_locus$Pop=="8",4])
NortonSound_odd_het <- mean(By_locus[By_locus$Pop=="9",4])
NortonSound_even_het <- mean(By_locus[By_locus$Pop=="10",4])
PugetSound_odd_het <- mean(By_locus[By_locus$Pop=="11",4])
PugetSound_even_het <- mean(By_locus[By_locus$Pop=="12",4])
Magadan_odd_het <- mean(By_locus[By_locus$Pop=="13",4])
Magadan_even_het <- mean(By_locus[By_locus$Pop=="14",4])

even_loc <- c(Amur_even_het,Kamchatka_even_het,PrinceWilliamSound_even_het,Hokkaido_even_het,
          NortonSound_even_het,PugetSound_even_het,Magadan_even_het)
odd_loc <- c(Amur_odd_het,Kamchatka_odd_het,PrinceWilliamSound_odd_het,Hokkaido_odd_het,
         NortonSound_odd_het,PugetSound_odd_het,Magadan_odd_het)

median(even_loc)
median(odd_loc)

mean(even_loc)
mean(odd_loc)

wilcox.test(even_loc, odd_loc, alternative = "two.sided", exact = TRUE )
wilcox.test(even_loc, odd_loc, alternative = "g", exact = TRUE )
wilcox.test(even_loc, odd_loc, alternative = "l", exact = TRUE )
wilcox.test(odd_loc, even_loc, alternative = "g", exact = TRUE )



#################################
#average expected het

even_He <-c(0.160,0.163,0.161,0.161,0.160,0.160,0.154)

odd_He<- c(0.170,0.172,0.171,0.171,0.172,0.165,0.156)

mean(even_He)
mean(odd_He)

even_Ne <- c(965,20731,6673,56866,16395,3990,1110)
odd_Ne <- c(1699,7051,20321,18847,5143,2565,5916)

mean(even_Ne)
mean(odd_Ne)




##########################################

All_data<- read.table("All_Data.txt", header = TRUE)
even_pops<- read.table("even_pops.txt", header = TRUE)
odd_pops <- read.table("odd_pops.txt", header = TRUE)

head(All_data)
Amur_even_het <- mean(All_data[All_data$pop=="Amur_even",4])
Amur_odd_het <- mean(All_data[All_data$pop=="Amur_odd",4])
Kamchatka_odd_het <- mean(All_data[All_data$pop=="Kamchatka_odd",4])
Kamchatka_even_het <- mean(All_data[All_data$pop=="Kamchatka_even",4])
PrinceWilliamSound_odd_het <- mean(All_data[All_data$pop=="PrinceWilliamSound_odd",4])
PrinceWilliamSound_even_het <- mean(All_data[All_data$pop=="PrinceWilliamSound_even",4])
Hokkaido_even_het <- mean(All_data[All_data$pop=="Hokkaido_even",4])
Hokkaido_odd_het <- mean(All_data[All_data$pop=="Hokkaido_odd",4])
NortonSound_odd_het <- mean(All_data[All_data$pop=="NortonSound_odd",4])
NortonSound_even_het <- mean(All_data[All_data$pop=="NortonSound_even",4])
PugetSound_odd_het <- mean(All_data[All_data$pop=="PugetSound_odd",4])
PugetSound_even_het <- mean(All_data[All_data$pop=="PugetSound_even",4])
Magadan_odd_het <- mean(All_data[All_data$pop=="Magadan_odd",4])
Magadan_even_het <- mean(All_data[All_data$pop=="Magadan_even",4])

even <- c(Amur_even_het,Kamchatka_even_het,PrinceWilliamSound_even_het,Hokkaido_even_het,
          NortonSound_even_het,PugetSound_even_het,Magadan_even_het)
odd <- c(Amur_odd_het,Kamchatka_odd_het,PrinceWilliamSound_odd_het,Hokkaido_odd_het,
         NortonSound_odd_het,PugetSound_odd_het,Magadan_odd_het)

median(even)
median(odd)

mean(even)
mean(odd)

wilcox.test(even, odd, alternative = "two.sided", exact = TRUE )
wilcox.test(even, odd, alternative = "g", exact = TRUE )
wilcox.test(even, odd, alternative = "l", exact = TRUE )
wilcox.test(odd, even, alternative = "g", exact = TRUE )

Ber_even <- c(Amur_even_het,Kamchatka_even_het,Hokkaido_even_het,
          NortonSound_even_het,Magadan_even_het)
N_even <- c(PrinceWilliamSound_even_het,PugetSound_even_het)

wilcox.test(Ber_even, N_even, alternative = "two.sided", exact = TRUE )
wilcox.test(Ber_even, N_even, alternative = "g", exact = TRUE )
wilcox.test(Ber_even, N_even, alternative = "l", exact = TRUE )


Ber_odd <- c(Amur_odd_het,Kamchatka_odd_het,Hokkaido_odd_het,
         NortonSound_odd_het,Magadan_odd_het)
N_odd <- c(PrinceWilliamSound_odd_het,PugetSound_odd_het)

wilcox.test(Ber_odd, N_odd, alternative = "two.sided", exact = TRUE )
wilcox.test(Ber_odd, N_odd, alternative = "g", exact = TRUE )
wilcox.test(Ber_odd, N_odd, alternative = "l", exact = TRUE )
