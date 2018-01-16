### Calculating F statistics 
###    nonhierarchically with diveRsity
### Carolyn Tarpey | April 2016
### ---------------------------------------

#install.packages("diveRsity")

library(diveRsity)

citation("diveRsity")


setwd('G:/Analysis/Pop_analysis/Populations_b3_may/diVeRsity')

All <- readGenepop("batch3_ALL_GENEPOP.txt")
Amur <- readGenepop("batch3_Amur_GENEPOP.txt")
Hok <- readGenepop("batch3_Hokkaido_GENEPOP.txt")
Mag <- readGenepop("batch3_Magadan_GENEPOP.txt")
Kam <- readGenepop("batch3_Kamchatka_GENEPOP.txt")
Nor <- readGenepop("batch3_Norton_GENEPOP.txt")
PWS <- readGenepop("batch3_PWS_GENEPOP.txt")
Pug <- readGenepop("batch3_Puget_GENEPOP.txt")
Even <- readGenepop("batch3_Even_GENEPOP.txt")
Odd <- readGenepop("batch3_Odd_GENEPOP.txt")