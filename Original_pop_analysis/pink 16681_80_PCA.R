### Plotting the PCA from PLINK
###   Several variations
### Carolyn Tarpey | July 2015 
### ---------------------------------------



library("RColorBrewer")
library(ggplot2)


###############DATA

pop_key <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/POPINFO.txt", header = TRUE, sep = '\t')
pop_key2 <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/POPINFO2.txt", header = TRUE, sep = '\t')

ALL_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/16681_80_pca.eigenvec")
ALL_tt = merge(ALL_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')

Beringia_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Beringia_pops_pca.eigenvec")
Beringia_tt = merge(Beringia_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')


Asia_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Asia_pops_pca.eigenvec")
Asia_tt = merge(Asia_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')

NorthAmerica_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/NorthAmerica_pops_pca.eigenvec")
NorthAmerica_tt = merge(NorthAmerica_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')


Beringia_Odd_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Beringia_Odd_pops_pca.eigenvec")
Beringia_Odd_tt = merge(Beringia_Odd_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')

Beringia_Even_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Beringia_Even_pops_pca.eigenvec")
Beringia_Even_tt = merge(Beringia_Even_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')

Odd_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Odd_pops_pca.eigenvec")
Odd_tt = merge(Odd_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')

Even_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Even_pops_pca.eigenvec")
Even_tt = merge(Even_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')


names(ALL_tt)

#####COLORS

cols <-c(AMUR10 = "#008000", AMUR11 = "#008000", HAYLY09 = "#10FF10", HAYLY10 = "#10FF10", 
         KOPPE96 = "#FF20B5", KOPPE91 = "#FF20B5", KUSHI06 = "#1A330D", KUSHI07 = "#1A330D",
         NOME91 = "#F1B6DA", NOME94 = "#F1B6DA", SNOH03 = "#990066", SNOH96 = "#990066", TAUY09 = "#3CB371", TAUY12 = "#3CB371")

lins <-c(Even = "#FF9933", Odd = "#CC66FF")
lins_NEW <-c(Even = "#D73C5A", Odd = "#4C7DE7")


blue_pink <-c(AMUR10 = #D73C5A, AMUR11 = "#4E70C8", HAYLY09 = "#6E9D19", HAYLY10 = "#A12787", 
                KOPPE96 = "#A46CD0", KOPPE91 = "#6BDE40", KUSHI06 = "#AD2868", KUSHI07 = "#465C8C",
              NOME91 = "#F1B6DA", NOME94 = "#451B29", SNOH03 = "#990066", SNOH96 = "#E63D25", TAUY09 = "#3E4758", TAUY12 = "#DF9E39")

cols_new <-c(AMUR10 = "#008000", AMUR11 = "#008000", HAYLY09 = "#10FF10", HAYLY10 = "#10FF10", 
             KOPPE96 = "#FF20B5", KOPPE91 = "#FF20B5", KUSHI06 = "#1A330D", KUSHI07 = "#1A330D",
             NOME91 = "#A46CD0", NOME94 = "#A46CD0", SNOH03 = "#990066", SNOH96 = "#990066", TAUY09 = "#3CB371", TAUY12 = "#3CB371")



###TESTERS title = "Individual based PCA with 16,681 SNPs",
ggplot(data = ALL_tt) + geom_point(aes(x = V3, y = V4, shape = LINEAGE, color = POPNAME), alpha = .8, size = 4) + 
  scale_colour_manual(values = cols_OddEven, name ="Population", labels = PopulationNames) +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))


ggplot(data = ALL_tt) + geom_point(aes(x = V3, y = V4, color = LINEAGE), alpha = .8, size = 4) +  scale_colour_manual(values = lins_NEW) + theme_classic()


cols_OddEven <-c(AMUR10 = "#EE547D", AMUR11 = "#4E70C8", HAYLY09 = "#4C7DE7", HAYLY10 = "#ED2DBE", 
                 KOPPE96 = "#D186CF", KOPPE91 = "#B0B0D8", KUSHI06 = "#AD2868", KUSHI07 = "#354A7C",
                 NOME91 = "#95E0CA", NOME94 = "#E63D25", SNOH03 = "#708FEC", SNOH96 ="#a12787" , TAUY09 = "#3C92A8", TAUY12 = "#DF9E39")

cols_OddEven_BER <-c(AMUR10 = "#D73C5A", AMUR11 = "#4E70C8", HAYLY09 = "#4C7DE7", HAYLY10 = "#D24ED2", 
                  KUSHI06 = "#AD2868", KUSHI07 = "#354A7C", NOME91 = "#95E0CA", NOME94 = "#E63D25" , TAUY09 = "#3C92A8", TAUY12 = "#DF9E39")



###@@@@@@@@@@@@@@@@@@@@@@ FOR AFS
###ALL 

ggplot(data = ALL_tt) + geom_point(aes(x = V3, y = V4, shape = CONTINENT, color = POPNAME), alpha = .8, size = 4) + 
  scale_colour_manual(values = cols_OddEven, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))

###ALL 

ggplot(data = ALL_tt) + geom_point(aes(x = V3, y = V4,  color = POPNAME), alpha = .8, size = 4) + 
  scale_colour_manual(values = cols_OddEven, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))


###Beringia_tt
ggplot(data = Beringia_tt) + geom_point(aes(x = V3, y = V4, shape = CONTINENT, color = POPNAME), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = cols_OddEven, name ="Population") +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))

ggplot(data = Beringia_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = cols_OddEven, name ="Population") +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))


###Asia_tt

ggplot(data = Asia_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))

###NorthAmerica_tt

ggplot(data = NorthAmerica_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +  scale_x_continuous(trans="reverse")+
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))

###Beringia_Even_tt

ggplot(data = Beringia_Even_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))

###Beringia_ODD_tt

ggplot(data = Beringia_Odd_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))

###Even_tt

ggplot(data = Even_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .9, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))

###Odd_tt

ggplot(data = Odd_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))



++++++++++++++++++++++++++  ALL   ++++++++++++++++++++++++++
  
ALL_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/16681_80_pca.eigenvec")
ALL_tt = merge(ALL_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')

ggplot(data = ALL_tt) + geom_point(aes(x = V3, y = V4,color = POPNAME), alpha = .5, size = 4) +  scale_color_discrete()
  
ggplot(data = ALL_tt) + geom_point(aes(x = V3, y = V4, shape =LINEAGE, color = POPNAME), alpha = .5, size = 4) + scale_colour_manual(values = cols)
ggplot(data = ALL_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .5, size = 4) +  scale_colour_manual(values = cols)

ggplot(data = ALL_tt) + geom_point(aes(x = V3, y = V4, color = LINEAGE), alpha = .5, size = 4) +  scale_colour_manual(values = lins)

ggplot(data = ALL_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .5, size = 4)+  scale_colour_manual(values = blue_pink)
ggplot(data = ALL_tt) + geom_point(aes(x = V3, y = V4, color = LINEAGE), alpha = .5, size = 4)+  scale_colour_manual(values = lins)

++++++++++++++++++++++++++  Beringia   ++++++++++++++++++++++++++
  
Beringia_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Beringia_pops_pca.eigenvec")
Beringia_tt = merge(Beringia_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')

ggplot(data = Beringia_tt) + geom_point(aes(x = V3, y = V4, shape= LINEAGE, color = POPNAME), alpha = .5, size = 4) +  scale_colour_manual(values = cols)

ggplot(data = Beringia_tt) + geom_point(aes(x = V3, y = V4,color = LINEAGE), alpha = .5, size = 4)+  scale_colour_manual(values = lins)
ggplot(data = Beringia_tt) + geom_point(aes(x = V3, y = V4,color = LINEAGE), alpha = .5, size = 4)+  scale_colour_manual(values = lins)

ggplot(data = Beringia_tt) + geom_point(aes(x = V4, y = V5,color = POPNAME), alpha = .5, size = 4)+  scale_colour_manual(values = cols)
ggplot(data = Beringia_tt) + geom_point(aes(x = V3, y = V4,color = LINEAGE), alpha = .5, size = 4)+  scale_colour_manual(values = lins)

++++++++++++++++++++++++++  Asia   ++++++++++++++++++++++++++
  
Asia_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Asia_pops_pca.eigenvec")
Asia_tt = merge(Asia_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')

ggplot(data = Asia_tt) + geom_point(aes(x = V3, y = V4,color = POPNAME), alpha = .5, size = 4) +  scale_colour_manual(values = cols)
ggplot(data = Asia_tt) + geom_point(aes(x = V3, y = V4,color = LINEAGE), alpha = .5, size = 4) +  scale_colour_manual(values = lins)

ggplot(data = Asia_tt) + geom_point(aes(x = V4, y = V5,color = POPNAME), alpha = .5, size = 4)+  scale_colour_manual(values = cols)
ggplot(data = Asia_tt) + geom_point(aes(x = V3, y = V4,color = LINEAGE), alpha = .5, size = 4)+  scale_colour_manual(values = lins)

++++++++++++++++++++++++++  NorthAmerica   ++++++++++++++++++++++++++

NorthAmerica_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/NorthAmerica_pops_pca.eigenvec")
NorthAmerica_tt = merge(NorthAmerica_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')

ggplot(data = NorthAmerica_tt) + geom_point(aes(x = V3, y = V4,shape = LINEAGE, color = POPNAME), alpha = .5, size = 4)+  scale_colour_manual(values = cols)

ggplot(data = NorthAmerica_tt) + geom_point(aes(x = V3, y = V4,color = LINEAGE), alpha = .5, size = 4)+  scale_colour_manual(values = lins)

ggplot(data = NorthAmerica_tt) + geom_point(aes(x = V4, y = V5,shape = LINEAGE, color = POPNAME), alpha = .5, size = 4)+  scale_colour_manual(values = cols)
ggplot(data = NorthAmerica_tt) + geom_point(aes(x = V3, y = V4,color = LINEAGE), alpha = .5, size = 4)+  scale_colour_manual(values = lins)


++++++++++++++++++++++++++  Beringia Odd   ++++++++++++++++++++++++++

 
Beringia_Odd_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Beringia_Odd_pops_pca.eigenvec")
Beringia_Odd_tt = merge(Beringia_Odd_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')

ggplot(data = Beringia_Odd_tt) + geom_point(aes(x = V3, y = V4,color = POPNAME), alpha = .5, size = 4) +  scale_colour_manual(values = cols)
ggplot(data = Beringia_Odd_tt) + geom_point(aes(x = V3, y = V4,color = LINEAGE), alpha = .5, size = 4) +  scale_colour_manual(values = lins)

ggplot(data = Beringia_Odd_tt) + geom_point(aes(x = V4, y = V5,color = POPNAME), alpha = .5, size = 4)+  scale_colour_manual(values = cols)
ggplot(data = Beringia_Odd_tt) + geom_point(aes(x = V3, y = V4,color = LINEAGE), alpha = .5, size = 4)+  scale_colour_manual(values = lins)

++++++++++++++++++++++++++  Beringia Even   ++++++++++++++++++++++++++

Beringia_Even_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Beringia_Even_pops_pca.eigenvec")
Beringia_Even_tt = merge(Beringia_Even_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')

ggplot(data = Beringia_Even_tt) + geom_point(aes(x = V3, y = V4,color = POPNAME), alpha = .5, size = 4)+  scale_colour_manual(values = cols)
ggplot(data = Beringia_Even_tt) + geom_point(aes(x = V3, y = V4,color = LINEAGE), alpha = .5, size = 4)+  scale_colour_manual(values = lins)

ggplot(data = Beringia_Even_tt) + geom_point(aes(x = V4, y = V5,color = POPNAME), alpha = .5, size = 4)+  scale_colour_manual(values = cols)
ggplot(data = Beringia_Even_tt) + geom_point(aes(x = V3, y = V4,color = LINEAGE), alpha = .5, size = 4)+  scale_colour_manual(values = lins)

++++++++++++++++++++++++++  ODD   ++++++++++++++++++++++++++

Odd_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Odd_pops_pca.eigenvec")
Odd_tt = merge(Odd_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')

ggplot(data = Odd_tt) + geom_point(aes(x = V3, y = V4,color = POPNAME), alpha = .5, size = 4) +  scale_colour_manual(values = cols)
ggplot(data = Odd_tt) + geom_point(aes(x = V3, y = V4,color = LINEAGE), alpha = .5, size = 4) +  scale_colour_manual(values = lins)

ggplot(data = Odd_tt) + geom_point(aes(x = V4, y = V5,color = POPNAME), alpha = .5, size = 4)+  scale_colour_manual(values = cols)
ggplot(data = Odd_tt) + geom_point(aes(x = V3, y = V4,color = LINEAGE), alpha = .5, size = 4)+  scale_colour_manual(values = lins)

++++++++++++++++++++++++++  Even   ++++++++++++++++++++++++++

Even_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Even_pops_pca.eigenvec")
Even_tt = merge(Even_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')

ggplot(data = Even_tt) + geom_point(aes(x = V3, y = V4,color = POPNAME), alpha = .5, size = 4) +  scale_colour_manual(values = cols)
ggplot(data = Even_tt) + geom_point(aes(x = V3, y = V4,color = LINEAGE), alpha = .5, size = 4) +  scale_colour_manual(values = lins)

ggplot(data = Even_tt) + geom_point(aes(x = V4, y = V5,color = POPNAME), alpha = .5, size = 4)+  scale_colour_manual(values = cols)
ggplot(data = Even_tt) + geom_point(aes(x = V3, y = V4,color = LINEAGE), alpha = .5, size = 4)+  scale_colour_manual(values = lins)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  #BeringiaColors = scale_color_manual(values = "AMUR10", "AMUR11", "HAYLY09","HAYLY10","KUSHI06","KUSHI07","NOME91","NOME94","TAUY09","TAUY12")
  #NorthAmericaColors = scale_color_manual(values = "KOPPE96","KOPPE91","NOME91","NOME94","SNOH03","SNOH96")
  #  +  scale_color_discrete()
  # ggsave()
  # KOPPE91 = "#4C51B7"
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  ##I want HUE colors: RAWMATERIAL
  
  blue_odd <- c("#4E70C8","#6E9D19","#AFE469","#6BDE40","#465C8C","#8BDDCA","#3E4758")

periwinkles: "#75C5E1","#7394AB","#708FEC","#B0B0D8" ,"#8586EC"

pink_even <- c("#D73C5A","#A12787","#A46CD0","#AD2868","#451B29","#E63D25","#DF9E39")

salmons: "#EE547D",  "#F15854", "#D85EE6"



#PopulationNames <- c("AMUR 2010", "AMUR 2011", "HAYLYLUYA 2009", "HAYLYLUYA 2010", "KOPPEN 1996", "KOPPEN 1991", "KUSHIRO 2006", "KUSHIRO 2007",
#                     "NOME 1991", "NOME 1994", "SNOHOMISH 2003","SNOHOMISH 1996", "TAUY 2009", "TAUY 2012")

PopulationNames <- c("Amur 2010", "Amur 2011", "Haylyluya 2009", "Haylyluya 2010", "Koppen 1996", "Koppen 1991", "Kushiro 2006", "Kushiro 2007",
                     "Nome 1991", "Nome 1994", "Snohomish 2003","Snohomish 1996", "Tauy 2009", "Tauy 2012")

PopulationNames_BER <- c("Amur 2010", "Amur 2011", "Haylyluya 2009", "Haylyluya 2010", "Kushiro 2006", "Kushiro 2007",
                         "Nome 1991", "Nome 1994", "Tauy 2009", "Tauy 2012")


PopulationNames_BER_odd <- c("Amur 2011", "Haylyluya 2009",  "Kushiro 2007","Nome 1991", "Tauy 2009")


PopulationNames_BER_even <- c("Amur 2010", "Haylyluya 2010", "Kushiro 2006", "Nome 1994",  "Tauy 2012")



PopulationNames_ASIA <- c("Amur 2010", "Amur 2011", "Haylyluya 2009", "Haylyluya 2010", "Kushiro 2006", "Kushiro 2007",
                          "Tauy 2009", "Tauy 2012")

PopulationNames_NA <- c("Koppen 1996", "Koppen 1991", "Nome 1991", "Nome 1994", "Snohomish 2003","Snohomish 1996")

PopulationNames_ODD <- c( "Amur 2011", "Haylyluya 2009", "Koppen 1991",  "Kushiro 2007", "Nome 1991", "Snohomish 2003","Tauy 2009")

PopulationNames_EVEN <- c("Amur 2010", "Haylyluya 2010", "Koppen 1996","Kushiro 2006", "Nome 1994", "Snohomish 1996", "Tauy 2012")

