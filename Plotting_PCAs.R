### Plotting the PCA from PLINK
###   Final versions 
### Carolyn Tarpey | April 2016 
### ---------------------------------------

library("RColorBrewer")
library(ggplot2)



###############DATA

pop_key <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/POPINFO.txt", header = TRUE, sep = '\t')

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

cols_OddEven <-c(AMUR10 = "#EE547D", AMUR11 = "#4E70C8", HAYLY09 = "#4C7DE7", HAYLY10 = "#ED2DBE", 
                 KOPPE96 = "#D186CF", KOPPE91 = "#B0B0D8", KUSHI06 = "#AD2868", KUSHI07 = "#354A7C",
                 NOME91 = "#95E0CA", NOME94 = "#E63D25", SNOH03 = "#708FEC", SNOH96 ="#a12787" , TAUY09 = "#3C92A8", TAUY12 = "#DF9E39")



###@@@@@@@@@@@@@@@@@@@@@@dimmension one and two

###ALL 
pdf("G:/Analysis/Pop_analysis/Populations_b3_may/PCA_R_images_plink/PDF of updated PCA plots.pdf", width = 9, height = 7)

ggplot(data = ALL_tt) + geom_point(aes(x = V3, y = V4,  color = POPNAME), alpha = .8, size = 4) + 
  scale_colour_manual(values = cols_OddEven, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Dimensions 1 & 2", x = "Dimension 1 (30.24%)", y = "Dimension 2 (13.83%)", size = 30))


###Beringia_tt

ggplot(data = Beringia_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = cols_OddEven, name ="Population") +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( title= "Beringia Dimensions 1 & 2", x = "Dimension 1 (33.75%)", y = "Dimension 2 (7.16%)", size = 30))


###Asia_tt

ggplot(data = Asia_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Asia Dimensions 1 & 2", x = "Dimension 1 (32.10%)", y = "Dimension 2 (6.83%)", size = 30))

###NorthAmerica_tt

ggplot(data = NorthAmerica_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +  scale_x_continuous(trans="reverse")+
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "North America Dimensions 1 & 2", x = "Dimension 1 (22.76%)", y = "Dimension 2 (10.36%)", size = 30)) 


###@@@@@@@@@@@@@@@@@@@@@@dimmension one and three

###ALL 

ggplot(data = ALL_tt) + geom_point(aes(x = V3, y = V5,  color = POPNAME), alpha = .8, size = 4) + 
  scale_colour_manual(values = cols_OddEven, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Dimensions 1 & 3", x = "Dimension 1 (30.24%)", y = "Dimension 3 (8.73%)", size = 30))


###Beringia_tt

ggplot(data = Beringia_tt) + geom_point(aes(x = V3, y = V5, color = POPNAME), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = cols_OddEven, name ="Population") +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "Beringia Dimensions 1 & 3", x = "Dimension 1 (33.75%)", y = "Dimension 3 (6.18%)", size = 30))


###Asia_tt

ggplot(data = Asia_tt) + geom_point(aes(x = V3, y = V5, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "Asia Dimensions 1 & 3", x = "Dimension 1 (32.10%)", y = "Dimension 3 (6.04%)", size = 30))

###NorthAmerica_tt

ggplot(data = NorthAmerica_tt) + geom_point(aes(x = V3, y = V5, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +  scale_x_continuous(trans="reverse")+
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "North America Dimensions 1 & 3", x = "Dimension 1 (22.76%)", y = "Dimension 3 (6.73%)", size = 30)) 


###@@@@@@@@@@@@@@@@@@@@@@dimmension two and three

###ALL 

ggplot(data = ALL_tt) + geom_point(aes(x = V4, y = V5,  color = POPNAME), alpha = .8, size = 4) + 
  scale_colour_manual(values = cols_OddEven, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Dimensions 2 & 3", x = "Dimension 2 (13.83%)", y = "Dimension 3 (8.73%)", size = 30))


###Beringia_tt

ggplot(data = Beringia_tt) + geom_point(aes(x = V4, y = V5, color = POPNAME), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = cols_OddEven, name ="Population") +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "Beringia Dimensions 2 & 3", x = "Dimension 2 (7.16%)", y = "Dimension 3 (6.18%)", size = 30))


###Asia_tt

ggplot(data = Asia_tt) + geom_point(aes(x = V4, y = V5, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "Asia Dimensions 2 & 3", x = "Dimension 2 (6.83%)", y = "Dimension 3 (6.04%)", size = 30))

###NorthAmerica_tt

ggplot(data = NorthAmerica_tt) + geom_point(aes(x = V4, y = V5, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +  scale_x_continuous(trans="reverse")+
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "North America Dimensions 2 & 3", x = "Dimension 2 (10.36%)", y = "Dimension 3 (6.73%)", size = 30)) 



###@@@@@@@@@@@@@@@@@@@@@@dimmension  two and four

###ALL 

ggplot(data = ALL_tt) + geom_point(aes(x = V4, y = V6,  color = POPNAME), alpha = .8, size = 4) + 
  scale_colour_manual(values = cols_OddEven, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Dimensions 2 & 4", x = "Dimension 2 (13.83%)", y = "Dimension 4 (5.13%)", size = 30))


###Beringia_tt

ggplot(data = Beringia_tt) + geom_point(aes(x = V4, y = V6, color = POPNAME), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = cols_OddEven, name ="Population") +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "Beringia Dimensions 2 & 4", x = "Dimension 2 (7.16%)", y = "Dimension 4 (3.54%)", size = 30))


###Asia_tt

ggplot(data = Asia_tt) + geom_point(aes(x = V4, y = V6, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "Asia Dimensions 2 & 4", x = "Dimension 2 (6.83%)", y = "Dimension 4 (3.40%)", size = 30))

###NorthAmerica_tt

ggplot(data = NorthAmerica_tt) + geom_point(aes(x = V4, y = V6, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +  scale_x_continuous(trans="reverse")+
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "North America Dimensions 2 & 4", x = "Dimension 2 (10.36%)", y = "Dimension 4 (5.88%)", size = 30)) 


###@@@@@@@@@@@@@@@@@@@@@@  dimmension three and four

###ALL 


ggplot(data = ALL_tt) + geom_point(aes(x = V5, y = V6,  color = POPNAME), alpha = .8, size = 4) + 
  scale_colour_manual(values = cols_OddEven, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Dimensions 3 & 4", x = "Dimension 3 (8.73%)", y = "Dimension 4 (5.13%)", size = 30))


###Beringia_tt

ggplot(data = Beringia_tt) + geom_point(aes(x = V5, y = V6, color = POPNAME), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = cols_OddEven, name ="Population") +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Beringia Dimensions 3 & 4", x = "Dimension 3 (6.18%)", y = "Dimension 4 (3.54%)", size = 30))


###Asia_tt

ggplot(data = Asia_tt) + geom_point(aes(x = V5, y = V6, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( x = "Dimension 3 (6.04%)", y = "Dimension 4 (3.40%)", size = 30))

###NorthAmerica_tt

ggplot(data = NorthAmerica_tt) + geom_point(aes(x = V5, y = V6, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +  scale_x_continuous(trans="reverse")+
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "North America Dimensions 3 & 4", x = "Dimension 3 (6.73%)", y = "Dimension 4 (5.88%)", size = 30)) 

dev.off()









###############################  Lineages

###Beringia_Even_tt

ggplot(data = Beringia_Even_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))


###Beringia_ODD_tt

ggplot(data = Beringia_Odd_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))

###Even_tt

ggplot(data = Even_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .9, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))

###Odd_tt

ggplot(data = Odd_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))


# 
# scale_shape_discrete(name ="Population", breaks=c("AMUR10", "AMUR11", "HAYLY09", "HAYLY10", 
#                                                   "KOPPE96", "KOPPE91", "KUSHI06", "KUSHI07","NOME91", "NOME94", "SNOH03", "SNOH96", "TAUY09", "TAUY12"),
#                      labels=c("Amur even", "Amur odd", "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", 
#                               "Prince William Sound even", "Hokkaido even", "Hokkaido odd", "Norton Sound odd", "Norton Sound even", 
#                               "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) +
