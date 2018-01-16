### Plotting PCA in 3D
###    
### Carolyn Tarpey | May 2016
### ---------------------------------------

#install.packages("scatterplot3d")

library(scatterplot3d)
library("RColorBrewer")
library(ggplot2)
library(colorspace)
library(plyr)
library(colorRamps)
library(rgl)

citation("scatterplot3d")

setwd('G:/Analysis/Pop_analysis/Populations_b3_may/3dPCA')

pop_key <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/POPINFO_original.txt", header = TRUE, sep = '\t')

ALL_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/16681_80_pca.eigenvec")
ALL_tt = merge(ALL_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')
ALL_geo <- ALL_tt[order(ALL_tt$Order_geo),] #sort by a geographical order, then odd then even 
ALL_lin <- ALL_tt[order(ALL_tt$Order_lin),] #sort by odd then even then geographical order 

Beringia_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Beringia_pops_pca.eigenvec")
Beringia_tt = merge(Beringia_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')
Beringia_geo <- Beringia_tt[order(Beringia_tt$Order_geo),] #sort by a geographical order, then odd then even 
Beringia_lin <- Beringia_tt[order(Beringia_tt$Order_lin),]  #sort by odd then even then geographical order  

Asia_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Asia_pops_pca.eigenvec")
Asia_tt = merge(Asia_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')
Asia_geo <- Asia_tt[order(Asia_tt$Order_geo),] #sort by a geographical order, then odd then even 
Asia_lin<- Asia_tt[order(Asia_tt$Order_lin),]  #sort by odd then even then geographical order  

NorthAmerica_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/NorthAmerica_pops_pca.eigenvec")
NorthAmerica_tt = merge(NorthAmerica_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')
NorthAmerica_geo <- NorthAmerica_tt[order(NorthAmerica_tt$Order_geo),] #sort by a geographical order, then odd then even 
NorthAmerica_lin <- NorthAmerica_tt[order(NorthAmerica_tt$Order_lin),]  #sort by odd then even then geographical order  

names(ALL_geo)
head(ALL_geo)
head(ALL_lin)

#####COLORS
cols_OddEven <-c(AMUR10 = "#EE547D", AMUR11 = "#4E70C8", HAYLY09 = "#4C7DE7", HAYLY10 = "#ED2DBE", 
                 KOPPE96 = "#D186CF", KOPPE91 = "#B0B0D8", KUSHI06 = "#AD2868", KUSHI07 = "#354A7C",
                 NOME91 = "#95E0CA", NOME94 = "#E63D25", SNOH03 = "#708FEC", SNOH96 ="#a12787" , TAUY09 = "#3C92A8", TAUY12 = "#DF9E39")

#pal<- choose_palette()

## These are pretty, but the mid colors are too light 
cols_ord7 <- c(AMUR10 = "#fc8d59", AMUR11 ="#fc8d59",  HAYLY09 = "#ffffbf", HAYLY10 = "#ffffbf", 
               KOPPE96 = "#91bfdb", KOPPE91 = "#91bfdb", KUSHI06 = "#d73027", KUSHI07 = "#d73027",
               NOME91 = "#e0f3f8", NOME94 = "#e0f3f8", SNOH03 = "#4575b4", SNOH96 ="#4575b4" , TAUY09 = "#fee090", TAUY12 = "#fee090")
# #d73027  #fc8d59  #fee090  #ffffbf  #e0f3f8  #91bfdb  #4575b4

## These are harsh and still not good colorRamps::blue2red(7) 
#cols_ord7 <- c(AMUR10 = "#0055FF", AMUR11 ="#0055FF",  HAYLY09 = "#00FFFF", HAYLY10 = "#00FFFF", 
#               KOPPE96 = "#FF5500", KOPPE91 = "#FF5500", KUSHI06 = "#0000FF", KUSHI07 = "#0000FF",
#               NOME91 = "#FFAA00", NOME94 = "#FFAA00", SNOH03 = "#FF0000", SNOH96 ="#FF0000" , TAUY09 = "#00AAFF", TAUY12 = "#00AAFF")

#"#0000FF" "#0055FF" "#00AAFF" "#00FFFF" "#FFAA00" "#FF5500" "#FF0000"

## These are too close together colorRampPalette(c("red", "blue"))(7)
# cols_ord7 <- c(AMUR10 = "#D4002A", AMUR11 ="#D4002A",  HAYLY09 = "#7F007F", HAYLY10 = "#7F007F", 
#                KOPPE96 = "#2A00D4", KOPPE91 = "#2A00D4", KUSHI06 = "#FF0000", KUSHI07 = "#FF0000",
#                NOME91 = "#5500AA", NOME94 = "#5500AA", SNOH03 = "#0000FF", SNOH96 ="#0000FF" , TAUY09 = "#AA0055", TAUY12 = "#AA0055")
#"#FF0000" "#D4002A" "#AA0055" "#7F007F" "#5500AA" "#2A00D4" "#0000FF"

#scale_shape_discrete(name ="Population", breaks=c("AMUR10", "AMUR11", "HAYLY09", "HAYLY10", 
#  "KOPPE96", "KOPPE91", "KUSHI06", "KUSHI07","NOME91", "NOME94", "SNOH03", "SNOH96", "TAUY09", "TAUY12"),

#labels=c("Amur even", "Amur odd", "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", 
#  "Prince William Sound even", "Hokkaido even", "Hokkaido odd", "Norton Sound odd", "Norton Sound even", 
#  "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")

### shapes
shapes_OddEven <- c(6, 1, 1, 6, 6, 1, 6, 1, 1, 6, 1, 6, 1, 6) #works for the ALL pops

###@@@@@@@@@@@@@@@@@@@@@@  PLOT 3D
#pdf("G:/Analysis/Pop_analysis/Populations_b3_may/PCA_R_images_plink/PCA_in_3D.pdf", width = 9, height = 7)

cloud(z~x+y, data = DF, pch= 19, col.point = DF$group, 
      key = list(points = list(pch = 19, col = seq_along(levels(DF$group))), 
                 text = list(levels(DF$group)), space = 'top', columns = nlevels(DF$group)))


# plot
with(df, scatterplot3d(PC1, PC2, PC3, color = as.numeric(POPULATION), pch=19, main="PCA 1000 Genomes Exome Data (3D)")) 

# add legend
legend("topleft", pch=19, col=c(AMUR10 = "#0055FF", AMUR11 ="#0055FF",  HAYLY09 = "#00FFFF", HAYLY10 = "#00FFFF", 
              KOPPE96 = "#FF5500", KOPPE91 = "#FF5500", KUSHI06 = "#0000FF", KUSHI07 = "#0000FF",
               NOME91 = "#FFAA00", NOME94 = "#FFAA00", SNOH03 = "#FF0000", SNOH96 ="#FF0000" , TAUY09 = "#00AAFF", TAUY12 = "#00AAFF"), 
       legend=c("Amur even", "Amur odd", "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
              "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even"))

# plot
with(ALL_geo, scatterplot3d(x = V3, y = V4, z = V5,  color = as.numeric(POPNAME), pch = as.numeric(LINEAGE), main="PCA All Populations (3D)"),
     main= "PCA All Populations (3D)", xlab= "Dimension 1", ylab= "Dimension 2", zlab= "Dimension 3") 
legend("topleft", pch=19, col=c(AMUR10 = "#0055FF", AMUR11 ="#0055FF",  HAYLY09 = "#00FFFF", HAYLY10 = "#00FFFF", 
                                KOPPE96 = "#FF5500", KOPPE91 = "#FF5500", KUSHI06 = "#0000FF", KUSHI07 = "#0000FF",
                                NOME91 = "#FFAA00", NOME94 = "#FFAA00", SNOH03 = "#FF0000", SNOH96 ="#FF0000" , TAUY09 = "#00AAFF", TAUY12 = "#00AAFF"), legend=c("AFR", "EUR", "SAS", "EAS"))

plot3d(ALL_geo[,3:5], main= "PCA All Populations (3D)", xlab= "Dimension 1", ylab= "Dimension 2", zlab= "Dimension 3", col= as.numeric(ALL_geo$POPNAME))

###ALL 

ggplot(data = ALL_geo) + geom_point(aes(x = V3, y = V4,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = cols_ord7, name ="Population", labels = c("Amur even", "Amur odd",
    "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
    "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) +
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Dimensions 1 & 2", x = "Dimension 1 (30.24%)", y = "Dimension 2 (13.83%)", size = 30))


with(ALL_geo, scatterplot3d(x = V3, y = V4, z = V5,  color = as.numeric(POPNAME), pch = as.numeric(LINEAGE), main="PCA All Populations (3D)"),
     main= "PCA All Populations (3D)", xlab= "Dimension 1 (30.24%)", ylab= "Dimension 2 (13.83%)", zlab= "Dimension 3 (8.73%)") 

plot3d(ALL_geo[,3:5], main= "PCA All Populations (3D)", xlab= "Dimension 1 (30.24%)", ylab= "Dimension 2 (13.83%)", zlab= "Dimension 3 (8.73%)", col= as.numeric(ALL_geo$POPNAME))


###Beringia_tt

ggplot(data = Beringia_geo) + geom_point(aes(x = V3, y = V4, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = cols_OddEven, name ="Population", labels = c("Amur even", "Amur odd", "Kamchatka odd", "Kamchatka even",
        "Hokkaido even","Hokkaido odd", "Norton Sound odd", "Norton Sound even",  "Magadan odd","Magadan even")) +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( title= "Beringia Dimensions 1 & 2", x = "Dimension 1 (33.75%)", y = "Dimension 2 (7.16%)", size = 30))

with(Beringia_geo, scatterplot3d(x = V3, y = V4, z = V5,  color = as.numeric(POPNAME), pch = as.numeric(LINEAGE), main="PCA Beringia Populations (3D)"),
     xlab= "Dimension 1 (33.75%)", ylab= "Dimension 2 (7.16%)", zlab= "Dimension 3 (6.18%)") 

plot3d(Beringia_geo[,3:5], main= "PCA Beringia Populations (3D)", xlab= "Dimension 1 (33.75%)", ylab= "Dimension 2 (7.16%)", zlab= "Dimension 3 (6.18%)", col= as.numeric(Beringia_geo$POPNAME))


###Asia_tt

ggplot(data = Asia_geo) + geom_point(aes(x = V3, y = V4, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population", labels = c("Amur even", "Amur odd","Kamchatka odd", "Kamchatka even", "Hokkaido even", "Hokkaido odd","Magadan odd","Magadan even")) + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(), legend.text = element_text(size= 17), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Asia Dimensions 1 & 2", x = "Dimension 1 (32.10%)", y = "Dimension 2 (6.83%)", size = 30))


with(Asia_geo, scatterplot3d(x = V3, y = V4, z = V5,  color = as.numeric(POPNAME), pch = as.numeric(LINEAGE), main="PCA Asian Populations (3D)"),
     xlab= "Dimension 1 (32.10%)", ylab= "Dimension 2 (6.83%)", zlab= "Dimension 3 (6.04%)") 

plot3d(Asia_geo[,3:5], main= "PCA Asian Populations (3D)", xlab= "Dimension 1 (32.10%)", ylab= "Dimension 2 (6.83%)", zlab= "Dimension 3 (6.04%)", col= as.numeric(Asia_geo$POPNAME))


 ###NorthAmerica_tt
 
ggplot(data = NorthAmerica_geo) + geom_point(aes(x = V3, y = V4, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population", labels = c("Prince William Sound odd", "Prince William Sound even", 
      "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even")) +  scale_x_continuous(trans="reverse")+
  axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "North America Dimensions 1 & 2", x = "Dimension 1 (22.76%)", y = "Dimension 2 (10.36%)", size = 30)) 
                     

with(NorthAmerica_geo, scatterplot3d(x = V3, y = V4, z = V5,  color = as.numeric(POPNAME), pch = as.numeric(LINEAGE), main="PCA North American Populations (3D)"),
     xlab= "Dimension 1 (32.10%)", ylab= "Dimension 2 (10.36%)", zlab= "Dimension 3 (6.73%)") 

plot3d(NorthAmerica_geo[,3:5], main= "PCA North America Populations (3D)", xlab= "Dimension 1 (22.76%)", ylab= "Dimension 2 (10.36%)", zlab= "Dimension 3 (6.73%)", col= as.numeric(NorthAmerica_geo$POPNAME))


                    