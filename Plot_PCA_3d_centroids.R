### Plotting PCA in 3D
###    USing the centroid of each pop
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
#install.packages("plot3D")
library(plot3D)
citation("scatterplot3d")

setwd('G:/Analysis/Pop_analysis/Populations_b3_may/3dPCA')

pop_key <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/POPINFO_original.txt", header = TRUE, sep = '\t')

ALL_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/16681_80_pca.eigenvec")
ALL_tt = merge(ALL_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')

pop_1<-ALL_tt[ALL_tt$V1=="pop_1",]
pop_2<-ALL_tt[ALL_tt$V1=="pop_2",]
pop_3<-ALL_tt[ALL_tt$V1=="pop_3",]
pop_4<-ALL_tt[ALL_tt$V1=="pop_4",]
pop_5<-ALL_tt[ALL_tt$V1=="pop_5",]
pop_6<-ALL_tt[ALL_tt$V1=="pop_6",]
pop_7<-ALL_tt[ALL_tt$V1=="pop_7",]
pop_8<-ALL_tt[ALL_tt$V1=="pop_8",]
pop_9<-ALL_tt[ALL_tt$V1=="pop_9",]
pop_10<-ALL_tt[ALL_tt$V1=="pop_10",]
pop_11<-ALL_tt[ALL_tt$V1=="pop_11",]
pop_12<-ALL_tt[ALL_tt$V1=="pop_12",]
pop_13<-ALL_tt[ALL_tt$V1=="pop_13",]
pop_14<-ALL_tt[ALL_tt$V1=="pop_14",]

pop_1_PCA1<-mean(pop_1$V3) 
pop_2_PCA1<-mean(pop_2$V3) 
pop_3_PCA1<-mean(pop_3$V3) 
pop_4_PCA1<-mean(pop_4$V3) 
pop_5_PCA1<-mean(pop_5$V3) 
pop_6_PCA1<-mean(pop_6$V3) 
pop_7_PCA1<-mean(pop_7$V3) 
pop_8_PCA1<-mean(pop_8$V3) 
pop_9_PCA1<-mean(pop_9$V3) 
pop_10_PCA1<-mean(pop_10$V3) 
pop_11_PCA1<-mean(pop_11$V3) 
pop_12_PCA1<-mean(pop_12$V3) 
pop_13_PCA1<-mean(pop_13$V3) 
pop_14_PCA1<-mean(pop_14$V3) 

pop_1_PCA2<-mean(pop_1$V4) 
pop_2_PCA2<-mean(pop_2$V4) 
pop_3_PCA2<-mean(pop_3$V4) 
pop_4_PCA2<-mean(pop_4$V4) 
pop_5_PCA2<-mean(pop_5$V4) 
pop_6_PCA2<-mean(pop_6$V4) 
pop_7_PCA2<-mean(pop_7$V4) 
pop_8_PCA2<-mean(pop_8$V4) 
pop_9_PCA2<-mean(pop_9$V4) 
pop_10_PCA2<-mean(pop_10$V4) 
pop_11_PCA2<-mean(pop_11$V4) 
pop_12_PCA2<-mean(pop_12$V4) 
pop_13_PCA2<-mean(pop_13$V4) 
pop_14_PCA2<-mean(pop_14$V4)  

pop_1_PCA3<-mean(pop_1$V5) 
pop_2_PCA3<-mean(pop_2$V5) 
pop_3_PCA3<-mean(pop_3$V5) 
pop_4_PCA3<-mean(pop_4$V5) 
pop_5_PCA3<-mean(pop_5$V5) 
pop_6_PCA3<-mean(pop_6$V5) 
pop_7_PCA3<-mean(pop_7$V5) 
pop_8_PCA3<-mean(pop_8$V5) 
pop_9_PCA3<-mean(pop_9$V5) 
pop_10_PCA3<-mean(pop_10$V5) 
pop_11_PCA3<-mean(pop_11$V5) 
pop_12_PCA3<-mean(pop_12$V5) 
pop_13_PCA3<-mean(pop_13$V5) 
pop_14_PCA3<-mean(pop_14$V5) 

ALL_cents_PCA1<- c(pop_1_PCA1,pop_2_PCA1,pop_3_PCA1, pop_4_PCA1,pop_5_PCA1, pop_6_PCA1,pop_7_PCA1, 
                   pop_8_PCA1,pop_9_PCA1, pop_10_PCA1, pop_11_PCA1, pop_12_PCA1, pop_13_PCA1, pop_14_PCA1) 

ALL_cents_PCA2<- c(pop_1_PCA2,pop_2_PCA2,pop_3_PCA2, pop_4_PCA2,pop_5_PCA2, pop_6_PCA2,pop_7_PCA2, 
                   pop_8_PCA2,pop_9_PCA2, pop_10_PCA2, pop_11_PCA2, pop_12_PCA2, pop_13_PCA2, pop_14_PCA2) 

ALL_cents_PCA3<- c(pop_1_PCA3,pop_2_PCA3,pop_3_PCA3, pop_4_PCA3,pop_5_PCA3, pop_6_PCA3,pop_7_PCA3, 
                   pop_8_PCA3,pop_9_PCA3, pop_10_PCA3, pop_11_PCA3, pop_12_PCA3, pop_13_PCA3, pop_14_PCA3) 

Beringia_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Beringia_pops_pca.eigenvec")
Beringia_tt = merge(Beringia_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')


pop_1<-Beringia_tt[Beringia_tt$V1=="pop_1",]
pop_2<-Beringia_tt[Beringia_tt$V1=="pop_2",]
pop_3<-Beringia_tt[Beringia_tt$V1=="pop_3",]
pop_4<-Beringia_tt[Beringia_tt$V1=="pop_4",]
pop_7<-Beringia_tt[Beringia_tt$V1=="pop_7",]
pop_8<-Beringia_tt[Beringia_tt$V1=="pop_8",]
pop_9<-Beringia_tt[Beringia_tt$V1=="pop_9",]
pop_10<-Beringia_tt[Beringia_tt$V1=="pop_10",]
pop_13<-Beringia_tt[Beringia_tt$V1=="pop_13",]
pop_14<-Beringia_tt[Beringia_tt$V1=="pop_14",]

pop_1_PCA1<-mean(pop_1$V3) 
pop_2_PCA1<-mean(pop_2$V3) 
pop_3_PCA1<-mean(pop_3$V3) 
pop_4_PCA1<-mean(pop_4$V3) 
pop_7_PCA1<-mean(pop_7$V3) 
pop_8_PCA1<-mean(pop_8$V3) 
pop_9_PCA1<-mean(pop_9$V3) 
pop_10_PCA1<-mean(pop_10$V3) 
pop_13_PCA1<-mean(pop_13$V3) 
pop_14_PCA1<-mean(pop_14$V3) 

pop_1_PCA2<-mean(pop_1$V4) 
pop_2_PCA2<-mean(pop_2$V4) 
pop_3_PCA2<-mean(pop_3$V4) 
pop_4_PCA2<-mean(pop_4$V4) 
pop_7_PCA2<-mean(pop_7$V4) 
pop_8_PCA2<-mean(pop_8$V4) 
pop_9_PCA2<-mean(pop_9$V4) 
pop_10_PCA2<-mean(pop_10$V4) 
pop_13_PCA2<-mean(pop_13$V4) 
pop_14_PCA2<-mean(pop_14$V4)  

pop_1_PCA3<-mean(pop_1$V5) 
pop_2_PCA3<-mean(pop_2$V5) 
pop_3_PCA3<-mean(pop_3$V5) 
pop_4_PCA3<-mean(pop_4$V5) 
pop_7_PCA3<-mean(pop_7$V5) 
pop_8_PCA3<-mean(pop_8$V5) 
pop_9_PCA3<-mean(pop_9$V5) 
pop_10_PCA3<-mean(pop_10$V5) 
pop_13_PCA3<-mean(pop_13$V5) 
pop_14_PCA3<-mean(pop_14$V5) 

Ber_cents_PCA1<- c(pop_1_PCA1,pop_2_PCA1,pop_3_PCA1, pop_4_PCA1, pop_7_PCA1, 
                   pop_8_PCA1,pop_9_PCA1, pop_10_PCA1, pop_13_PCA1, pop_14_PCA1) 

Ber_cents_PCA2<- c(pop_1_PCA2,pop_2_PCA2,pop_3_PCA2, pop_4_PCA2, pop_7_PCA2, 
                   pop_8_PCA2,pop_9_PCA2, pop_10_PCA2, pop_13_PCA2, pop_14_PCA2) 

Ber_cents_PCA3<- c(pop_1_PCA3,pop_2_PCA3,pop_3_PCA3, pop_4_PCA3, pop_7_PCA3, 
                   pop_8_PCA3,pop_9_PCA3, pop_10_PCA3, pop_13_PCA3, pop_14_PCA3) 

Asia_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Asia_pops_pca.eigenvec")
Asia_tt = merge(Asia_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')

pop_1<-Asia_tt[Asia_tt$V1=="pop_1",]
pop_2<-Asia_tt[Asia_tt$V1=="pop_2",]
pop_3<-Asia_tt[Asia_tt$V1=="pop_3",]
pop_4<-Asia_tt[Asia_tt$V1=="pop_4",]
pop_7<-Asia_tt[Asia_tt$V1=="pop_7",]
pop_8<-Asia_tt[Asia_tt$V1=="pop_8",]
pop_13<-Asia_tt[Asia_tt$V1=="pop_13",]
pop_14<-Asia_tt[Asia_tt$V1=="pop_14",]

pop_1_PCA1<-mean(pop_1$V3) 
pop_2_PCA1<-mean(pop_2$V3) 
pop_3_PCA1<-mean(pop_3$V3) 
pop_4_PCA1<-mean(pop_4$V3) 
pop_7_PCA1<-mean(pop_7$V3) 
pop_8_PCA1<-mean(pop_8$V3) 
pop_13_PCA1<-mean(pop_13$V3) 
pop_14_PCA1<-mean(pop_14$V3) 

pop_1_PCA2<-mean(pop_1$V4) 
pop_2_PCA2<-mean(pop_2$V4) 
pop_3_PCA2<-mean(pop_3$V4) 
pop_4_PCA2<-mean(pop_4$V4) 
pop_7_PCA2<-mean(pop_7$V4) 
pop_8_PCA2<-mean(pop_8$V4) 
pop_13_PCA2<-mean(pop_13$V4) 
pop_14_PCA2<-mean(pop_14$V4)  

pop_1_PCA3<-mean(pop_1$V5) 
pop_2_PCA3<-mean(pop_2$V5) 
pop_3_PCA3<-mean(pop_3$V5) 
pop_4_PCA3<-mean(pop_4$V5) 
pop_7_PCA3<-mean(pop_7$V5) 
pop_8_PCA3<-mean(pop_8$V5) 
pop_13_PCA3<-mean(pop_13$V5) 
pop_14_PCA3<-mean(pop_14$V5) 

Asia_cents_PCA1<- c(pop_1_PCA1,pop_2_PCA1,pop_3_PCA1, pop_4_PCA1,pop_7_PCA1, 
                   pop_8_PCA1,pop_13_PCA1, pop_14_PCA1) 

Asia_cents_PCA2<- c(pop_1_PCA2,pop_2_PCA2,pop_3_PCA2, pop_4_PCA2,pop_7_PCA2, 
                   pop_8_PCA2,pop_13_PCA2, pop_14_PCA2) 

Asia_cents_PCA3<- c(pop_1_PCA3,pop_2_PCA3,pop_3_PCA3, pop_4_PCA3,pop_7_PCA3, 
                   pop_8_PCA3,pop_13_PCA3, pop_14_PCA3) 

NorthAmerica_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/NorthAmerica_pops_pca.eigenvec")
NorthAmerica_tt = merge(NorthAmerica_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')

pop_5<-NorthAmerica_tt[NorthAmerica_tt$V1=="pop_5",]
pop_6<-NorthAmerica_tt[NorthAmerica_tt$V1=="pop_6",]
pop_9<-NorthAmerica_tt[NorthAmerica_tt$V1=="pop_9",]
pop_10<-NorthAmerica_tt[NorthAmerica_tt$V1=="pop_10",]
pop_11<-NorthAmerica_tt[NorthAmerica_tt$V1=="pop_11",]
pop_12<-NorthAmerica_tt[NorthAmerica_tt$V1=="pop_12",]

pop_5_PCA1<-mean(pop_5$V3) 
pop_6_PCA1<-mean(pop_6$V3) 
pop_9_PCA1<-mean(pop_9$V3) 
pop_10_PCA1<-mean(pop_10$V3) 
pop_11_PCA1<-mean(pop_11$V3) 
pop_12_PCA1<-mean(pop_12$V3) 

pop_5_PCA2<-mean(pop_5$V4) 
pop_6_PCA2<-mean(pop_6$V4) 
pop_9_PCA2<-mean(pop_9$V4) 
pop_10_PCA2<-mean(pop_10$V4) 
pop_11_PCA2<-mean(pop_11$V4) 
pop_12_PCA2<-mean(pop_12$V4) 

pop_5_PCA3<-mean(pop_5$V5) 
pop_6_PCA3<-mean(pop_6$V5) 
pop_9_PCA3<-mean(pop_9$V5) 
pop_10_PCA3<-mean(pop_10$V5) 
pop_11_PCA3<-mean(pop_11$V5) 
pop_12_PCA3<-mean(pop_12$V5) 

NA_cents_PCA1<- c(pop_5_PCA1, pop_6_PCA1, pop_9_PCA1, pop_10_PCA1, pop_11_PCA1, pop_12_PCA1) 
NA_cents_PCA2<- c(pop_5_PCA2, pop_6_PCA2, pop_9_PCA2, pop_10_PCA2, pop_11_PCA2, pop_12_PCA2) 
NA_cents_PCA3<- c(pop_5_PCA3, pop_6_PCA3, pop_9_PCA3, pop_10_PCA3, pop_11_PCA3, pop_12_PCA3) 

### shapes
ALL_shape<- c(19, 17, 17, 19, 17, 19, 19, 17, 17, 19, 17, 19, 17, 19) #works for the ALL pops
Asia_shape<- c(19, 17, 17, 19, 19, 17, 17, 19) 
Ber_shape<- c(17, 19, 19, 17, 17, 19, 19, 17, 19, 17) 
NA_shape<- c(17, 19, 17, 19, 17, 19) 

#####COLORS

ALL_col<- c("#cd5490","#cd5490","#c9cf4a","#c9cf4a","#7297ee","#7297ee","#61005e","#61005e", "#677f3e","#677f3e","#00158a","#00158a","#e0957e","#e0957e")
Asia_col<- c("#cd5490","#cd5490","#c9cf4a","#c9cf4a","#61005e","#61005e","#e0957e","#e0957e")
Ber_col<- c("#cd5490","#cd5490","#c9cf4a","#c9cf4a","#61005e","#61005e", "#677f3e","#677f3e","#e0957e","#e0957e")
NA_col<- c("#7297ee","#7297ee", "#677f3e","#677f3e","#00158a","#00158a")

# hokkaiod dk purple "#61005e"
# amur pink "#cd5490",
# magadan salmon "#e0957e",
# kamch light green "#c9cf4a",
# norton green "#677f3e",
# pws sky blue "#7297ee",
# puget navy "#00158a",

# pop_1	PAMUR10	"#cd5490"
# pop_2	PAMUR11	"#cd5490"
# pop_3	PHAYLY09 "#c9cf4a"
# pop_4	PHAYLY10 "#c9cf4a"
# pop_5	PKOPPE91 "#7297ee"
# pop_6	PKOPPE96 "#7297ee"
# pop_7	PKUSHI06 "#61005e"
# pop_8	PKUSHI07 "#61005e"
# pop_9	PNOME91	 "#677f3e"
# pop_10	PNOME94	 "#677f3e"
# pop_11	PSNOH03	"#00158a"
# pop_12	PSNOH96	"#00158a"
# pop_13	PTAUY09 "#e0957e"
# pop_14	PTAUY12	"#e0957e"


###@@@@@@@@@@@@@@@@@@@@@@  PLOT 3D

###ALL 

scatterplot3d(ALL_cents_PCA1, ALL_cents_PCA2, ALL_cents_PCA3, pch=ALL_shape, color = ALL_col, cex.symbols= 2, main= "PCA All Populations (3D)", 
              xlab= "Dimension 1 (30.24%)", ylab= "Dimension 2 (13.83%)", zlab= "Dimension 3 (8.73%)")

plot3d(ALL_cents_PCA1, ALL_cents_PCA2, ALL_cents_PCA3, type= "p", main= "PCA All Populations (3D)", xlab= "Dimension 1 (30.24%)", 
       ylab= "Dimension 2 (13.83%)", zlab= "Dimension 3 (8.73%)",  pch= ALL_shape, col = ALL_col, size = 20 )


###Beringia_tt
scatterplot3d(Ber_cents_PCA1, Ber_cents_PCA2, Ber_cents_PCA3, pch=Ber_shape, color = Ber_col, cex.symbols= 2, main="PCA Beringia Populations (3D)",
  xlab= "Dimension 1 (33.75%)", ylab= "Dimension 2 (7.16%)", zlab= "Dimension 3 (6.18%)")

plot3d(Ber_cents_PCA1, Ber_cents_PCA2, Ber_cents_PCA3, type= "p", main="PCA Beringia Populations (3D)",
       xlab= "Dimension 1 (33.75%)", ylab= "Dimension 2 (7.16%)", zlab= "Dimension 3 (6.18%)",  pch= Ber_shape, col = Ber_col, size = 20 )

###Asia_tt
scatterplot3d(Asia_cents_PCA1, Asia_cents_PCA2, Asia_cents_PCA3, pch=Asia_shape, color = Asia_col, cex.symbols= 2, main="PCA Asian Populations (3D)",
  xlab= "Dimension 1 (32.10%)", ylab= "Dimension 2 (6.83%)", zlab= "Dimension 3 (6.04%)")

plot3d(Asia_cents_PCA1, Asia_cents_PCA2, Asia_cents_PCA3, type= "p", main= "PCA Asian Populations (3D)", 
       xlab= "Dimension 1 (32.10%)", ylab= "Dimension 2 (6.83%)", zlab= "Dimension 3 (6.04%)",  pch= Asia_shape, col = Asia_col, size = 20 )


###NorthAmerica_tt
scatterplot3d(NA_cents_PCA1, NA_cents_PCA2, NA_cents_PCA3, pch=NA_shape, color = NA_col,  cex.symbols= 2, main="PCA North American Populations (3D)",
  xlab= "Dimension 1 (22.76%)", ylab= "Dimension 2 (10.36%)", zlab= "Dimension 3 (6.73%)")

plot3d(NA_cents_PCA1, NA_cents_PCA2, NA_cents_PCA3, type= "p", main="PCA North American Populations (3D)",
  xlab= "Dimension 1 (22.76%)", ylab= "Dimension 2 (10.36%)", zlab= "Dimension 3 (6.73%)",  pch= NA_shape, col = NA_col, size = 20)




