### Plotting the PCA from PLINK
###   QC of 16662 overlapping loci between the two data sets 
### Carolyn Tarpey | November 2017 
### ---------------------------------------

library("RColorBrewer")
library(ggplot2)
library(colorspace)
library(plyr)
library(colorRamps)


###############DATA

pop_key <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/POPINFO_updated.txt", header = TRUE, sep = '\t')
pop_key_OLD <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/POPINFO_OLD.txt", header = TRUE, sep = '\t')


ALL_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/All_pops_pca.eigenvec")
ALL_tt = merge(ALL_eigenvec_table, pop_key, by.x = 'V1', by.y = 'POP')
ALL_geo <- ALL_tt[order(ALL_tt$Order_geo),] #sort by a geographical order, then odd then even 

Beringia_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/Beringia_pops_pca.eigenvec")
Beringia_tt = merge(Beringia_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')
Beringia_geo <- Beringia_tt[order(Beringia_tt$Order_geo),] #sort by a geographical order, then odd then even 

NorthAmerica_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/NorthAmerica_pops_pca.eigenvec")
NorthAmerica_tt = merge(NorthAmerica_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')
NorthAmerica_geo <- NorthAmerica_tt[order(NorthAmerica_tt$Order_geo),] #sort by a geographical order, then odd then even 

Odd_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/Odd_pops_pca.eigenvec")
Odd_tt = merge(Odd_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')
odd_geo <- Odd_tt[order(Odd_tt$Order_geo),] #sort by a geographical order, then odd then even 

Even_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/Even_pops_pca.eigenvec")
Even_tt = merge(Even_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')
Even_geo <- Even_tt[order(Even_tt$Order_geo),] #sort by a geographical order, then odd then even 

Old_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/Old_pops_new_genes_pca.eigenvec")
Old_tt = merge(Old_eigenvec_table,pop_key_OLD, by.x = 'V1', by.y = 'POP')
#Old_geo <- Old_tt[order(Old_tt$Order_geo),] #sort by a geographical order, then odd then even 

New_pops_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/New_pops_pca.eigenvec")
New_pops_tt = merge(New_pops_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')
New_pops_geo <- New_pops_tt[order(New_pops_tt$Order_geo),] #sort by a geographical order, then odd then even 

NoSusitna_pops_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/ALL_noSusitna_pca.eigenvec")
NoSusitna_pops_tt = merge(NoSusitna_pops_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')
NoSusitna_pops_geo <- NoSusitna_pops_tt[order(NoSusitna_pops_tt$Order_geo),] #sort by a geographical order, then odd then even 


##########################By population 

AMUR10_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/AMUR10_pca.eigenvec")
AMUR10_tt = merge(AMUR10_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')

AMUR11_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/AMUR11_pca.eigenvec")
AMUR11_tt = merge(AMUR11_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')

HAYLY09_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/HAYLY09_pca.eigenvec")
HAYLY09_tt = merge(HAYLY09_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')

HAYLY10_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/HAYLY10_pca.eigenvec")
HAYLY10_tt = merge(HAYLY10_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')

KOPPE91_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/KOPPE91_pca.eigenvec")
KOPPE91_tt = merge(KOPPE91_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')

KOPPE96_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/KOPPE96_pca.eigenvec")
KOPPE96_tt = merge(KOPPE96_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')

KUSHI06_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/KUSHI06_pca.eigenvec")
KUSHI06_tt = merge(KUSHI06_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')

KUSHI07_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/KUSHI07_pca.eigenvec")
KUSHI07_tt = merge(KUSHI07_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')

LAKEL06_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/LAKEL06_pca.eigenvec")
LAKEL06_tt = merge(LAKEL06_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')

LAKEL07_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/LAKEL07_pca.eigenvec")
LAKEL07_tt = merge(LAKEL07_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')

NOME91_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/NOME91_pca.eigenvec")
NOME91_tt = merge(NOME91_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')

NOME94_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/NOME94_pca.eigenvec")
NOME94_tt = merge(NOME94_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')

SNOH03_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/SNOH03_pca.eigenvec")
SNOH03_tt = merge(SNOH03_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')

SNOH96_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/SNOH96_pca.eigenvec")
SNOH96_tt = merge(SNOH96_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')

SUSIT13_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/SUSIT13_pca.eigenvec")
SUSIT13_tt = merge(SUSIT13_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')

SUSIT14_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/SUSIT14_pca.eigenvec")
SUSIT14_tt = merge(SUSIT14_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')

TAUY09_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/TAUY09_pca.eigenvec")
TAUY09_tt = merge(TAUY09_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')

TAUY12_eigenvec_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/tauy12_edit_pca.eigenvec")
TAUY12_tt = merge(TAUY12_eigenvec_table,pop_key, by.x = 'V1', by.y = 'POP')


names(ALL_tt)

#####COLORS
ALL_col <-c(AMUR10 = "#cd5490", AMUR11 = "#cd5490", SUSIT13= "#563e7f", HAYLY09 = "#c9cf4a", HAYLY10 = "#c9cf4a", 
              KOPPE96 = "#7297ee", KOPPE91 = "#7297ee", KUSHI06 = "#61005e", KUSHI07 = "#61005e", LAKEL06 = "#2171b5",LAKEL07 = "#2171b5",
              NOME91 = "#677f3e", NOME94 = "#677f3e", SNOH03 = "#00158a", SNOH96 ="#00158a", SUSIT14 ="#563e7f",TAUY09 = "#e0957e", TAUY12 = "#e0957e")

ALL_col_noSusitna <-c(AMUR10 = "#cd5490", AMUR11 = "#cd5490", HAYLY09 = "#c9cf4a", HAYLY10 = "#c9cf4a", 
                      KOPPE96 = "#7297ee", KOPPE91 = "#7297ee", KUSHI06 = "#61005e", KUSHI07 = "#61005e", LAKEL06 = "#2171b5",LAKEL07 = "#2171b5",
                      NOME91 = "#677f3e", NOME94 = "#677f3e", SNOH03 = "#00158a", SNOH96 ="#00158a", TAUY09 = "#e0957e", TAUY12 = "#e0957e")

Ber_col <-c(AMUR10 = "#cd5490", AMUR11 = "#cd5490",  HAYLY09 = "#c9cf4a", HAYLY10 = "#c9cf4a", 
             KUSHI06 = "#61005e", KUSHI07 = "#61005e",
            NOME91 = "#677f3e", NOME94 = "#677f3e", TAUY09 = "#e0957e", TAUY12 = "#e0957e")

NA_col <-c(SUSIT13= "#563e7f", 
            KOPPE96 = "#7297ee", KOPPE91 = "#7297ee",LAKEL06 = "#2171b5",LAKEL07 = "#2171b5",
            NOME91 = "#677f3e", NOME94 = "#677f3e", SNOH03 = "#00158a", SNOH96 ="#00158a", SUSIT14 ="#563e7f")

OLD_col <-c(AMUR10 = "#cd5490", AMUR11 = "#cd5490",  HAYLY09 = "#c9cf4a", HAYLY10 = "#c9cf4a", 
            KOPPE96 = "#7297ee", KOPPE91 = "#7297ee", KUSHI06 = "#61005e", KUSHI07 = "#61005e",
            NOME91 = "#677f3e", NOME94 = "#677f3e", SNOH03 = "#00158a", SNOH96 ="#00158a", TAUY09 = "#e0957e", TAUY12 = "#e0957e")

NEW_col <-c(SUSIT13= "#563e7f", LAKEL06 = "#2171b5",LAKEL07 = "#2171b5", SUSIT14 ="#563e7f")

Odd_col <-c(AMUR11 = "#cd5490", SUSIT13= "#563e7f", HAYLY09 = "#c9cf4a", 
            KOPPE91 = "#7297ee", KUSHI07 = "#61005e", LAKEL07 = "#2171b5",
            NOME91 = "#677f3e", SNOH03 = "#00158a", TAUY09 = "#e0957e")

Even_col <-c(AMUR10 = "#cd5490", HAYLY10 = "#c9cf4a", 
            KOPPE96 = "#7297ee",  KUSHI06 = "#61005e",LAKEL06 = "#2171b5",
             NOME94 = "#677f3e", SNOH96 ="#00158a", SUSIT14 ="#563e7f", TAUY12 = "#e0957e")


#ALL_col<- c("#cd5490","#cd5490","#c9cf4a","#c9cf4a","#7297ee","#7297ee","#61005e","#61005e", "#677f3e","#677f3e","#00158a","#00158a","#e0957e","#e0957e","#2171b5","#2171b5","#563e7f", "#563e7f")
#ALL_col_l <- rep(1, length=length(ALL_col)) 
#names(ALL_col_l) <- ALL_col
#pie(ALL_col_l, col=ALL_col, cex=.75, main = "ALL_COL")
# 
# Asia_col<- c("#cd5490","#cd5490","#c9cf4a","#c9cf4a","#61005e","#61005e","#e0957e","#e0957e")
# Asia_col_l <- rep(1, length=length(Asia_col))
# names(Asia_col_l) <- Asia_col
# pie(Asia_col_l, col=Asia_col, cex=.75, main = "Asia_Col")
# 
# Ber_col<- c("#cd5490","#cd5490","#c9cf4a","#c9cf4a","#61005e","#61005e", "#677f3e","#677f3e","#e0957e","#e0957e")
# Ber_col_l <- rep(1, length=length(Ber_col))
# names(Ber_col_l) <- Ber_col
# pie(Ber_col_l, col=Ber_col, cex=.75, main = "Ber_Col")
# 
# NA_col<- c("#7297ee","#7297ee", "#677f3e","#677f3e","#00158a","#00158a")
# NA_col_l <- rep(1, length=length(NA_col))
# names(NA_col_l) <- NA_col
# pie(NA_col_l, col=NA_col, cex=.75, main = "NA_Col")
# 
# NA_ext_col<- c("#7297ee","#7297ee", "#677f3e","#677f3e","#2171b5","#2171b5","#00158a","#00158a","#563e7f", "#563e7f")
# NA_ext_col_l <- rep(1, length=length(NA_ext_col))
# names(NA_ext_col_l) <- NA_ext_col
# pie(NA_ext_col_l, col=NA_ext_col, cex=.75, main = "NA_Ext_Col")
# 


####ALL lineage test for Lakel06 and 07
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/PCAimages/Lakel06_07_lineageTest.pdf", width = 9, height = 7)

ALL_col_test <-c(AMUR10 = "#cd5490", AMUR11 = "#cd5490", SUSIT13= "#563e7f", HAYLY09 = "#c9cf4a", HAYLY10 = "#c9cf4a", 
                 KOPPE96 = "#7297ee", KOPPE91 = "#7297ee", KUSHI06 = "#61005e", KUSHI07 = "#61005e", LAKEL06 = "black", LAKEL07 = "red",
                 NOME91 = "#677f3e", NOME94 = "#677f3e", SNOH03 = "#00158a", SNOH96 ="#00158a", SUSIT14 ="#563e7f",TAUY09 = "#e0957e", TAUY12 = "#e0957e")

ggplot(data = ALL_tt) + geom_point(aes(x = V3, y = V4,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = ALL_col_test, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations TEST ", x = "Principal Component 1", y = "Principal Component 2", size = 30)) 

dev.off()

############Filtering done in plink- these are the 31,485 loci and some filterign: 
###############DATA- these are the results of the filtering done in PLINk

pop_key2 <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/Summary_Stats/POPINFO.txt", header = TRUE, sep = '\t')

ALL_plink_filter <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/PCA/filter_inds_pca.eigenvec")
head(ALL_plink_filter)
ALL_plink_filter_tt = merge(ALL_plink_filter,pop_key2, by.x = 'V1', by.y = 'CLUSTER')
ALL_plink_filter_geo <- ALL_plink_filter_tt[order(ALL_plink_filter_tt$Order_geo),] #sort by a geographical order, then odd then even 
head(ALL_plink_filter_geo)

##PLOT

ggplot(data= ALL_plink_filter_geo) + geom_point(aes(x=V3, y=V4, color= ALL_plink_filter_geo$POPNAME, shape =ALL_plink_filter_geo$LINEAGE), alpha=.8, size=4) +
  theme_classic() + theme(text = element_text(size= 15), legend.title = element_blank(),legend.text = element_text(size= 14),
                          axis.line.x = element_line(), axis.line.y = element_line()) +scale_colour_manual(values = ALL_col_test, name ="Population") 




######################################################################  These are the 16662 loci and some editing of TAUY for inds

###ALL 
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/PCAimages/ALLpops_16662_tauy12edited.pdf", width = 9, height = 7)

ggplot(data = ALL_tt) + geom_point(aes(x = V3, y = V4,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = ALL_col, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Principal Components 1 & 2", x = "Principal Component 1", y = "Principal Component 2", size = 30)) 

dev.off()


##@@@@@@@@@@@@@@@@@@@@@@  one and three

###ALL 

ggplot(data = ALL_tt) + geom_point(aes(x = V3, y = V5,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = ALL_col, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Principal Components 1 & 2", x = "Principal Component 1", y = "Principal Component 3", size = 30)) 

##@@@@@@@@@@@@@@@@@@@@@@ two and three

###ALL 

ggplot(data = ALL_geo) + geom_point(aes(x = V4, y = V5,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = ALL_col, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Principal Components 2 & 3", x = "Principal Component 2 ", y = "Principal Component 3 ", size = 30))

#dev.off()

#####ALL NoSusitna_pops_tt

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/PCAimages/ALLpops_16662_NOSusitna.pdf", width = 9, height = 7)


ggplot(data = NoSusitna_pops_tt) + geom_point(aes(x = V3, y = V4,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = ALL_col_noSusitna, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations, NO Susitna DIM 1 & 2", x = "Principal Component 1", y = "Principal Component 2", size = 30)) 


dev.off()


###################################################

###@@@@@@@@@@@@@@@@@@@@@@dimmension one and two

###OLD 
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/PCAimages/OLD_pops_new_genes.pdf", width = 9, height = 7)

ggplot(data = Old_tt) + geom_point(aes(x = V3, y = V4,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = OLD_col, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "OLD Populations Dimensions 1 & 2", x = "Dimension 1 ", y = "Dimension 2 ", size = 30))

dev.off()


###NEW
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/PCAimages/NEW_pops_pca.pdf", width = 9, height = 7)


ggplot(data = New_pops_tt) + geom_point(aes(x = V3, y = V4,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = ALL_col, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "New Populations Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))

dev.off()

###Beringia_tt

ggplot(data = Beringia_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = Ber_col, name ="Population") +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( title= "Beringia Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))


###NorthAmerica_tt
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/PCAimages/NEW_NA_pops.pdf", width = 9, height = 7)


ggplot(data = NorthAmerica_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_color_manual(values = NA_col, name ="Population") +  scale_x_continuous(trans="reverse")+
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "North America Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30)) 

dev.off()

###@@@@@@@@@@@@@@@@@@@@@@dimmension one and three
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/PCAimages/PCA_plots_different_dimensions.pdf", width = 9, height = 7)

###ALL 

ggplot(data = ALL_geo) + geom_point(aes(x = V3, y = V5,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = ALL_col, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Dimensions 1 & 3", x = "Dimension 1", y = "Dimension 3", size = 30))


###Beringia_tt

ggplot(data = Beringia_geo) + geom_point(aes(x = V3, y = V5, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = Ber_col, name ="Population") +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(), legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "Beringia Dimensions 1 & 3", x = "Dimension 1", y = "Dimension 3", size = 30))


###NorthAmerica_tt

ggplot(data = NorthAmerica_geo) + geom_point(aes(x = V3, y = V5, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_color_manual(values = NA_col, name ="Population") +  scale_x_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "North America Dimensions 1 & 3", x = "Dimension 1", y = "Dimension 3", size = 30)) 


###@@@@@@@@@@@@@@@@@@@@@@dimmension two and three

###ALL 

ggplot(data = ALL_geo) + geom_point(aes(x = V4, y = V5,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = ALL_col, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Dimensions 2 & 3", x = "Dimension 2", y = "Dimension 3", size = 30))


###Beringia_tt

ggplot(data = Beringia_geo) + geom_point(aes(x = V4, y = V5, color = POPNAME), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = Ber_col, name ="Population") +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "Beringia Dimensions 2 & 3", x = "Dimension 2", y = "Dimension 3", size = 30))


###NorthAmerica_tt

ggplot(data = NorthAmerica_geo) + geom_point(aes(x = V4, y = V5, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_color_manual(values = NA_col, name ="Population") +  scale_x_continuous(trans="reverse")+
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "North America Dimensions 2 & 3", x = "Dimension 2", y = "Dimension 3", size = 30)) 



###@@@@@@@@@@@@@@@@@@@@@@dimmension  two and four

###ALL 

ggplot(data = ALL_tt) + geom_point(aes(x = V4, y = V6,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = ALL_col, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Dimensions 2 & 4", x = "Dimension 2", y = "Dimension 4", size = 30))


###Beringia_tt

ggplot(data = Beringia_tt) + geom_point(aes(x = V4, y = V6, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = Ber_col, name ="Population") +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "Beringia Dimensions 2 & 4", x = "Dimension 2", y = "Dimension 4", size = 30))


###NorthAmerica_tt

ggplot(data = NorthAmerica_tt) + geom_point(aes(x = V4, y = V6, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_color_manual(values = NA_col, name ="Population") +  scale_x_continuous(trans="reverse")+
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "North America Dimensions 2 & 4", x = "Dimension 2", y = "Dimension 4", size = 30)) 


###@@@@@@@@@@@@@@@@@@@@@@  dimmension three and four


###ALL 


ggplot(data = ALL_tt) + geom_point(aes(x = V5, y = V6,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = ALL_col, name = "Population") +
  scale_x_continuous(trans = "reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Dimensions 3 & 4", x = "Dimension 3", y = "Dimension 4", size = 30))


###Beringia_tt

ggplot(data = Beringia_tt) + geom_point(aes(x = V5, y = V6, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = Ber_col, name ="Population") +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Beringia Dimensions 3 & 4", x = "Dimension 3", y = "Dimension 4", size = 30))

###NorthAmerica_tt

ggplot(data = NorthAmerica_tt) + geom_point(aes(x = V5, y = V6, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_color_manual(values = NA_col, name ="Population") +  scale_x_continuous(trans="reverse")+
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "North America Dimensions 3 & 4", x = "Dimension 3", y = "Dimension 4", size = 30)) 

dev.off()


###############################  POPULATIONS

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/TEST_PCA/PCAimages/INDV_pops.pdf", width = 9, height = 7)


###AMUR10 

ggplot(data = AMUR10_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "AMUR10 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))


###AMUR11

ggplot(data = AMUR11_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population") +
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "AMUR11 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))


###HAYLY09

ggplot(data = HAYLY09_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population")+ 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "HAYLY09 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))

###HAYLY10

ggplot(data = HAYLY10_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "HAYLY10 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))

###KOPPE91

ggplot(data = KOPPE91_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "KOPPEN91 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))

###KOPPE96_tt

ggplot(data = KOPPE96_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "KOPPEN 96 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))

###LAKEL06_tt

ggplot(data = LAKEL06_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "LAKEL06 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))


###LAKEL07_tt

ggplot(data = LAKEL07_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "LAKEL07 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))

###NOME91_tt

ggplot(data = NOME91_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "NOME91 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))

###NOME94_tt

ggplot(data = NOME94_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "NOME94 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))

###SNOH03_tt

ggplot(data = SNOH03_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "SNOH03 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))


###SNOH96_tt

ggplot(data = SNOH96_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .9, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "SNOH96 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))


###SUSIT13_tt

ggplot(data = SUSIT13_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "SUSIT13 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))

###SUSIT14_tt

ggplot(data = SUSIT14_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "SUSIT14 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))

###TAUY09_tt

ggplot(data = TAUY09_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "TAUY09 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))


###TAUY12_tt

ggplot(data = TAUY12_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .9, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "TAUY12 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))

###KUSHI06_tt

ggplot(data = KUSHI06_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "KUSHI06 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))

###KUSHI07_tt

ggplot(data = KUSHI07_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = ALL_col, name ="Population") + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "KUSHI07 Dimensions 1 & 2", x = "Dimension 1", y = "Dimension 2", size = 30))

dev.off()
