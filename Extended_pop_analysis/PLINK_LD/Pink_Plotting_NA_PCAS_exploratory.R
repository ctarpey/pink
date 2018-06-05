### Plotting The North American Extended Populations PCA from PLINK 
###   Data exploration to see what clusters should have LD run together
###
### Carolyn Tarpey | June 2018
### ---------------------------------------

library(RColorBrewer)
library(ggplot2)
library(colorspace)
library(plyr)
library(colorRamps)

###############Import the pop info for the populations (IT IS UPDATED!!)
pop_key <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/POPINFO.txt", header = TRUE, sep = '\t')

###############Import the eigenvalues that were run in PLINK and merge them with the Pop info
NA_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAsofNA/na_all_pca.eigenvec")
NA_eigen_table <- merge(NA_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
NA_eigen_table <- NA_eigen_table[order(NA_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(NA_eigen_table)

EVEN_NA_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAsofNA/even_na_pca.eigenvec")
EVEN_NA_eigen_table <- merge(EVEN_NA_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
EVEN_NA_eigen_table <- EVEN_NA_eigen_table[order(EVEN_NA_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(EVEN_NA_eigen_table)

ODD_NA_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAsofNA/odd_na_pca.eigenvec")
ODD_NA_eigen_table <- merge(ODD_NA_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
ODD_NA_eigen_table <- ODD_NA_eigen_table[order(ODD_NA_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(ODD_NA_eigen_table)


#####COLORS
#cols_OddEven <-c(AMUR10 = "#EE547D", AMUR11 = "#4E70C8", HAYLY09 = "#4C7DE7", HAYLY10 = "#ED2DBE", 
#                 KOPPE96 = "#D186CF", KOPPE91 = "#B0B0D8", KUSHI06 = "#AD2868", KUSHI07 = "#354A7C",
#                 NOME91 = "#95E0CA", NOME94 = "#E63D25", SNOH03 = "#708FEC", SNOH96 ="#a12787" , TAUY09 = "#3C92A8", TAUY12 = "#DF9E39")

cols_OddEven <-c( SUSIT13 = "#4C7DE7", SUSIT14 = "#ED2DBE", LAKEL07 = "#D186CF", LAKEL06 = "#B0B0D8", 
                  KOPPE91 = "#95E0CA", KOPPE96 = "#E63D25", SNOH03 = "#708FEC", SNOH96 ="#a12787" )


NA_col<- c("#cd5490","#cd5490","#c9cf4a","#c9cf4a","#61005e","#61005e","#e0957e","#e0957e")
NA_col_4<- c("#cd5490","#c9cf4a","#61005e","#e0957e")


pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAsofNA/ALLNA_pcas.pdf", width = 9, height = 7)
###ALL one and two
ggplot(data = NA_eigen_table) + geom_point(aes(x = V3, y = V4,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = NA_col, name ="Population") +  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All North American Populations Principal Components 1 & 2", x = "Principal Component 1 (11.51%)", y = "Principal Component 2 (7.38%)", size = 20))

###ALL one and three
ggplot(data = NA_eigen_table) + geom_point(aes(x = V3, y = V5,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = NA_col, name ="Population")+   scale_x_continuous(trans="reverse") +  
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All North American Populations Principal Components 1 & 3", x = "Principal Component 1 (11.51%)", y = "Principal Component 3 (4.48%)", size = 20))

###ALL two and three
ggplot(data = NA_eigen_table) + geom_point(aes(x = V4, y = V5,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = NA_col, name ="Population") +  
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All North American Populations Principal Components 2 & 3", x = "Principal Component 2 (7.38%)", y = "Principal Component 3 (4.48%)", size = 20))
dev.off()


pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAsofNA/EVENNA_pcas.pdf", width = 9, height = 7)
###Even NA one and two
ggplot(data = EVEN_NA_eigen_table) + geom_point(aes(x = V3, y = V4,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = NA_col_4, name ="Population")  + scale_x_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Even North American Populations Principal Components 1 & 2", x = "Principal Component 1 (4.00%)", y = "Principal Component 2 (1.94%)", size = 20))

###Even NA one and three
ggplot(data = EVEN_NA_eigen_table) + geom_point(aes(x = V3, y = V5,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = NA_col_4, name ="Population") + scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Even North American Populations Principal Components 1 & 3", x = "Principal Component 1 (4.00%)", y = "Principal Component 3 (1.29%)", size = 20))

###Even NA two and three
ggplot(data = EVEN_NA_eigen_table) + geom_point(aes(x = V4, y = V5,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = NA_col_4, name ="Population") +  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Even North American Populations Principal Components 2 & 3", x = "Principal Component 2 (1.94%)", y = "Principal Component 3 (1.29%)", size = 20))
dev.off()

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAsofNA/ODDNA_pcas.pdf", width = 9, height = 7)
###odd NA one and two
ggplot(data = ODD_NA_eigen_table) + geom_point(aes(x = V3, y = V4,  color = POPNAME), shape = 17, alpha = .8, size = 4) + 
  scale_colour_manual(values = NA_col_4, name ="Population")  + scale_x_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd North American Populations Principal Components 1 & 2", x = "Principal Component 1 (6.83%)", y = "Principal Component 2 (3.00%)", size = 20))

###Odd NA one and three
ggplot(data = ODD_NA_eigen_table) + geom_point(aes(x = V3, y = V5,  color = POPNAME), shape = 17, alpha = .8, size = 4) + 
  scale_colour_manual(values = NA_col_4, name ="Population") + scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd North American Populations Principal Components 1 & 3", x = "Principal Component 1 (6.83%)", y = "Principal Component 3 (1.53%)", size = 20))

###ODd NA two and three
ggplot(data = ODD_NA_eigen_table) + geom_point(aes(x = V4, y = V5,  color = POPNAME), shape = 17, alpha = .8, size = 4) + 
  scale_colour_manual(values = NA_col_4, name ="Population") +  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd North American Populations Principal Components 2 & 3", x = "Principal Component 2 (3.00%)", y = "Principal Component 3 (1.53%)", size = 20))
dev.off()