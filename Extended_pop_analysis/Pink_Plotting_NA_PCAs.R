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
pop_key <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/POPINFO_LS_susitna.txt", header = TRUE, sep = '\t')
head(pop_key)

############Eigen tables #############
###Import the eigenvalues that were run in PLINK and merge them with the Pop info

All_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/all_pops_pca.eigenvec")
All_eigen_table <- merge(All_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
All_eigen_table <- All_eigen_table[order(All_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(All_eigen_table)

All_even_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/even_pops_pca.eigenvec")
All_even_eigen_table <- merge(All_even_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
All_even_eigen_table <- All_even_eigen_table[order(All_even_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(All_even_eigen_table)

All_odd_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/odd_pops_pca.eigenvec")
All_odd_eigen_table <- merge(All_odd_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
All_odd_eigen_table <- All_odd_eigen_table[order(All_odd_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(All_odd_eigen_table)

##North America
NA_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/na_all_pca.eigenvec")
NA_eigen_table <- merge(NA_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
NA_eigen_table <- NA_eigen_table[order(NA_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(NA_eigen_table)

EVEN_NA_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/even_na_pca.eigenvec")
EVEN_NA_eigen_table <- merge(EVEN_NA_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
EVEN_NA_eigen_table <- EVEN_NA_eigen_table[order(EVEN_NA_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(EVEN_NA_eigen_table)

ODD_NA_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/odd_na_pca.eigenvec")
ODD_NA_eigen_table <- merge(ODD_NA_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
ODD_NA_eigen_table <- ODD_NA_eigen_table[order(ODD_NA_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(ODD_NA_eigen_table)

NA_nome_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/na_nome_pops_pca.eigenvec")
NA_nome_eigen_table <- merge(NA_nome_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
NA_nome_eigen_table <- NA_nome_eigen_table[order(NA_nome_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(NA_nome_eigen_table)

EVEN_NA_nome_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/e_na_nome_pca.eigenvec")
EVEN_NA_nome_eigen_table <- merge(EVEN_NA_nome_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
EVEN_NA_nome_eigen_table <- EVEN_NA_nome_eigen_table[order(EVEN_NA_nome_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(EVEN_NA_nome_eigen_table)

ODD_NA_nome_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/o_na_nome_pca.eigenvec")
ODD_NA_nome_eigen_table <- merge(ODD_NA_nome_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
ODD_NA_nome_eigen_table <- ODD_NA_nome_eigen_table[order(ODD_NA_nome_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(ODD_NA_nome_eigen_table)

NA_nome_NS_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/na_nome_NS_pops_pca.eigenvec")
NA_nome_NS_eigen_table <- merge(NA_nome_NS_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
NA_nome_NS_eigen_table <- NA_nome_NS_eigen_table[order(NA_nome_NS_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(NA_nome_NS_eigen_table)

## Asia 
Asia_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/a_all_pca.eigenvec")
Asia_eigen_table <- merge(Asia_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
Asia_eigen_table <- Asia_eigen_table[order(Asia_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(Asia_eigen_table)

EVEN_A_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/even_a_pca.eigenvec")
EVEN_A_eigen_table <- merge(EVEN_A_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
EVEN_A_eigen_table <- EVEN_A_eigen_table[order(EVEN_A_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(EVEN_A_eigen_table)

ODD_A_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/odd_a_pca.eigenvec")
ODD_A_eigen_table <- merge(ODD_A_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
ODD_A_eigen_table <- ODD_A_eigen_table[order(ODD_A_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(ODD_A_eigen_table)


Asia_Susitna_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/Asia_nome_susitna_pops_pca.eigenvec")
Asia_Susitna_eigen_table <- merge(Asia_Susitna_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
Asia_Susitna_eigen_table <- Asia_Susitna_eigen_table[order(Asia_Susitna_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(Asia_Susitna_eigen_table)

#############COLORS################################################
#plot(rep(1,7),col=ALL_col,pch=19,cex=4)

OG_col<- c("Hokkaido"="#61005e","Amur"="#cd5490","Magadan"="#e0957e","Kamchatka"="#c9cf4a","Norton Sound"="#677f3e", "Prince William Sound"="#7297ee","Puget Sound"="#00158a")
All_col<- c("Hokkaido"="#61005e","Amur"="#cd5490","Magadan"="#e0957e","Kamchatka"="#c9cf4a","Norton Sound"="#677f3e","Susitna" = "#599091", 
            "Prince William Sound"="#7297ee","Skeena" = "#2171b5", "Puget Sound"="#00158a")

NA_col_4<-  c("Susitna" = "#599091", "Prince William Sound"="#7297ee","Skeena" = "#2171b5", "Puget Sound"="#00158a")
NA_col_5<-  c("Norton Sound"="#677f3e","Susitna" = "#599091", "Prince William Sound"="#7297ee","Skeena" = "#2171b5", "Puget Sound"="#00158a")
NA_col_5_NS<-  c("Norton Sound"="#677f3e", "Prince William Sound"="#7297ee","Skeena" = "#2171b5", "Puget Sound"="#00158a")
Asia_col <- c("Norton Sound"="#677f3e","Kamchatka"="#c9cf4a","Magadan"="#e0957e","Amur"="#cd5490","Hokkaido"="#61005e")

##THese have susitna in a grey
OG_col<- c("Hokkaido"="#61005e","Amur"="#cd5490","Magadan"="#e0957e","Kamchatka"="#c9cf4a","Norton Sound"="#677f3e", "Prince William Sound"="#7297ee","Puget Sound"="#00158a")
All_col<- c("Hokkaido"="#61005e","Amur"="#cd5490","Magadan"="#e0957e","Kamchatka"="#c9cf4a","Norton Sound"="#677f3e","Susitna" = "#8FA9B7", 
            "Prince William Sound"="#7297ee","Skeena" = "#2171b5", "Puget Sound"="#00158a")

NA_col_4<-  c("Susitna" = "#8FA9B7", "Prince William Sound"="#7297ee","Skeena" = "#2171b5", "Puget Sound"="#00158a")
NA_col_5<-  c("Norton Sound"="#677f3e","Susitna" = "#8FA9B7", "Prince William Sound"="#7297ee","Skeena" = "#2171b5", "Puget Sound"="#00158a")
NA_col_5_NS<-  c("Norton Sound"="#677f3e", "Prince William Sound"="#7297ee","Skeena" = "#2171b5", "Puget Sound"="#00158a")
Asia_col <- c("Norton Sound"="#677f3e","Kamchatka"="#c9cf4a","Magadan"="#e0957e","Amur"="#cd5490","Hokkaido"="#61005e")
Asia_sus_col <- c("Hokkaido"="#61005e","Amur"="#cd5490","Magadan"="#e0957e","Kamchatka"="#c9cf4a","Norton Sound"="#677f3e","Susitna" = "#8FA9B7")


#################  All pops, ALL LOCI ###############################
###ALL one and two
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/ALL_pcas_grey.pdf", width = 9, height = 7)

ggplot(data = All_eigen_table) + geom_point(aes(x =V3, y = V4,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) + scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") + 
  labs(list(title= "All Populations PC #s 1 & 2", x = "Principal Component 1 (21.30%) ", y = "Principal Component 2 (9.97%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Susitna", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) + guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))

ggplot(data = All_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line())+ scale_x_continuous(trans="reverse") +
  labs(list(title= "All Populations PC #s 1 & 3", x = "Principal Component 1 (21.30%)", y = "Principal Component 3 (9.06%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Susitna", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) + guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))

ggplot(data = All_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) +  labs(list(title= "All Populations PC #s 2 & 3", 
  x = "Principal Component 2 (9.97%)", y = "Principal Component 3 (9.06%)", size = 20))+  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan",
  "Kamchatka","Norton Sound","Susitna", "Prince William Sound","Skeena", "Puget Sound"), values = All_col) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))
dev.off()

### ALL Even one and two
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/All_even_pcas.pdf", width = 9, height = 7)

ggplot(data = All_even_eigen_table) + geom_point(aes(x =V3, y = V4,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) +  scale_y_continuous(trans="reverse") + 
  labs(list(title= "All Even Populations PC #s 1 & 2", x = "Principal Component 1 (5.93%)", y = "Principal Component 2 (4.44%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Susitna", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) +  guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))

ggplot(data = All_even_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line())+scale_y_continuous(trans="reverse") +
  labs(list(title= "All Even Populations PC #s 1 & 3", x = "Principal Component 1 (5.93%)", y = "Principal Component 3 (2.51%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Susitna", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) +  guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))

ggplot(data = All_even_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) + scale_y_continuous(trans="reverse") + 
  labs(list(title= "All Even Populations PC #s 2 & 3", x = "Principal Component 2 (4.44%)", y = "Principal Component 3 (2.51%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Susitna", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) +  guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))
dev.off()

###ALL Odd 
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/All_odd_pcas.pdf", width = 9, height = 7)

ggplot(data = All_odd_eigen_table) + geom_point(aes(x =V3, y = V4,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) + scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") + 
  labs(list(title= "All Odd Populations PC #s 1 & 2", x = "Principal Component 1 (8.73%)", y = "Principal Component 2 (7.21%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Susitna", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) +  guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))

ggplot(data = All_odd_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line())+ scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") + 
  labs(list(title= "All Odd Populations PC #s 1 & 3", x = "Principal Component 1 (8.73%)", y = "Principal Component 3 (3.78%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Susitna", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) +  guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))

ggplot(data = All_odd_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) + scale_y_continuous(trans="reverse") +
  labs(list(title= "All Odd Populations PC #s 2 & 3", x = "Principal Component 2 (7.21%)", y = "Principal Component 3 (3.78%) ", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Susitna", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) +  guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))
dev.off()


################### NA NO NOME #############################


pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/NA_NoNome_pcas.pdf", width = 9, height = 7)
###ALL one and two
ggplot(data = NA_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All North American Populations PC #s 1 & 2", x = "Principal Component 1 (11.51%)", y = "Principal Component 2 (7.38%)", size = 20)) + 
  scale_color_manual(breaks = c("Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))

###ALL one and three
ggplot(data = NA_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) + scale_x_continuous(trans="reverse") +  
  labs(list(title= "All North American Populations PC #s 1 & 3", x = "Principal Component 1 (11.51%)", y = "Principal Component 3 (4.48%)", size = 20))+ 
  scale_color_manual(breaks = c("Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))

###ALL two and three
ggplot(data = NA_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) + 
  labs(list(title= "All North American Populations PC #s 2 & 3", x = "Principal Component 2 (7.38%)", y = "Principal Component 3 (4.48%)", size = 20))+ 
  scale_color_manual(breaks = c("Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))
dev.off()


pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/EVEN_NA_pcas.pdf", width = 9, height = 7)
###Even NA one and two
ggplot(data = EVEN_NA_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
  legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Even North American Populations PC #s 1 & 2", x = "Principal Component 1 (4.00%)", y = "Principal Component 2 (1.94%)", size = 20))+ 
  scale_color_manual(breaks = c("Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))

###Even NA one and three
ggplot(data = EVEN_NA_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
  legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Even North American Populations PC #s 1 & 3", x = "Principal Component 1 (4.00%)", y = "Principal Component 3 (1.29%)", size = 20))+ 
  scale_color_manual(breaks = c("Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))

###Even NA two and three
ggplot(data = EVEN_NA_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) +
  scale_y_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
  legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) + scale_x_continuous(trans="reverse")+
  labs(list(title= "Even North American Populations PC #s 2 & 3", x = "Principal Component 2 (1.94%%)", y = "Principal Component 3 (1.29%)", size = 20))+ 
  scale_color_manual(breaks = c("Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))
dev.off()

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/ODD_NA_pcas.pdf", width = 9, height = 7)
###odd NA one and two
ggplot(data = ODD_NA_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
  legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd North American Populations PC #s 1 & 2", x = "Principal Component 1 (6.83%)", y = "Principal Component 2 (3.00%)", size = 20))+ 
  scale_color_manual(breaks = c("Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))

###Odd NA one and three
ggplot(data = ODD_NA_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd North American Populations PC #s 1 & 3", x = "Principal Component 1 (6.83%)", y = "Principal Component 3 (1.54%)", size = 20))+ 
  scale_color_manual(breaks = c("Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))

###ODd NA two and three
ggplot(data = ODD_NA_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14),axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd North American Populations PC #s 2 & 3", x = "Principal Component 2 (3.00%)", y = "Principal Component 3 (1.54%)", size = 20))+ 
  scale_color_manual(breaks = c("Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))
dev.off()

 
################### NA NOME #############################

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/NA_nome_pcas_grey.pdf", width = 9, height = 7)
###All NA Nome one and two
ggplot(data = NA_nome_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION, shape = LINEAGE),  alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
  legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) + scale_y_continuous(trans="reverse")+
  labs(list(title= "North American Populations PC #s 1 & 2", x = "Principal Component 1 (12.73%)", y = "Principal Component 2 (7.56%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5)))

###All NA Nome one and three
ggplot(data = NA_nome_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "North American Populations PC #s 1 & 3", x = "Principal Component 1 (12.73%) ", y = "Principal Component 3 (4.64%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5)))

###All NA Nome two and three
ggplot(data = NA_nome_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "North American Populations PC #s 2 & 3", x = "Principal Component 2 (7.56%) ", y = "Principal Component 3 (4.64%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) +
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5)))
dev.off()


pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/EVEN_NA_nome_pcas.pdf", width = 9, height = 7)
###Even NA Nome one and two
ggplot(data = EVEN_NA_nome_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") +  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
  legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Even North American Populations PC #s 1 & 2", x = "Principal Component 1 (4.15%)", y = "Principal Component 2 (2.27%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5))) + scale_y_continuous(trans="reverse") 

###Even NA Nome one and three
ggplot(data = EVEN_NA_nome_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +   theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Even North American Populations PC #s 1 & 3", x = "Principal Component 1 (4.15%) ", y = "Principal Component 3 (1.78%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5)))

###Even NA Nome two and three
ggplot(data = EVEN_NA_nome_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  scale_y_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Even North American Populations PC #s 2 & 3", x = "Principal Component 2 (2.27%) ", y = "Principal Component 3 (1.78%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5)))
dev.off()

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/ODD_NA_nome_pcas.pdf", width = 9, height = 7)
###odd NA Nome one and two
ggplot(data = ODD_NA_nome_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION ), shape = 17, alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") +  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
  legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd North American Populations PC #s 1 & 2", x = "Principal Component 1 (7.05%)", y = "Principal Component 2 (3.65%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5))) + scale_y_continuous(trans="reverse") 

###Odd NA Nome one and three
ggplot(data = ODD_NA_nome_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") +  theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14),axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd North American Populations PC #s 1 & 3", x = "Principal Component 1 (7.05%) ", y = "Principal Component 3 (2.90%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5)))

###ODd NA Nome two and three
ggplot(data = ODD_NA_nome_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
   theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd North American Populations PC #s 2 & 3", x = "Principal Component 2 (3.65%) ", y = "Principal Component 3 (2.90%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Susitna", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5)))
dev.off()

########### NA Nome, NO Susitna #################

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/NA_nome_NS_pcas.pdf", width = 9, height = 7)

###odd NA Nome one and two
ggplot(data = NA_nome_NS_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") +  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
  legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "North American Pops, +Nome -Susitna PC#s 1 & 2", x = "Principal Component 1 (10.57%)", y = "Principal Component 2 (4.10%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5_NS) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5_NS)))  

###Odd NA Nome one and three
ggplot(data = NA_nome_NS_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") +  theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14),axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "North American Pops, +Nome -Susitna PC#s 1 & 3", x = "Principal Component 1 (10.57%) ", y = "Principal Component 3 (3.55%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5_NS) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5_NS)))

###ODd NA Nome two and three
ggplot(data = NA_nome_NS_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "North American Pops, +Nome -Susitna PC#s 2 & 3", x = "Principal Component 2 (4.10%) ", y = "Principal Component 3 (3.55%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5_NS) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5_NS)))

dev.off()

########### Asia  #################
Asia_col

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/Even_Asia_pcas.pdf", width = 9, height = 7)

###Asian one and two
ggplot(data = Asia_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") +  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
  legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Asian Populations PC#s 1 & 2", x = "Principal Component 1 (13.04%)", y = "Principal Component 2 (3.57%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Kamchatka", "Magadan", "Amur", "Hokkaido"), values = Asia_col) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=Asia_col)))  

###Asian one and three
ggplot(data = Asia_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") +  theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14),axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Asian Populations PC#s 1 & 3", x = "Principal Component 1 (13.04%) ", y = "Principal Component 3 (2.66%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Kamchatka", "Magadan", "Amur", "Hokkaido"), values = Asia_col) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=Asia_col))) + scale_y_continuous(trans="reverse") 

###Asian two and three
ggplot(data = Asia_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
   theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Asian Populations PC#s 2 & 3", x = "Principal Component 2 (3.57%) ", y = "Principal Component 3 (2.66%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Kamchatka", "Magadan", "Amur", "Hokkaido"), values = Asia_col) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=Asia_col))) + scale_y_continuous(trans="reverse") 

dev.off()


pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/Even_Asia_pcas.pdf", width = 9, height = 7)

###Even Asian one and two
ggplot(data = EVEN_A_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") +  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
                                                                 legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Even Lineage Asian Populations PC#s 1 & 2", x = "Principal Component 1 (2.19%)", y = "Principal Component 2 (1.33%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Kamchatka", "Magadan", "Amur", "Hokkaido"), values = Asia_col) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=Asia_col)))  

###Even Asian one and three
ggplot(data = EVEN_A_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") +  theme_classic() + theme(text = element_text(size= 17), 
                                                                 legend.title = element_blank(),legend.text = element_text(size= 14),axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Even Lineage Asian Populations PC#s 1 & 3", x = "Principal Component 1 (2.19%) ", y = "Principal Component 3 (1.07%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Kamchatka", "Magadan", "Amur", "Hokkaido"), values = Asia_col) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=Asia_col)))

###Even Asian two and three
ggplot(data = EVEN_A_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), 
                                                                legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Even Lineage Asian Populations PC#s 2 & 3", x = "Principal Component 2 (1.33%) ", y = "Principal Component 3 (1.07%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Kamchatka", "Magadan", "Amur", "Hokkaido"), values = Asia_col) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=Asia_col)))

dev.off()


pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/Odd_Asia_pcas.pdf", width = 9, height = 7)

###odd Asian one and two
ggplot(data = ODD_A_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") +  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
                                                                 legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd Lineage Asian Populations PC#s 1 & 2", x = "Principal Component 1 (2.55%)", y = "Principal Component 2 (1.75%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Kamchatka", "Magadan", "Amur", "Hokkaido"), values = Asia_col) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=Asia_col)))  

###Odd Asian one and three
ggplot(data = ODD_A_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") +  theme_classic() + theme(text = element_text(size= 17), 
                                                                 legend.title = element_blank(),legend.text = element_text(size= 14),axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd Lineage Asian Populations PC#s 1 & 3", x = "Principal Component 1 (2.55%) ", y = "Principal Component 3 (1.35%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Kamchatka", "Magadan", "Amur", "Hokkaido"), values = Asia_col) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=Asia_col)))

###Odd Asian two and three
ggplot(data = ODD_A_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd Lineage Asian Populations PC#s 2 & 3", x = "Principal Component 2 (1.75%) ", y = "Principal Component 3 (1.35%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Kamchatka", "Magadan", "Amur", "Hokkaido"), values = Asia_col) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=Asia_col)))

dev.off()

##########Asia with Susitna################


pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/Asia_Susitna_pcas.pdf", width = 9, height = 7)

###odd Asian one and two
ggplot(data = Asia_Susitna_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") +  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
  legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Asian Populations With Susitna PC#s 1 & 2", x = "Principal Component 1 (14.88%)", y = "Principal Component 2 (8.34%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Susitna"), values = Asia_sus_col) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=Asia_sus_col)))  

###Odd Asian one and three
ggplot(data = Asia_Susitna_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") +  theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14),axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Asian Populations With Susitna PC#s 1 & 3", x = "Principal Component 1 (14.88%) ", y = "Principal Component 3 (4.80%)", size = 20))+ 
  scale_color_manual(breaks = c("Susitna","Norton Sound","Kamchatka", "Magadan", "Amur", "Hokkaido"), values = Asia_sus_col) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=Asia_sus_col)))

###Odd Asian two and three
ggplot(data = Asia_Susitna_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Asian Populations With Susitna PC#s 2 & 3", x = "Principal Component 2 (8.34%) ", y = "Principal Component 3 (4.80%)", size = 20))+ 
  scale_color_manual(breaks = c("Susitna","Norton Sound","Kamchatka", "Magadan", "Amur", "Hokkaido"), values = Asia_sus_col) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=Asia_sus_col)))

dev.off()
