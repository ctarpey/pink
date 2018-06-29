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
pop_key <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/POPINFO_LS.txt", header = TRUE, sep = '\t')
head(pop_key)

############Eigen tables #############
###Import the eigenvalues that were run in PLINK and merge them with the Pop info
#OG pops for QC
OG_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/og_pops_pca.eigenvec")
OG_eigen_table <- merge(OG_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
OG_eigen_table <- OG_eigen_table[order(OG_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(OG_eigen_table)

OG_16662_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/og_16662_pops_pca.eigenvec")
OG_16662_eigen_table <- merge(OG_16662_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
OG_16662_eigen_table <- OG_16662_eigen_table[order(OG_16662_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(OG_16662_eigen_table)

OG_newloci_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/og_newloci_pops_pca.eigenvec")
OG_newloci_eigen_table <- merge(OG_newloci_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
OG_newloci_eigen_table <- OG_newloci_eigen_table[order(OG_newloci_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(OG_newloci_eigen_table)

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

All_31485_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs31485loci/All_31485_pops_pca.eigenvec")
All_31485_eigen_table <- merge(All_31485_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
All_31485_eigen_table <- All_31485_eigen_table[order(All_31485_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(All_31485_eigen_table)

NA_nome_NS_eigen_table <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/na_nome_NS_pops_pca.eigenvec")
NA_nome_NS_eigen_table <- merge(NA_nome_NS_eigen_table, pop_key, by.x = 'V1', by.y = 'POP')
NA_nome_NS_eigen_table <- NA_nome_NS_eigen_table[order(NA_nome_NS_eigen_table$Order_geo),] #sort by a geographical order, then odd then even 
head(NA_nome_NS_eigen_table)


#############COLORS################################################
#plot(rep(1,7),col=ALL_col,pch=19,cex=4)

OG_col<- c("Hokkaido"="#61005e","Amur"="#cd5490","Magadan"="#e0957e","Kamchatka"="#c9cf4a","Norton Sound"="#677f3e", "Prince William Sound"="#7297ee","Puget Sound"="#00158a")
All_col<- c("Hokkaido"="#61005e","Amur"="#cd5490","Magadan"="#e0957e","Kamchatka"="#c9cf4a","Norton Sound"="#677f3e","Cook Inlet" = "#563e7f", 
            "Prince William Sound"="#7297ee","Skeena" = "#2171b5", "Puget Sound"="#00158a")

NA_col_4<-  c("Cook Inlet" = "#563e7f", "Prince William Sound"="#7297ee","Skeena" = "#2171b5", "Puget Sound"="#00158a")
NA_col_5<-  c("Norton Sound"="#677f3e","Cook Inlet" = "#563e7f", "Prince William Sound"="#7297ee","Skeena" = "#2171b5", "Puget Sound"="#00158a")
NA_col_5_NS<-  c("Norton Sound"="#677f3e", "Prince William Sound"="#7297ee","Skeena" = "#2171b5", "Puget Sound"="#00158a")
Asia_col <- c("Hokkaido"="#61005e","Amur"="#cd5490","Magadan"="#e0957e","Kamchatka"="#c9cf4a","Norton Sound"="#677f3e")

#################  All pops, 31485 LOCI ###############################
###ALL one and two
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/ALL_31485_pcas.pdf", width = 9, height = 7)

ggplot(data = All_31485_eigen_table) + geom_point(aes(x =V3, y = V4,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) + scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") + 
  labs(list(title= "All Populations PC #s 1 & 2", x = "Principal Component 1 ", y = "Principal Component 2 ", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Cook Inlet", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) + guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))

ggplot(data = All_31485_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line())+ scale_x_continuous(trans="reverse") +
  labs(list(title= "All Populations PC #s 1 & 3", x = "Principal Component 1 ", y = "Principal Component 3", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Cook Inlet", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) + guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))

ggplot(data = All_31485_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) +  labs(list(title= "All Populations PC #s 2 & 3", 
  x = "Principal Component 2", y = "Principal Component 3", size = 20))+  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan",
  "Kamchatka","Norton Sound","Cook Inlet", "Prince William Sound","Skeena", "Puget Sound"), values = All_col) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))

dev.off()
#################  All pops, ALL LOCI ###############################
###ALL one and two
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/ALL_pcas.pdf", width = 9, height = 7)

ggplot(data = All_eigen_table) + geom_point(aes(x =V3, y = V4,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) + scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") + 
  labs(list(title= "All Populations PC #s 1 & 2", x = "Principal Component 1 (21.30%) ", y = "Principal Component 2 (9.97%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Cook Inlet", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) + guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))

ggplot(data = All_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line())+ scale_x_continuous(trans="reverse") +
  labs(list(title= "All Populations PC #s 1 & 3", x = "Principal Component 1 (21.30%)", y = "Principal Component 3 (9.06%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Cook Inlet", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) + guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))

ggplot(data = All_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) +  labs(list(title= "All Populations PC #s 2 & 3", 
  x = "Principal Component 2 (9.97%)", y = "Principal Component 3 (9.06%)", size = 20))+  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan",
  "Kamchatka","Norton Sound","Cook Inlet", "Prince William Sound","Skeena", "Puget Sound"), values = All_col) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))
dev.off()

### ALL Even one and two
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/All_even_pcas.pdf", width = 9, height = 7)

ggplot(data = All_even_eigen_table) + geom_point(aes(x =V3, y = V4,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) +  scale_y_continuous(trans="reverse") + 
  labs(list(title= "All Even Populations PC #s 1 & 2", x = "Principal Component 1 (5.93%)", y = "Principal Component 2 (4.44%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Cook Inlet", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) +  guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))

ggplot(data = All_even_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line())+scale_y_continuous(trans="reverse") +
  labs(list(title= "All Even Populations PC #s 1 & 3", x = "Principal Component 1 (5.93%)", y = "Principal Component 3 (2.51%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Cook Inlet", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) +  guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))

ggplot(data = All_even_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) + scale_y_continuous(trans="reverse") + 
  labs(list(title= "All Even Populations PC #s 2 & 3", x = "Principal Component 2 (4.44%)", y = "Principal Component 3 (2.51%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Cook Inlet", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) +  guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))
dev.off()

###ALL Odd 
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/All_odd_pcas.pdf", width = 9, height = 7)

ggplot(data = All_odd_eigen_table) + geom_point(aes(x =V3, y = V4,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) + scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") + 
  labs(list(title= "All Odd Populations PC #s 1 & 2", x = "Principal Component 1 (8.73%)", y = "Principal Component 2 (7.21%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Cook Inlet", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) +  guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))

ggplot(data = All_odd_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line())+ scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") + 
  labs(list(title= "All Odd Populations PC #s 1 & 3", x = "Principal Component 1 (8.73%)", y = "Principal Component 3 (3.78%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Cook Inlet", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) +  guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))

ggplot(data = All_odd_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) + scale_y_continuous(trans="reverse") +
  labs(list(title= "All Odd Populations PC #s 2 & 3", x = "Principal Component 2 (7.21%)", y = "Principal Component 3 (3.78%) ", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound","Cook Inlet", "Prince William Sound",
  "Skeena", "Puget Sound"), values = All_col) +  guides(color = guide_legend(override.aes = list(shape = 22, fill=All_col)))
dev.off()

################# OG All LOCI ###############################
###OG ALL one and two
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/OG_pcas.pdf", width = 9, height = 7)

ggplot(data = OG_eigen_table) + geom_point(aes(x =V3, y = V4,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) + scale_x_continuous(trans="reverse") + 
  labs(list(title= "All OG Populations PC #s 1 & 2", x = "Principal Component 1 (17.24%)", y = "Principal Component 2 (8.62%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound", "Prince William Sound","Puget Sound"), values = OG_col)+
  guides(color = guide_legend(override.aes = list(shape = 22, fill=OG_col)))

ggplot(data = OG_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
                          axis.line.x = element_line(), axis.line.y = element_line())+ scale_x_continuous(trans="reverse") +
  labs(list(title= "All OG Populations PC #s 1 & 3", x = "Principal Component 1 (17.24%)", y = "Principal Component 3 (5.12%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound", "Prince William Sound","Puget Sound"), values = OG_col)+
  guides(color = guide_legend(override.aes = list(shape = 22, fill=OG_col)))

ggplot(data = OG_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All OG Populations PC #s 2 & 3", x = "Principal Component 2 (8.62%)", y = "Principal Component 3 (5.12%)", 4size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound", "Prince William Sound","Puget Sound"), values = OG_col)+
  guides(color = guide_legend(override.aes = list(shape = 22, fill=OG_col)))

dev.off()



#################OG 16662 (13106 IRL)###############################
###OG 16662 (13106 IRL) ALL one and two
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/OG_16662_pcas.pdf", width = 9, height = 7)

ggplot(data = OG_16662_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) + scale_y_continuous(trans="reverse") + scale_x_continuous(trans="reverse") +
  labs(list(title= "All OG 16662 loci Populations PC #s 1 & 2", x = "Principal Component 1 (18.56%)", y = "Principal Component 2 (8.37%)", size = 20)) + 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound", "Prince William Sound","Puget Sound"), values = OG_col)+
  guides(color = guide_legend(override.aes = list(shape = 22, fill=OG_col)))

ggplot(data = OG_16662_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line())+ scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  labs(list(title= "All OG Populations 16662 loci  PC #s 1 & 3", x = "Principal Component 1 (18.56%)", y = "Principal Component 3 (5.00%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound", "Prince William Sound","Puget Sound"), values = OG_col)+
  guides(color = guide_legend(override.aes = list(shape = 22, fill=OG_col)))

ggplot(data = OG_16662_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) +  scale_y_continuous(trans="reverse")+ scale_x_continuous(trans="reverse")+
  labs(list(title= "All OG Populations 16662 loci  PC #s 2 & 3", x = "Principal Component 2 (8.37%)", y = "Principal Component 3 (5.00%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound", "Prince William Sound","Puget Sound"), values = OG_col)+
  guides(color = guide_legend(override.aes = list(shape = 22, fill=OG_col)))
dev.off()


#################### OG POPS NEW LOCI 10653############################
###OG New LOci only , 10653  ALL one and two
pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/OG_newloci_pcas.pdf", width = 9, height = 7)

ggplot(data = OG_newloci_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All OG Only New loci Populations PC #s 1 & 2", x = "Principal Component 1 (15.82%)", y = "Principal Component 2 (10.45%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound", "Prince William Sound","Puget Sound"), values = ALL_col)+
  guides(color = guide_legend(override.aes = list(shape = 22, fill=ALL_col)))

ggplot(data = OG_newloci_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line())+ scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  labs(list(title= "All OG Populations Only New loci  PC #s 1 & 3", x = "Principal Component 1 (15.82%)", y = "Principal Component 3 (5.18%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound", "Prince William Sound","Puget Sound"), values = ALL_col)+
  guides(color = guide_legend(override.aes = list(shape = 22, fill=ALL_col)))

ggplot(data = OG_newloci_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All OG Populations Only New loci  PC #s 2 & 3", x = "Principal Component 2 (10.45%) ", y = "Principal Component 3 (5.18%)", size = 20))+ 
  scale_color_manual(breaks = c("Hokkaido","Amur","Magadan","Kamchatka","Norton Sound", "Prince William Sound","Puget Sound"), values = ALL_col)+
  guides(color = guide_legend(override.aes = list(shape = 22, fill=ALL_col)))

dev.off()

################### NA NO NOME #############################


pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/NA_NoNome_pcas.pdf", width = 9, height = 7)
###ALL one and two
ggplot(data = NA_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All North American Populations PC #s 1 & 2", x = "Principal Component 1 (11.51%)", y = "Principal Component 2 (7.38%)", size = 20)) + 
  scale_color_manual(breaks = c("Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))

###ALL one and three
ggplot(data = NA_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) + scale_x_continuous(trans="reverse") +  
  labs(list(title= "All North American Populations PC #s 1 & 3", x = "Principal Component 1 (11.51%)", y = "Principal Component 3 (4.48%)", size = 20))+ 
  scale_color_manual(breaks = c("Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))

###ALL two and three
ggplot(data = NA_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),legend.text = element_text(size= 14),
  axis.line.x = element_line(), axis.line.y = element_line()) + 
  labs(list(title= "All North American Populations PC #s 2 & 3", x = "Principal Component 2 (7.38%)", y = "Principal Component 3 (4.48%)", size = 20))+ 
  scale_color_manual(breaks = c("Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))
dev.off()


pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/EVEN_NA_pcas.pdf", width = 9, height = 7)
###Even NA one and two
ggplot(data = EVEN_NA_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
  legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Even North American Populations PC #s 1 & 2", x = "Principal Component 1 (4.00%)", y = "Principal Component 2 (1.94%)", size = 20))+ 
  scale_color_manual(breaks = c("Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))

###Even NA one and three
ggplot(data = EVEN_NA_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
  legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Even North American Populations PC #s 1 & 3", x = "Principal Component 1 (4.00%)", y = "Principal Component 3 (1.29%)", size = 20))+ 
  scale_color_manual(breaks = c("Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))

###Even NA two and three
ggplot(data = EVEN_NA_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) +
  scale_y_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
  legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) + scale_x_continuous(trans="reverse")+
  labs(list(title= "Even North American Populations PC #s 2 & 3", x = "Principal Component 2 (1.94%%)", y = "Principal Component 3 (1.29%)", size = 20))+ 
  scale_color_manual(breaks = c("Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))
dev.off()

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/ODD_NA_pcas.pdf", width = 9, height = 7)
###odd NA one and two
ggplot(data = ODD_NA_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
  legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd North American Populations PC #s 1 & 2", x = "Principal Component 1 (6.83%)", y = "Principal Component 2 (3.00%)", size = 20))+ 
  scale_color_manual(breaks = c("Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))

###Odd NA one and three
ggplot(data = ODD_NA_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd North American Populations PC #s 1 & 3", x = "Principal Component 1 (6.83%)", y = "Principal Component 3 (1.54%)", size = 20))+ 
  scale_color_manual(breaks = c("Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))

###ODd NA two and three
ggplot(data = ODD_NA_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14),axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd North American Populations PC #s 2 & 3", x = "Principal Component 2 (3.00%)", y = "Principal Component 3 (1.54%)", size = 20))+ 
  scale_color_manual(breaks = c("Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_4) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_4)))
dev.off()

 
################### NA NOME #############################

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/NA_nome_pcas.pdf", width = 9, height = 7)
###All NA Nome one and two
ggplot(data = NA_nome_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION, shape = LINEAGE),  alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
  legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) + scale_y_continuous(trans="reverse")+
  labs(list(title= "North American Populations PC #s 1 & 2", x = "Principal Component 1 (12.73%)", y = "Principal Component 2 (7.56%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5)))

###All NA Nome one and three
ggplot(data = NA_nome_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "North American Populations PC #s 1 & 3", x = "Principal Component 1 (12.73%) ", y = "Principal Component 3 (4.64%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5)))

###All NA Nome two and three
ggplot(data = NA_nome_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "North American Populations PC #s 2 & 3", x = "Principal Component 2 (7.56%) ", y = "Principal Component 3 (4.64%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) +
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5)))
dev.off()


pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/EVEN_NA_nome_pcas.pdf", width = 9, height = 7)
###Even NA Nome one and two
ggplot(data = EVEN_NA_nome_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") +  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
  legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Even North American Populations PC #s 1 & 2", x = "Principal Component 1 (4.15%)", y = "Principal Component 2 (2.27%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5))) + scale_y_continuous(trans="reverse") 

###Even NA Nome one and three
ggplot(data = EVEN_NA_nome_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +   theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Even North American Populations PC #s 1 & 3", x = "Principal Component 1 (4.15%) ", y = "Principal Component 3 (1.78%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5)))

###Even NA Nome two and three
ggplot(data = EVEN_NA_nome_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  scale_y_continuous(trans="reverse") + theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Even North American Populations PC #s 2 & 3", x = "Principal Component 2 (2.27%) ", y = "Principal Component 3 (1.78%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5)))
dev.off()

pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/OneTagSNP/PLINK/PCAs/ODD_NA_nome_pcas.pdf", width = 9, height = 7)
###odd NA Nome one and two
ggplot(data = ODD_NA_nome_eigen_table) + geom_point(aes(x = V3, y = V4,  color = LOCATION ), shape = 17, alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") +  theme_classic() + theme(text = element_text(size= 17), legend.title = element_blank(),
  legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd North American Populations PC #s 1 & 2", x = "Principal Component 1 (7.05%)", y = "Principal Component 2 (3.65%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5))) + scale_y_continuous(trans="reverse") 

###Odd NA Nome one and three
ggplot(data = ODD_NA_nome_eigen_table) + geom_point(aes(x = V3, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
  scale_x_continuous(trans="reverse") +  theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14),axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd North American Populations PC #s 1 & 3", x = "Principal Component 1 (7.05%) ", y = "Principal Component 3 (2.90%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) + 
  guides(color = guide_legend(override.aes = list(shape = 22, fill=NA_col_5)))

###ODd NA Nome two and three
ggplot(data = ODD_NA_nome_eigen_table) + geom_point(aes(x = V4, y = V5,  color = LOCATION), shape = 17, alpha = .8, size = 4) + 
   theme_classic() + theme(text = element_text(size= 17), 
  legend.title = element_blank(),legend.text = element_text(size= 14), axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd North American Populations PC #s 2 & 3", x = "Principal Component 2 (3.65%) ", y = "Principal Component 3 (2.90%)", size = 20))+ 
  scale_color_manual(breaks = c("Norton Sound","Cook Inlet", "Prince William Sound", "Skeena", "Puget Sound"), values = NA_col_5) + 
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



# outputFile <- file("C:/Users/Carolyn/Desktop/PLINK_PED.txt", "wb")
# write.table(PLINK_PED,outputFile,quote=FALSE,row.names=FALSE,col.names=FALSE,eol="\n")
# close(outputFile)
# 
# PLINK_PED
