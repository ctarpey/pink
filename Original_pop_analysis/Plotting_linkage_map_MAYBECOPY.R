### Plotting the Linkage map from LepMap
###    plots for talks, with lines and centromeres
### Carolyn Tarpey | January 2016; updated October 2106 (at end of code)
### ---------------------------------------

library(ggplot2)
library(stringr)
library(plyr)


################################################# PLOT FOR ASP TALK NOV 2015

#LinkageMap.txt  this is made with ReOrder_LG17_LOD17_MAP.txt and the # LG24_LG13_ManForceCut_N24D_LD17_MAP.txt file for the rest of the loci.

LinkageMap <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\LinkageMap.txt", sep = "\t", header = TRUE)
Lines <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\Lines_behind.txt", sep = "\t", header = TRUE)

ASP_color <- c(yes = "#FD3C06", no ="blue")

names(LinkageMap) 
names(Lines)

ggplot(data = LinkageMap) + geom_point(aes(x = LepLG, y = Position, color = Paralog), alpha = .5, size = 4) + 
  scale_colour_manual(values = ASP_color) + theme_classic() + scale_x_continuous(breaks=1:26) + ylim(0, 150)

ggplot(data = LinkageMap) + geom_segment(data = Lines, (aes(x = x, y = y, xend = xend, yend = yend, alpha = .3))) +
  geom_point(aes(x = LepLG, y = Position , color = Paralog), alpha = .5, size = 4) + scale_colour_manual(values = ASP_color) + 
  theme_classic() + scale_x_continuous(breaks=c(1,3,5,7,9,11,13,15,17,19,21,23,26))

ggplot(data = LinkageMap) + geom_segment(data = Lines, mapping= aes(x = x, y = y, xend = xend, yend = yend), size =.5, color = "black") +
  geom_point(aes(x = LepLG, y = Position , color = Paralog), alpha = .4, size = 3.5) + scale_colour_manual(values = ASP_color) +
  scale_x_continuous(breaks=c(2,4,6,8,10,12,14,16,18,20,22,24,26)) + theme_classic() + xlab("Linkage Group") + ylab("Genetic Distance (cM)") +
  theme(axis.title.x = element_text(color = "black", size = 20), axis.text.x = element_text(vjust = 0.5, size= 14)) +
  theme(axis.title.y = element_text(color = "black", size = 20), axis.text.y = element_text(vjust = 0.5, size= 14)) +
  theme(legend.title = element_text(color = "black", size = 16), legend.text = element_text(vjust = 0.5, size= 14))

ggsave(filename = "G:\\Analysis\\Mapping\\AllHaps\\LepMap\\ASP_plots\\ASP_pink_cmt_8x6_sm.pdf", width = 8, height= 6)
ggsave(filename = "G:\\Analysis\\Mapping\\AllHaps\\LepMap\\ASP_plots\\ASP_pink_cmt_8x8_sm.pdf", width = 8, height= 8)


################################################# Corrected LINKAGE MAP December 2015

#LinkageMap.txt this is made with reordered LG 17, LOD16  and the # LG24_LG13_ManForceCut_N24D_LD17_MAP.txt file for the rest of the loci.

LinkageMap <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\LinkageMap.txt", sep = "\t", header = TRUE)
Lines <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\Lines_behind.txt", sep = "\t", header = TRUE)

ASP_color <- c(yes = "#FD3C06", no ="blue")

names(LinkageMap) 
names(Lines)

ggplot(data = LinkageMap) + geom_point(aes(x = LepLG, y = Position, color = Paralog), alpha = .5, size = 4) + 
  scale_colour_manual(values = ASP_color) + theme_classic() + scale_x_continuous(breaks=1:26) + ylim(0, 150)

ggplot(data = LinkageMap) + geom_segment(data = Lines, (aes(x = x, y = y, xend = xend, yend = yend, alpha = .3))) +
  geom_point(aes(x = LepLG, y = Position , color = Paralog), alpha = .5, size = 4) + scale_colour_manual(values = ASP_color) + 
  theme_classic() + scale_x_continuous(breaks=c(1,3,5,7,9,11,13,15,17,19,21,23,26))

ggplot(data = LinkageMap) + geom_segment(data = Lines, mapping= aes(x = x, y = y, xend = xend, yend = yend), size =.5, color = "black") +
  geom_point(aes(x = LepLG, y = Position , color = Paralog), alpha = .4, size = 3.5) + scale_colour_manual(values = ASP_color) +
  scale_x_continuous(breaks=c(2,4,6,8,10,12,14,16,18,20,22,24,26)) + theme_classic() + xlab("Linkage Group") + ylab("Genetic Distance (cM)") +
  theme(axis.title.x = element_text(color = "black", size = 20), axis.text.x = element_text(vjust = 0.5, size= 14)) +
  theme(axis.title.y = element_text(color = "black", size = 20), axis.text.y = element_text(vjust = 0.5, size= 14)) +
  theme(legend.title = element_text(color = "black", size = 16), legend.text = element_text(vjust = 0.5, size= 14))

ggsave(filename = "G:\\Analysis\\Mapping\\AllHaps\\LepMap\\Updated_LG17_8x6_sm.pdf", width = 8, height= 6)
ggsave(filename = "G:\\Analysis\\Mapping\\AllHaps\\LepMap\\Updated_LG17_8x8_sm.pdf", width = 8, height= 8)


#################################################  LINKAGE MAP With Centromeres, January 2016

#LinkageMap.txt this is made with reordered LG 17, LOD16  and the # LG24_LG13_ManForceCut_N24D_LD17_MAP.txt file for the rest of the loci.

LinkageMap <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LinkageMap.txt", sep = "\t", header = TRUE)
Lines <- read.table("G:\\Analysis\\Mapping\\AllHaps\\Lines_behind.txt", sep = "\t", header = TRUE)
Centromeres <- read.table("G:\\Analysis\\Mapping\\AllHaps\\Centromeres_4plots.txt", sep = "\t", header = TRUE)

ASP_color <- c(yes = "#FD3C06", no ="blue")

names(LinkageMap) 
names(Lines)
names(Centromeres)

# ggplot(data = LinkageMap) + geom_point(aes(x = NewLG, y = LepMap_Position, color = Paralog), alpha = .5, size = 4) + 
#   scale_colour_manual(values = ASP_color) + theme_classic() + scale_x_continuous(breaks=1:26) + ylim(0, 150)
# 
# ggplot(data = LinkageMap) + geom_segment(data = Lines, (aes(x = x, y = y, xend = xend, yend = yend, alpha = .3))) +
#   geom_point(aes(x = NewLG, y = LepMap_Position , color = Paralog), alpha = .5, size = 4) + scale_colour_manual(values = ASP_color) + 
#   theme_classic() + scale_x_continuous(breaks=c(1,3,5,7,9,11,13,15,17,19,21,23,26))

ggplot(data = LinkageMap) + geom_segment(data = Lines, mapping= aes(x = x, y = y, xend = xend, yend = yend), size =.5, color = "black") +
  geom_segment(data = Centromeres, mapping= aes(x = x, y = y, xend = xend, yend = yend), size =4.9, color = "limegreen") +
  geom_point(aes(x = NewLG, y = LepMap_Position , color = Paralog), alpha = .4, size = 3.5) + scale_colour_manual(values = ASP_color) +
  scale_x_continuous(breaks=c(2,4,6,8,10,12,14,16,18,20,22,24,26)) + theme_classic() + xlab("Linkage Group") + ylab("Genetic Distance (cM)") +
  theme(axis.title.x = element_text(color = "black", size = 20), axis.text.x = element_text(vjust = 0.5, size= 14)) +
  theme(axis.title.y = element_text(color = "black", size = 20), axis.text.y = element_text(vjust = 0.5, size= 14)) +
  theme(legend.title = element_text(color = "black", size = 16), legend.text = element_text(vjust = 0.5, size= 14))

ggsave(filename = "G:\\Analysis\\Mapping\\AllHaps\\Linkage_cents_8x6_sm.pdf", width = 8, height= 6)
ggsave(filename = "G:\\Analysis\\Mapping\\AllHaps\\Linkage_cents_8x8_sm.pdf", width = 8, height= 8)



################################################# Updated plot with new LG from GARRETT

#LINKAGEmap_GU.txt

LinkageMap <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LINKAGEmap_GU.txt", sep = "\t", header = TRUE)
Lines <- read.table("G:\\Analysis\\Mapping\\AllHaps\\Lines_behind_GU.txt", sep = "\t", header = TRUE)
Centromeres <- read.table("G:\\Analysis\\Mapping\\AllHaps\\Centromeres_4plots_GU.txt", sep = "\t", header = TRUE)

ASP_color <- c(yes = "#FD3C06", no ="blue")

names(LinkageMap) 
names(Lines)
names(Centromeres)


ggplot(data = LinkageMap) + geom_segment(data = Lines, mapping= aes(x = x, y = y, xend = xend, yend = yend), size =.5, color = "black") +
  geom_segment(data = Centromeres, mapping= aes(x = x, y = y, xend = xend, yend = yend), size =4.9, color = "limegreen") +
  geom_point(aes(x = NewLG, y = LepMap_Position , color = Paralog), alpha = .4, size = 3.5) + scale_colour_manual(values = ASP_color) +
  scale_x_continuous(breaks=c(2,4,6,8,10,12,14,16,18,20,22,24,26)) + theme_classic() + xlab("Linkage Group") + ylab("Genetic Distance (cM)") +
  theme(axis.title.x = element_text(color = "black", size = 20), axis.text.x = element_text(vjust = 0.5, size= 14)) +
  theme(axis.title.y = element_text(color = "black", size = 20), axis.text.y = element_text(vjust = 0.5, size= 14)) +
  theme(legend.title = element_text(color = "black", size = 16), legend.text = element_text(vjust = 0.5, size= 14))


ggsave(filename = "G:\\Analysis\\Mapping\\AllHaps\\LepMap\\PINK_garrett_updated_8x6_sm.pdf", width = 8, height= 6)
ggsave(filename = "G:\\Analysis\\Mapping\\AllHaps\\LepMap\\PINK_garrett_updated__8x8_sm.pdf", width = 8, height= 8)


################################################# Updated plot with new LG from GARRETT_ without Centromeres for Report

#LINKAGEmap_GU.txt

LinkageMap <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LINKAGEmap_GU.txt", sep = "\t", header = TRUE)
Lines <- read.table("G:\\Analysis\\Mapping\\AllHaps\\Lines_behind_GU.txt", sep = "\t", header = TRUE)
Centromeres <- read.table("G:\\Analysis\\Mapping\\AllHaps\\Centromeres_4plots_GU.txt", sep = "\t", header = TRUE)

ASP_color <- c(yes = "#FD3C06", no ="blue")

names(LinkageMap) 
names(Lines)
names(Centromeres)


ggplot(data = LinkageMap) + geom_segment(data = Lines, mapping= aes(x = x, y = y, xend = xend, yend = yend), size =.5, color = "black") +
  #geom_segment(data = Centromeres, mapping= aes(x = x, y = y, xend = xend, yend = yend), size =4.9, color = "limegreen") +
  geom_point(aes(x = NewLG, y = LepMap_Position , color = Paralog), alpha = .4, size = 3.5) + scale_colour_manual(values = ASP_color) +
  scale_x_continuous(breaks=c(2,4,6,8,10,12,14,16,18,20,22,24,26)) + theme_classic() + xlab("Linkage Group") + ylab("Genetic Distance (cM)") +
  theme(axis.title.x = element_text(color = "black", size = 20), axis.text.x = element_text(vjust = 0.5, size= 14)) +
  theme(axis.title.y = element_text(color = "black", size = 20), axis.text.y = element_text(vjust = 0.5, size= 14)) +
  theme(legend.title = element_text(color = "black", size = 16), legend.text = element_text(vjust = 0.5, size= 14))


ggsave(filename = "G:\\Analysis\\Mapping\\AllHaps\\LepMap\\PINK_garrett_updated_nocents_8x6_sm.pdf", width = 8, height= 6)
ggsave(filename = "G:\\Analysis\\Mapping\\AllHaps\\LepMap\\PINK_garrett_updated_nocents_8x8_sm.pdf", width = 8, height= 8)
