### Plotting the Linkage map from LepMap
###    all families together, several different runs
### Carolyn Tarpey | August 2015 
### ---------------------------------------


install.packages("lattice", dependencies = TRUE)
library(lattice)

library(ggplot2)
library(stringr)
library(plyr)



###rwaples/chum_populations/linkage_map/LEPmap/with_paralogs/plotting_LGs.R
#####@rwaples rwaples on Mar 13 consensus map

##plots the counts of the number of markers in a LG

library(ggplot2)

collapsed <- 'G:/Analysis/Mapping/AllHaps/LepMap/ReOrder_LG17_LOD17_MAP.txt'

LINKAGEmap <-"G:\\Analysis\\Mapping\\AllHaps\\LepMap\\LinkageMap.txt"


plot_LG_hist <- function(path_to_LEPmap_chr_file){
  chr_assigments = read.table(path_to_LEPmap_chr_file, header = TRUE)
  assigned_to_LG <- subset(chr_assigments, LG > 0)
  ggplot(data = assigned_to_LG, aes(x = LG)) + geom_histogram(binwidth = 1, color = 'darkgreen', fill = 'grey') + 
    geom_text(label = paste("Unassigned: ", as.character(nrow(chr_assigments)-nrow(assigned_to_LG))), x = 30, y = max(table(assigned_to_LG$V1))-40) + theme_minimal()
} 


plot_LG_hist(collapsed)


########################################


####colors: saffron gold "#FF9933", carolina meets cornflower "#6699FF",  black navy "#000f2d", coral "#ff6666", baby blue "#c6d9ff", 
##slightly bluer bubblegum pink "#ff66ff", pale ballet pink "#ffd9f2", bubble gum "#ff66cc"

dup_color <-c(yes = "#FF9933", no = "#6699FF")
dup_color <-c(yes = "#000f2d", no = "#6699FF")
dup_color <-c(yes = "#ff6666", no = "#6699FF")
dup_color <-c(yes = "#c6d9ff", no = "#ff66ff")
dup_color <-c(yes = "#000f2d", no = "#ff66ff")
dup_color <-c(yes = "#000f2d", no = "#ffd9f2")




map_positions <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\MAP.txt", sep = "\t", header = TRUE)
duplicates <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\MAPdups.txt", sep = "\t", header = TRUE)

map_positions_n13 <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\MAPreLG13.txt", sep = "\t", header = TRUE)
duplicates_n13 <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\MAPreLG13dups.txt", sep = "\t", header = TRUE)
nonduplicates_n13 <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\MAPreLG13NONdups.txt", sep = "\t", header = TRUE)

ggplot(data = map_positions) + geom_point(aes(x = LG, y = position, color = paralog), alpha = .5, size = 4, color = "blue") 
ggplot(data = map_positions) + geom_point(aes(x = LG, y = position, color = paralog), alpha = .5, size = 4) + scale_colour_manual(values = dup_color)
ggplot(data = map_positions) + geom_point(aes(x = LG, y = position, color = paralog), alpha = .5, size = 3) + scale_colour_manual(values = dup_color) + theme_classic()


### show duplicates in the data

dup_color <-c(yes = "#000f2d", no = "maroon1")

ggplot(data = map_positions) + geom_point(aes(x = LG, y = position, color = paralog), alpha = .5, size = 3) + scale_colour_manual(values = dup_color)+ theme_classic()


##plotting only the duplicates
ggplot(data = duplicates) + geom_point(aes(x = LG, y = position), alpha = .5, size = 3, color = "blue") 


##this jitters the duplicates so you can see them 
ggplot(data = duplicates, aes(x = LG, y = position)) + geom_point(position = position_jitter(w = 0.3, h = 0.3), alpha = .5, size = 3, color = "firebrick") 

ggplot(data = duplicates) + geom_jitter(aes(x = LG, y = position), alpha = .5, size = 3, color = "blue", factor = 1.1) 



####################################################################################
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX LOD 16XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
####################################################################################

MAP_lg13_ydups_n24_ndupunsLOD16 <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\MAP_lg13_ydups_n24_ndupunsLOD16.txt", sep = "\t", header = TRUE)
MAP_lg13_ydups_n24_LOD16 <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\MAP_lg13_ydups_n24_LOD16.txt", sep = "\t", header = TRUE)


color <-c(yes = "limegreen", no = "mediumslateblue")
ggplot(data = MAP_lg13_ydups_n24_ndupunsLOD16) + geom_point(aes(x = LG, y = position, color = paralog), alpha = .5, size = 3) + scale_colour_manual(values = color)+ theme_classic()

color <-c(yes = "limegreen", no = "mediumslateblue")
ggplot(data = MAP_lg13_ydups_n24_LOD16) + geom_point(aes(x = LG, y = position, color = paralog), alpha = .5, size = 3) + scale_colour_manual(values = color)+ theme_classic()








####################################################################################
#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX LOD 17 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
####################################################################################


##################         NEW REORDERED 13 XXXXXXXXXXX

###this one works ok to show duplicates in the data
dup_color <-c(yes = "limegreen", no = "mediumslateblue")
ggplot(data = map_positions_n13) + geom_point(aes(x = LG, y = position, color = paralog), alpha = .5, size = 3) + scale_colour_manual(values = dup_color)+ theme_classic()


##plotting only the duplicates
ggplot(data = duplicates_n13) + geom_point(aes(x = LG, y = position), alpha = .5, size = 3, color = "blue") 

##plotting only the nonduplicates
ggplot(data = nonduplicates_n13) + geom_point(aes(x = LG, y = position), alpha = .5, size = 3, color = "blue")


##this jitters the duplicates so you can see them 
ggplot(data = duplicates_n13, aes(x = LG, y = position)) + geom_point(position = position_jitter(w = 0.3, h = 0.3), alpha = .5, size = 3, color = "firebrick") 

ggplot(data = duplicates_n13) + geom_jitter(aes(x = LG, y = position), alpha = .5, size = 3, color = "blue", factor = 1.1) 

##this jitters the non duplicates so you can see them 
ggplot(data = nonduplicates_n13, aes(x = LG, y = position)) + geom_point(position = position_jitter(w = 0.3, h = 0.3), alpha = .2, size = 3, color = "plum2") 

ggplot(data = nonduplicates_n13) + geom_jitter(aes(x = LG, y = position), alpha = .05, size = 3, color = "darkslateblue", factor = 1.1) 

#######     ONLY Markers IN FAM 01, colored on the map of  NEW REORDERED 13 XXXXXXXXXXX

fam01_n13 <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\LG13reorderedFam01.txt", sep = "\t", header = TRUE)

###
fam_color <-c(yes = "limegreen", no = "orchid")
ggplot(data = fam01_n13) + geom_point(aes(x = LG, y = Pos, color = fam_01), alpha = .2, size = 3) + scale_colour_manual(values = fam_color)+ theme_classic()

ggplot(data = fam01_n13) + geom_jitter(aes(x = LG, Pos, color = fam_01), alpha = .5, size = 3, factor = 1.01) + theme_classic()

##################  Forced Reorder of 13 LG13_01as13 XXXXXXXXXXX

map_LG13_01as13 <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\LG13_01as13_MAP.txt", sep = "\t", header = TRUE)

color <-c(yes = "limegreen", no = "mediumslateblue")
ggplot(data = map_LG13_01as13) + geom_point(aes(x = LG, y = position, color = paralog), alpha = .5, size = 3) + scale_colour_manual(values = color)+ theme_classic()


##################  Forced Reorder of 13 LG13_01as0_MAP XXXXXXXXXXX

map_LG13_01as0 <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\LG13_01as0_MAP.txt", sep = "\t", header = TRUE)

color <-c(yes = "limegreen", no = "mediumslateblue")
ggplot(data = map_LG13_01as0) + geom_point(aes(x = LG, y = position, color = paralog), alpha = .5, size = 3) + scale_colour_manual(values = color)+ theme_classic()


################## NONDuplicated loci REORDERED 13 XXXXXXXXXXX

NonDupLG13_MAP <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\NonDupLG13_MAP.txt", sep = "\t", header = TRUE)


color <-c(yes = "limegreen", no = "mediumslateblue")

ggplot(data = NonDupLG13_MAP) + geom_point(aes(x = LG, y = position, color = paralog), alpha = .8, size = 4) + scale_colour_manual(values = color) + theme_classic()

#axis.text.y = "cM", axis.text.x = "Linkage Group Number", + scale_y_discrete(name = "cM") + theme(element_text(size = 14), element_text(size = 14) ) +  scale_x_discrete(breaks=c("1","3","5","7","9","11","13","15","17","19","21","23","25"), labels=c("1","3","5","7","9","11","13","15","17","19","21","23","25"))


p <- ggplot(data = NonDupLG13_MAP) + geom_point(aes(x = LG, y = position, color = paralog), alpha = .8, size = 4) + 
  scale_colour_manual(values = color) + theme_classic()

p + labs(x = "Linkage Group Number", element_text(size = 14), element_text(size = 14,) y = "cM") 


p + theme(axis.line=element_blank(),axis.text.x= "Linkage Group Number",
          axis.text.y= "cM", axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())


plt.scatter(x = NonDupLG13_MAP.LG, y = NonDupLG13_MAP.position, alpha = .7, s = 30, c = 'blue')
plt.xlim(0); plt.ylim(-5); plt.xlabel('Linkage Group Number',fontsize = 24); plt.ylabel("cM",fontsize = 24)
plt.xticks(range(1, max(NonDupLG13_MAP.LG)+1), fontsize = 16)
plt.show()


#################################################
# LG24_LG13_Ydups_N24_Lod16_MAP
LG24_LG13_Ydups_N24_Lod16_MAP <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\LG24_LG13_Ydups_N24_Lod16_MAP.txt", sep = "\t", header = TRUE)


color <-c(yes = "limegreen", no = "mediumslateblue")

ggplot(data = LG24_LG13_Ydups_N24_Lod16_MAP) + geom_point(aes(x = LG, y = position, color = paralog), alpha = .8, size = 4) + scale_colour_manual(values = color) + theme_classic()



#################################################
# LG24_LG13_manforce_24Nodup_LD17_MAP
LG24_LG13_manforce_24Nodup_LD17_MAP <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\LG24_LG13_manforce_24Nodup_LD17_MAP.txt", sep = "\t", header = TRUE)

color <-c(yes = "limegreen", no = "mediumslateblue")

ggplot(data = LG24_LG13_manforce_24Nodup_LD17_MAP) + geom_point(aes(x = LG, y = position, color = paralog), alpha = .8, size = 4) + scale_colour_manual(values = color) + theme_classic()



#################################################
# LG24_LG13_ManForceCut_N24D_LD17_MAP.txt
LG24_LG13_ManForceCut_N24D_LD17_MAP <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\LG24_LG13_ManForceCut_N24D_LD17_MAP.txt", sep = "\t", header = TRUE)

color <-c(yes = "#8EDB33", no = "darkslateblue")
color <-c(yes = "#57982F", no = "darkslateblue")

ggplot(data = LG24_LG13_ManForceCut_N24D_LD17_MAP) + geom_point(aes(x = LG, y = position, color = paralog), alpha = .8, size = 6) + scale_colour_manual(values = color) + theme_classic()


 #################################################
# LOD17_LG24.map
LOD17_LG24_MAP <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\LOD17_LG24_MAP.txt", sep = "\t", header = TRUE)

color <-c(yes = "limegreen", no = "mediumslateblue")

ggplot(data = LOD17_LG24_MAP) + geom_point(aes(x = LG, y = position, color = paralog), alpha = .8, size = 4) + scale_colour_manual(values = color) + theme_classic()

#################################################
# LG13_LG24_110H_65_MAP.txt
LG13_LG24_110H_65_MAP <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\LG13_LG24_110H_65_MAP.txt", sep = "\t", header = TRUE)

color <-c(yes = "limegreen", no = "mediumslateblue")

ggplot(data = LG13_LG24_110H_65_MAP) + geom_point(aes(x = LG, y = position, color = paralog), alpha = .8, size = 4) + scale_colour_manual(values = color) + theme_classic()

#################################################
# LG13_LG24_01H_MAP.txt
LG13_LG24_01H_MAP <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\LG13_LG24_01H_MAP.txt", sep = "\t", header = TRUE)

color <-c(yes = "limegreen", no = "#708090")
color <-c(yes = "limegreen", no = "darkslateblue")

ggplot(data = LG13_LG24_01H_MAP) +geom_line(aes(x= LG, y= position)) + geom_point(aes(x = LG, y = position, color = paralog), alpha = .6, size = 4) + scale_colour_manual(values = color) + theme_classic()
######


#################################################
# ReOrder_LG17_LOD17_MAP.txt this is made with the # LG24_LG13_ManForceCut_N24D_LD17_MAP.txt file for the rest of the loci.

ReOrder_LG17_LOD17_MAP <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\ReOrder_LG17_LOD17_MAP.txt", sep = "\t", header = TRUE)

color <-c(yes = "limegreen", no = "#708090")
color <-c(yes = "#57982F", no = "darkslateblue")

ggplot(data = ReOrder_LG17_LOD17_MAP) + geom_line(aes(x= LG, y= position))  + geom_point(aes(x = LG, y = position, color = paralog), alpha = .6, size = 4) + scale_colour_manual(values = color) + theme_classic()


#################################################
#LinkageMap.txt  this is made with ReOrder_LG17_LOD17_MAP.txt and the # LG24_LG13_ManForceCut_N24D_LD17_MAP.txt file for the rest of the loci.

LinkageMap <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\LinkageMap.txt", sep = "\t", header = TRUE)

color <-c(yes = "limegreen", no = "#708090")
color <-c(yes = "#57982F", no = "darkslateblue")
ASP_color <- c(yes = "orange", no ="blue" )

names(LinkageMap)  


ggplot(data = LinkageMap) +  geom_point(aes(x = LepLG, y = Position, color = Paralog), alpha = .6, size = 4) + scale_colour_manual(values = color) + theme_classic()
  
ggplot(data = LinkageMap) + geom_point(aes(x = LepLG, y = Position), color = 'blue', alpha = .5) +
  geom_segment(aes(x = LepLG, y = Position, xend =  LepLG  , yend = Position), color = 'black', alpha = .2) +
  theme_classic()

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

LinkageMap <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\LinkageMap.txt", sep = "\t", header = TRUE)
Lines <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\Lines_behind.txt", sep = "\t", header = TRUE)
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




