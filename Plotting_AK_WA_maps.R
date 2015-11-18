### Plotting the WA and AK linkage maps from LepMap
###    Alaska and WA families on seperate maps
### Carolyn Tarpey | September 2015 
### ---------------------------------------

install.packages("plyr")
install.packages("ggplot2")

library(lattice)
library(ggplot2)
library(stringr)
library(plyr)


AK_map <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\Alaska\\AK_map.txt", sep = "\t", header = TRUE)

WA_map <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\Washington\\WA_map.txt", sep = "\t", header = TRUE)

names(WA_map)  

other_color <-c(yes = "limegreen", no = "#708090")
color <-c(yes = "#57982F", no = "darkslateblue")
lighter_color <-c(yes = "limegreen", no = "mediumslateblue")


################################## Alaska


ggplot(data = AK_map) + geom_point(aes(x = LG, y = position, color = duplicate), alpha = .5, size = 3) + scale_colour_manual(values = color) + theme_classic()


p <- ggplot(data = AK_map) + geom_point(aes(x = LG, y = position, color = duplicate), alpha = .5, size = 3) + scale_colour_manual(values = color) + theme_classic()


p + xlab("Linkage Group Number") + xlab("cM")


p + labs(title= "Alaska odd-year pink linkage map", x = "Linkage Group Number", y = "cM") + theme(text = element_text(size = 18))


# +geom_segment(aes(x = LG, y = position,  yend = c(0, max(position)), size = .1))

#p + theme(axis.line=element_blank(),axis.title.x =element_text("Linkage Group Number") ,
#          axis.title.y = element_text("cM") , axis.ticks=element_blank(),legend.position="top",
#          panel.background=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(),
#          panel.grid.minor=element_blank(), plot.background=element_blank())


#plt.scatter(x = LG, y = position, alpha = .7, s = 30, c = 'blue')
#plt.xlim(1); plt.ylim(-1); plt.xlabel('Linkage Group Number',fontsize = 24); plt.ylabel("cM",fontsize = 24)
#plt.xticks(range(1, max(LG)+1), fontsize = 16)
#plt.show()


ggplot(data = AK_map) +  geom_point(aes(x = LG, y = position, color = duplicate), alpha = .6, size = 4) + scale_colour_manual(values = color) + theme_classic()

ggplot(data = AK_map) + geom_point(aes(x = LG, y = position), color = 'blue', alpha = .5) +
  geom_segment(aes(x = LG, y = position, xend =  LG  , yend = position), color = 'black', alpha = .2) +
  theme_classic()


################################## Washington 



ggplot(data = WA_map) + geom_point(aes(x = LG, y = position, color = duplicate), alpha = .5, size = 3) + scale_colour_manual(values = color) + theme_classic()



p <- ggplot(data = WA_map) + geom_point(aes(x = LG, y = position, color = duplicate), alpha = .5, size = 3) + scale_colour_manual(values = color) + theme_classic()


p + xlab("Linkage Group Number") + xlab("cM") + labs(title= "Washington odd-year pink linkage map", x = "Linkage Group Number", y = "cM") + theme(text = element_text(size = 18))


# +geom_segment(aes(x = LG, y = position,  yend = c(0, max(position)), size = .1))

#p + theme(axis.line=element_blank(),axis.title.x =element_text("Linkage Group Number") ,
#          axis.title.y = element_text("cM") , axis.ticks=element_blank(),legend.position="top",
#          panel.background=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(),
#          panel.grid.minor=element_blank(), plot.background=element_blank())


#plt.scatter(x = LG, y = position, alpha = .7, s = 30, c = 'blue')
#plt.xlim(1); plt.ylim(-1); plt.xlabel('Linkage Group Number',fontsize = 24); plt.ylabel("cM",fontsize = 24)
#plt.xticks(range(1, max(LG)+1), fontsize = 16)
#plt.show()


ggplot(data = WA_map) +  geom_point(aes(x = LG, y = position, color = duplicate), alpha = .6, size = 4) + scale_colour_manual(values = color) + theme_classic()

ggplot(data = WA_map) + geom_point(aes(x = LG, y = position), color = 'blue', alpha = .5) +
  geom_segment(aes(x = LG, y = position, xend =  LG  , yend = position), color = 'black', alpha = .2) +
  theme_classic()



