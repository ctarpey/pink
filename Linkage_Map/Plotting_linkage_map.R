### Plotting the Linkage map from LepMap
###    all families together, several different runs
### Carolyn Tarpey | August 2015 
### ---------------------------------------

install.packages("plyr")
install.packages("ggplot2")


library(lattice)
library(ggplot2)
library(stringr)
library(plyr)


LinkageMap2 <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\LinkageMap.txt", sep = "\t", header = TRUE)


LinkageMap <- read.table("G:\\Analysis\\Mapping\\AllHaps\\LepMap\\ComparingLepNo01to01MSTMap.txt", sep = "\t", header = TRUE)


names(LinkageMap)  

other_color <-c(yes = "limegreen", no = "#708090")
color <-c(yes = "#57982F", no = "darkslateblue")
pink_black_color <-c(yes = "#000f2d", no = "maroon1")
lighter_color <-c(yes = "limegreen", no = "mediumslateblue")



ggplot(data = LinkageMap) + geom_point(aes(x = LepLG, y = LepMapPosition, color = Paralog), alpha = .5, size = 3) + scale_colour_manual(values = color) + theme_classic()

ggplot(data = LinkageMap) + geom_point(aes(x = LepLgNew13_24, y = newLepPos, color = Paralog), alpha = .5, size = 3) + scale_colour_manual(values = color) + theme_classic()


p <- ggplot(data = LinkageMap) + geom_point(aes(x = LepLG, y = Position, color = Paralog), alpha = .5, size = 4) + scale_colour_manual(values = color) + theme_classic()


p + xlab("Linkage Group Number") + xlab("cM")


p + labs(x = "Linkage Group Number", y = "cM") + theme(text = element_text(size = 18))
 

# +geom_segment(aes(x = LepLG, y = Position,  yend = c(0, max(Position)), size = .1))

p + theme(axis.line=element_blank(),axis.title.x =element_text("Linkage Group Number") ,
          axis.title.y = element_text("cM") , axis.ticks=element_blank(),legend.position="top",
          panel.background=element_blank(), panel.border=element_blank(), panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(), plot.background=element_blank())


plt.scatter(x = NonDupLG13_MAP.LG, y = NonDupLG13_MAP.position, alpha = .7, s = 30, c = 'blue')
plt.xlim(1); plt.ylim(-1); plt.xlabel('Linkage Group Number',fontsize = 24); plt.ylabel("cM",fontsize = 24)
plt.xticks(range(1, max(NonDupLG13_MAP.LG)+1), fontsize = 16)
plt.show()


ggplot(data = LinkageMap) +  geom_point(aes(x = LepLG, y = Position, color = Paralog), alpha = .6, size = 4) + scale_colour_manual(values = color) + theme_classic()

ggplot(data = LinkageMap) + geom_point(aes(x = LepLG, y = Position), color = 'blue', alpha = .5) +
  geom_segment(aes(x = LepLG, y = Position, xend =  LepLG  , yend = Position), color = 'black', alpha = .2) +
  theme_classic()


