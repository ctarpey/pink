### Isolation by Distance calculations for populations
###    Distances in Kilometers, from ArcGIS FST standardized
### Carolyn Tarpey | June 2016
### ---------------------------------------

library(ape)
library(ggplot2)

setwd('G:/Analysis/Pop_analysis/Populations_b3_may/IBD')

#################################################################
#########  Kerry Run         ###################################
#################################################################

#making matrices of 16681 Fst 1-fst
all16881 <- read.table ("ALL_1-fst.txt", sep = "\t")
all16881 <- as.matrix(all16881)  

ALLO <- read.table ("ALLo_1-fst.txt", sep = "\t")
ALLO <- as.matrix(ALLO)  

ALLE <- read.table ("ALLe_1-fst.txt", sep = "\t")
ALLE <- as.matrix(ALLE) 

AsiaO <- read.table ("ASIAO_1-fst.txt", sep = "\t")
AsiaO <- as.matrix(AsiaO)  

AsiaE <- read.table ("ASIAE_1-fst.txt", sep = "\t")
AsiaE <- as.matrix(AsiaE)  

USAO <- read.table ("USAO_1-fst.txt", sep = "\t")
USAO <- as.matrix(USAO)  

USAE <- read.table ("USAE_1-fst.txt", sep = "\t")
USAE <- as.matrix(USAE) 



####coastline
allgeo_coast <- read.table ("ALL_coastal.txt", sep = "\t")
allgeo_coast  <- as.matrix(allgeo_coast)  

allsmall_coast <- read.table ("ALLo_coastal.txt", sep = "\t")
allsmall_coast  <- as.matrix(allsmall_coast) 

ASIAgeo_coast  <- read.table ("ASIA_coastal.txt", sep = "\t")
ASIAgeo_coast  <- as.matrix(ASIAgeo_coast)  

USAgeo_coast  <- read.table ("USA_coastal.txt", sep = "\t")
USAgeo_coast  <- as.matrix(USAgeo_coast)  

  
####mantel tests of the fst matrices against the geographic distance

ASIAoddMantel <- mantel.test(AsiaO, ASIAgeo_coast, nperm = 1000 )
ASIAevenMantel <- mantel.test(AsiaE, ASIAgeo_coast, nperm = 1000)
USAoddMantel <- mantel.test(USAO, USAgeo_coast, nperm = 1000)
USAevenMantel <- mantel.test(USAE, USAgeo_coast, nperm = 1000)
ALLMantel <- mantel.test(all16881, allgeo_coast, nperm = 1000)

ALLoMantel <- mantel.test(ALLO, allsmall_coast, nperm = 1000 )
ALLeMantel <- mantel.test(ALLE, allsmall_coast, nperm = 1000)


# Results
ASIAoddMantel 
ASIAevenMantel
USAoddMantel
USAevenMantel
ALLMantel
ALLoMantel
ALLeMantel

#################################################################
######### plotting the relationship ###################################
#################################################################

#plotting IBD 

#cols_OddEven <- c(O = "navy", L = "darkgrey", E = "magenta") 
#cols_group <- c(Beringia = "darkorange", America = "purple") 

melted_fst_dist <- read.table ("Plotting_IBD_ALL_edit.txt", header = TRUE, sep = "\t")
melted_fst_dist_origial <- read.table ("Plotting_IBD_ALL.txt", header = TRUE, sep = "\t")
head(melted_fst_dist)



#subset to dif version LFMM_e_LA<-Fst_mapped[Fst_mapped$even_LA=="yes",]
USAe <- melted_fst_dist[melted_fst_dist$USAgroup=="USAe",]
USAo <- melted_fst_dist[melted_fst_dist$USAgroup=="USAo",]
Asiae <- melted_fst_dist[melted_fst_dist$AsiaGroup=="AsiaE",]
Asiao <- melted_fst_dist[melted_fst_dist$AsiaGroup=="AsiaO",]
All <- melted_fst_dist[melted_fst_dist$ALL=="ALL",]

All_e <- melted_fst_dist[melted_fst_dist$Lineage=="E",]
All_o <- melted_fst_dist[melted_fst_dist$Lineage=="O",]

ggplot(data = melted_fst_dist, aes(x= coast, y=FST)) +
  geom_line(data= USAe, aes(x= coast, y=FST, color= "green" )) +
  geom_line(data= USAo, aes(x= coast, y=FST, color= "orange")) +
  geom_line(data= Asiao, aes(x= coast, y=FST, color= "blue")) +
  geom_line(data= Asiae, aes(x= coast, y=FST, color= "red")) +
  geom_point(size = 4, shape = 21)


ggplot(data = melted_fst_dist, aes(x= coast, y=FST)) +
  geom_point(data = melted_fst_dist, size = 4, shape = 21, color = "black") +
  geom_smooth(data= melted_fst_dist, method = lm, se= FALSE, color = "black") +
  geom_smooth(data= USAe, method = lm, se= FALSE, color = "blue") +
  geom_point(data= USAe, size = 4, shape = 21, color = "blue") +
  geom_smooth(data= USAo, method = lm, se= FALSE, color = "lightblue") +
  geom_point(data= USAo, size = 4, shape = 21, color = "lightblue") +
  geom_smooth(data= Asiao, method = lm, se= FALSE, color = "magenta") +
  geom_point(data= Asiao, size = 4, shape = 21, color = "magenta") +
  geom_smooth(data= Asiae, method = lm, se= FALSE, color = "orange") +
  geom_point(data= Asiae, size = 4, shape = 21, color = "orange") 
  

ggplot(data = All, aes(x= coast, y=FST)) +
  geom_smooth(data= All_e, method = lm, se= FALSE, color = "#00CC66") +
  geom_point(data= All_e, size = 4, shape = 21, color = "#00CC66") +
  geom_smooth(data= All_o, method = lm, se= FALSE, color = "purple") +
  geom_point(data= All_o, size = 4, shape = 21, color = "purple") +
  geom_smooth(data= USAe, method = lm, se= FALSE, color = "blue") +
  geom_point(data= USAe, size = 4, shape = 21, color = "blue") +
  geom_smooth(data= USAo, method = lm, se= FALSE, color = "#00CCCC") +
  geom_point(data= USAo, size = 4, shape = 21, color = "#00CCCC") +
  geom_smooth(data= Asiao, method = lm, se= FALSE, color = "magenta") +
  geom_point(data= Asiao, size = 4, shape = 21, color = "magenta") +
  geom_smooth(data= Asiae, method = lm, se= FALSE, color = "orange") +
  geom_point(data= Asiae, size = 4, shape = 21, color = "orange") +
  theme_classic() +
  geom_text(size = 6, color = "#00CC66", aes(175000, 0.06, label = "All even")) +
  geom_text(size = 6, color = "purple", aes(175000, 0.085, label = "All odd")) +
  geom_text(size = 6, color = "blue", aes(25000, 0.05, label = "USA even")) +
  geom_text(size = 6, color = "#00CCCC", aes(25000, 0.07, label = "USA odd")) +
  geom_text(size = 6, color = "magenta", aes(60000, 0.015, label = "Asia even")) +
  geom_text(size = 6, color = "orange", aes(60000, 0.025, label = "Asia odd")) +
  theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
        axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( x = "Distance between populations along the coast (km)", y = "Fst/(1-Fst)", size = 30))

###################################################################


#equation for linear regression 
# http://stackoverflow.com/questions/7549694/ggplot2-adding-regression-line-equation-and-r2-on-graph

##################################################################

lm_eqn <- function(df){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));                 
}

p1 <- p + geom_text(x = 25, y = 300, label = lm_eqn(df), parse = TRUE)

#broken up into one plot for each group
#plot dimensions set by all odd, the largest group


p6<-ggplot(data = All, aes(x= coast, y=FST)) +
  geom_smooth(data= All_e, method = lm, se= FALSE, color = "#00CC66") +
  geom_point(data= All_e, size = 4, shape = 21, color = "#00CC66") +
  theme_classic()+
  scale_x_continuous(limits = c(0, 280000)) +  scale_y_continuous(limits = c(0, 0.12)) +
  ggtitle("All even") + theme(plot.title = element_text(size = 15, color = "#00CC66" )) +
  theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
        axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( x = "", y = "")) + geom_text(x = 200000, y = 0.125, label = lm_eqn(All_e), parse = TRUE)
  

p3<-ggplot(data = All, aes(x= coast, y=FST)) +  
  geom_smooth(data= All_o, method = lm, se= FALSE, color = "purple") +
  geom_point(data= All_o, size = 4, shape = 24, color = "purple") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 280000)) +  scale_y_continuous(limits = c(0, 0.12)) +
  ggtitle("All odd") + theme(plot.title = element_text(size = 15, color = "purple" )) +
  theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
        axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( x = "", y = "Fst/(1-Fst)", size = 15))

p4<-ggplot(data = All, aes(x= coast, y=FST)) +  
  geom_smooth(data= USAe, method = lm, se= FALSE, color = "blue") +
  geom_point(data= USAe, size = 4, shape = 21, color = "blue") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 280000)) +  scale_y_continuous(limits = c(0, 0.12)) +
  ggtitle("USA even") + theme(plot.title = element_text(size = 15, color = "blue" )) +
  theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
        axis.line.x = element_line(), axis.line.y = element_line())  +
  labs(list( x = "", y = ""))

p1<-ggplot(data = All, aes(x= coast, y=FST)) +
  geom_smooth(data= USAo, method = lm, se= FALSE, color = "#00CCCC") +
  geom_point(data= USAo, size = 4, shape = 24, color = "#00CCCC") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 280000)) +  scale_y_continuous(limits = c(0, 0.12)) +
  ggtitle("USA odd") + theme(plot.title = element_text(size = 15, color = "#00CCCC" )) +
  theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
        axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( x = "", y = "Fst/(1-Fst)", size = 15))
  
p2<-ggplot(data = All, aes(x= coast, y=FST)) +
  geom_smooth(data= Asiao, method = lm, se= FALSE, color = "magenta") +
  geom_point(data= Asiao, size = 4, shape = 24, color = "magenta") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 280000)) +  scale_y_continuous(limits = c(0, 0.12)) +
  ggtitle("Asia odd") + theme(plot.title = element_text(size = 15, color = "magenta" )) +
  theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
        axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( x = "", y = "Fst/(1-Fst)", size = 15))

  
p5<-ggplot(data = All, aes(x= coast, y=FST)) +  
  geom_smooth(data= Asiae, method = lm, se= FALSE, color = "orange") +
  geom_point(data= Asiae, size = 4, shape = 21, color = "orange") +
  theme_classic() +
  scale_x_continuous(limits = c(0, 280000)) +  scale_y_continuous(limits = c(0, 0.12)) +
  ggtitle("Asia even") + theme(plot.title = element_text(size = 15, color = "orange" )) +
  theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
        axis.line.x = element_line(), axis.line.y = element_line())  +
  labs(list( x = "", y = ""))


############################################################################################  
############################################################################################
#     Multiple plots on one page
#   from: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
############################################################################################  
############################################################################################  

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


multiplot( p1, p2, p3, p4, p5, p6, cols=2)


###################################################################
#CODE TO PLOT LINNEAR FITS 
#https://susanejohnston.wordpress.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/

#linear fits
ALL_line <- lm(FST ~ coast, data = All)
USA_odd_line <- lm(FST ~ coast, data = USAo)
USA_even_line <- lm(FST ~ coast, data = USAe)
Aisa_odd_line <- lm(FST ~ coast, data = Asiao)
Asia_even_line <- lm(FST ~ coast, data = Asiae)

ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       " Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}

ggplotRegression(ALL_line)
ggplotRegression(USA_odd_line)
ggplotRegression(USA_even_line)
ggplotRegression(Aisa_odd_line)
ggplotRegression(Asia_even_line)

###################################################################

#################################################################
#########      Run for thesis ###################################
#################################################################

#making matrices of 16681 Fst
all16881 <- read.table ("ALL_1-fst_ordered_MATRIX.txt", sep = "\t")
all16881 <- as.matrix(all16881)  

odd16881 <- read.table ("Odd_all_1-fst_ordered_MATRIX.txt", sep = "\t")
odd16881 <- as.matrix(odd16881)  
even16881 <- read.table ("even_all_1-fst_ordered_MATRIX.txt", sep = "\t")
even16881 <- as.matrix(even16881)  

oddBER <- read.table ("odd_Beringia_1-fst_ordered_MATRIX.txt", sep = "\t")
oddBER <- as.matrix(oddBER)  
evenBER <- read.table ("even_Beringia_1-fst_ordered_MATRIX.txt", sep = "\t")
evenBER <- as.matrix(evenBER) 


#making geography matrices
####crosses the ocean
allgeo_ocean <- read.table ("geography_MATRIX_ocean.txt", sep = "\t")
allgeo_ocean  <- as.matrix(allgeo_ocean)  
oddgeo_ocean  <- read.table ("geography_MATRIX_odd_ocean.txt", sep = "\t")
oddgeo_ocean  <- as.matrix(oddgeo_ocean)  
evengeo_ocean  <- read.table ("geography_MATRIX_odd_ocean.txt", sep = "\t")
evengeo_ocean  <- as.matrix(evengeo_ocean)  
BER_ocean  <- read.table ("geography_MATRIX_BER_ocean.txt", sep = "\t")
BER_ocean  <- as.matrix(BER_ocean)


####coastline
allgeo_coast <- read.table ("geography_MATRIX_coast.txt", sep = "\t")
allgeo_coast  <- as.matrix(allgeo_coast)  
oddgeo_coast  <- read.table ("geography_MATRIX_odd_coast.txt", sep = "\t")
oddgeo_coast  <- as.matrix(oddgeo_coast)  
evengeo_coast  <- read.table ("geography_MATRIX_odd_coast.txt", sep = "\t")
evengeo_coast  <- as.matrix(evengeo_coast)  
BER_coast  <- read.table ("geography_MATRIX_BER_coast.txt", sep = "\t")
BER_coast  <- as.matrix(BER_coast)  

####mantel tests of the fst matrices against the geographic distance
####crosses the ocean
odd16881Mantel <- mantel.test(odd16881, oddgeo_ocean, nperm = 1000 )
even16881Mantel <- mantel.test(even16881, evengeo_ocean, nperm = 1000)
all16881Mantel <- mantel.test(all16881, allgeo_ocean, nperm = 1000)
oddBERMantel <- mantel.test(oddBER, BER_ocean, nperm = 1000)
evenBERMantel <- mantel.test(evenBER, BER_ocean, nperm = 1000)

# Results
odd16881Mantel 
even16881Mantel
all16881Mantel
oddBERMantel
evenBERMantel

####coastline
odd16881Mantel_coast <- mantel.test(odd16881, oddgeo_coast, nperm = 1000 )
even16881Mantel_coast <- mantel.test(even16881, evengeo_coast, nperm = 1000)
all16881Mantel_coast <- mantel.test(all16881, allgeo_coast, nperm = 1000)
oddBERMantel_coast <- mantel.test(oddBER, BER_coast, nperm = 1000)
evenBERMantel_coast <- mantel.test(evenBER, BER_coast, nperm = 1000)

# Results
odd16881Mantel_coast 
even16881Mantel_coast
all16881Mantel_coast
oddBERMantel_coast
evenBERMantel_coast

######################################################################################################
#plotting IBD 

cols_OddEven <- c(O = "navy", L = "darkgrey", E = "magenta") 
cols_group <- c(Beringia = "darkorange", America = "purple") 

melted_fst_dist <- read.table ("PlottingR/IBD_linearFST_coast_ocean.txt", header = TRUE, sep = "\t")
head(melted_fst_dist)

#subset to dif version LFMM_e_LA<-Fst_mapped[Fst_mapped$even_LA=="yes",]
beringia <- melted_fst_dist[melted_fst_dist$location=="Beringia",]
odd <- melted_fst_dist[melted_fst_dist$Symbol=="O",]
even <- melted_fst_dist[melted_fst_dist$Symbol=="E",]

Ber_odd <- beringia[beringia$Symbol=="O",]
Ber_even <- beringia[beringia$Symbol=="E",]

######################################################################################################
#CODE TO PLOT LINNEAR FITS 
#https://susanejohnston.wordpress.com/2012/08/09/a-quick-and-easy-function-to-plot-lm-results-in-r/

#linear fits
beringia_line_coast <- lm(FST ~ coast, data = beringia)
odd_line_coast <- lm(FST ~ coast, data = odd)
even_line_coast <- lm(FST ~ coast, data = even)
Ber_odd_line_coast <- lm(FST ~ coast, data = Ber_odd)
Ber_even_line_coast <- lm(FST ~ coast, data = Ber_even)

beringia_line_ocean <- lm(FST ~ ocean, data = beringia)
odd_line_ocean <- lm(FST ~ ocean, data = odd)
even_line_ocean <- lm(FST ~ ocean, data = even)
Ber_odd_line_ocean <- lm(FST ~ ocean, data = Ber_odd)
Ber_even_line_ocean <- lm(FST ~ ocean, data = Ber_even)

ggplotRegression <- function (fit) {
  
  require(ggplot2)
  
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "red") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                       " Intercept =",signif(fit$coef[[1]],5 ),
                       " Slope =",signif(fit$coef[[2]], 5),
                       " P =",signif(summary(fit)$coef[2,4], 5)))
}


ggplotRegression(beringia_line_coast)
ggplotRegression(odd_line_coast)
ggplotRegression(even_line_coast)
ggplotRegression(Ber_odd_line_coast)
ggplotRegression(Ber_even_line_coast)

ggplotRegression(beringia_line_ocean)
ggplotRegression(odd_line_ocean)
ggplotRegression(even_line_ocean)
ggplotRegression(Ber_odd_line_ocean)
ggplotRegression(Ber_even_line_ocean)

######################################################################################################

#test for nice plot with lm 


p <- ggplot(data = melted_fst_dist) + geom_point(aes(x = coast, y = FST,  color = Symbol), alpha = .8, size = 4) + 
  scale_colour_manual(values = cols_OddEven, name = "Symbol") +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations IBD", x = "Coastline Distance (km)", y = "Pairwise Fst/(1-Fst)", size = 30))

p2 <- p + geom_smooth(data = melted_fst_dist[melted_fst_dist$Symbol=="O",], method = "lm")

######################################################################################################

ggplot(data = melted_fst_dist) + geom_point(aes(x = coast, y = FST,  color = Symbol), alpha = .8, size = 4) + 
  scale_colour_manual(values = cols_OddEven, name = "Symbol") +  
  #geom_smooth(data = melted_fst_dist[melted_fst_dist$Symbol=="O",], method = "lm") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
 labs(list(title= "All Populations IBD", x = "Coastline Distance (km)", y = "Pairwise Fst/(1-Fst)", size = 30))


ggplot(data = beringia) + geom_point(aes(x = coast, y = FST,  color = Symbol), alpha = .8, size = 4) + 
  scale_colour_manual(values = cols_OddEven, name = "Symbol") +  
  #geom_smooth(data = melted_fst_dist[melted_fst_dist$Symbol=="O",], method = "lm", se= FALSE) +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Beringia Populations IBD", x = "Coastline Distance (km)", y = "Pairwise Fst/(1-Fst)", size = 30))

ggplot(data = odd) + geom_point(aes(x = coast, y = FST,  color = location), alpha = .8, size = 4) + 
  scale_colour_manual(values = cols_group, name = "location") +  
  #geom_smooth(data = melted_fst_dist[melted_fst_dist$Symbol=="O",], method = "lm", se= FALSE) +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Odd Populations IBD", x = "Coastline Distance (km)", y = "Pairwise Fst/(1-Fst)", size = 30))


ggplot(data = even) + geom_point(aes(x = coast, y = FST,  color = location), alpha = .8, size = 4) + 
  scale_colour_manual(values = cols_group, name = "location") +  
  #geom_smooth(data = melted_fst_dist[melted_fst_dist$Symbol=="O",], method = "lm", se= FALSE) +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Even Populations IBD", x = "Coastline Distance (km)", y = "Pairwise Fst/(1-Fst)", size = 30))



















######################################################################################################


# ####waterways
# allgeo_water <- read.table ("geography_MATRIX_water.txt", sep = "\t")
# allgeo_water  <- as.matrix(allgeo_water)  
# oddgeo_water  <- read.table ("geography_MATRIX_odd_water.txt", sep = "\t")
# oddgeo_water  <- as.matrix(oddgeo_water)  
# evengeo_water  <- read.table ("geography_MATRIX_odd_water.txt", sep = "\t")
# evengeo_water  <- as.matrix(evengeo_water)


# ####WATERWAy
# odd16881Mantel_water <- mantel.test(odd16881, oddgeo_water, nperm = 999 )
# even16881Mantel_water <- mantel.test(even16881, evengeo_water, nperm = 999)
# all16881Mantel_water <- mantel.test(all16881, allgeo_water, nperm = 999)
# 
# # Results
# odd16881Mantel_water
# even16881Mantel_water
# all16881Mantel_water

######################### no duplicated loci

#making matrices of 15996 Fst
odd15996 <- read.table ("15996_fst_MATRIX_odd.txt", sep = "\t")
odd15996 <- as.matrix(odd15996)  
even15996 <- read.table ("15996_fst_MATRIX_even.txt", sep = "\t")
even15996 <- as.matrix(even15996)  
all15996 <- read.table ("15996_fst_MATRIX.txt", sep = "\t")
all15996 <- as.matrix(all15996)  


####mantel tests of the fst matrices against the geographic distance
####crosses the ocean
odd15996Mantel <- mantel.test(odd15996, oddgeo_ocean, nperm = 999 )
even15996Mantel <- mantel.test(even15996, evengeo_ocean, nperm = 999)
all15996Mantel <- mantel.test(all15996, allgeo_ocean, nperm = 999)

# Results
odd15996Mantel
even15996Mantel
all15996Mantel

####coastline
odd15996Mantel_coast <- mantel.test(odd15996, oddgeo_coast, nperm = 999 )
even15996Mantel_coast <- mantel.test(even15996, evengeo_coast, nperm = 999)
all15996Mantel_coast <- mantel.test(all15996, allgeo_coast, nperm = 999)

# Results
odd15996Mantel_coast
even15996Mantel_coast
all15996Mantel_coast

####WATERWAy
odd15996Mantel_water <- mantel.test(odd15996, oddgeo_water, nperm = 999 )
even15996Mantel_water <- mantel.test(even15996, evengeo_water, nperm = 999)
all15996Mantel_water <- mantel.test(all15996, allgeo_water, nperm = 999)

# Results
odd15996Mantel_water
even15996Mantel_water
all15996Mantel_water

