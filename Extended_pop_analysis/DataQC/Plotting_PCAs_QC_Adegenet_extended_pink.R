### Plotting a PCA for each pop with Adegenet
###   Data not completely filtered
### Carolyn Tarpey | November 2017
### ---------------------------------------


library("RColorBrewer")
library(ggplot2)
library(colorspace)
library(plyr)
library(colorRamps)
library(adegenet)
library(ade4)
library(genetics)
install.packages("genetics")
#install.packages("colorRamps")

###############DATA

#pop_key <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/Summary_Stats/POPINFO.txt", header = TRUE, sep = '\t')
#All_Genepop <- read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/FilteringGenotypes/PLINK/PCA/filter_inds_pca.eigenvec")
#head(ALL_eigenvec_table)
#ALL_tt = merge(ALL_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')
#ALL_geo <- ALL_tt[order(ALL_tt$Order_geo),] #sort by a geographical order, then odd then even 
#head(ALL_tt)


#PLOT INDIVIDUAL PCAs- Written by Wes

#########PCA functions
get_pca=function(filename="Something.gen")
{
  genepop_file <- read.genepop(file=filename)
  trans_data_no_missing <- scaleGen(genepop_file,NA.method="mean")
  ind_pca <- dudi.pca(trans_data_no_missing,cent=FALSE,scale=FALSE,scannf=FALSE,nf=4)
  return(ind_pca)    
}

All_Pops_PCA <- get_pca(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/Adegenet/batch_4_16662.genepop.gen")

Even_Pops_PCA <- get_pca(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/Adegenet/batch_4_16662_Even.genepop.gen")

Odd_Pops_PCA <- get_pca(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/Adegenet/batch_4_16662_Odd.genepop.gen")

NA_Pops_PCA <- get_pca(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/Adegenet/batch_4_16662_NA.genepop.gen")

Ber_Pops_PCA <- get_pca(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/Adegenet/batch_4_16662_BER.genepop.gen")

Old_Pops_PCA <- get_pca(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/Adegenet/batch_4_16662_Old.genepop.gen")

Old_383_Pops_PCA <- get_pca(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/Adegenet/batch_4_16662_Old_383.genepop.gen")

New_Pops_PCA <- get_pca(filename="Z:/WORK/TARPEY/Exp_Pink_Pops/Adegenet/batch_4_16662_New.genepop.gen")



#function to plot pca
plot_pca=function(pca_obj=ind_pca_WR_all_snps,pop_names=c("Anvil B","Yako B","Little Togiak R","Agulowak R","Teal C"," Hansen C")
                  , colors=c("dodgerblue","blue","forestgreen","green","firebrick1","firebrick4","purple","pink","gold1"), 
                  pop_samp_sizes=c(45,47,47,48,44,45),pcs_to_plot=c(1,2),xlim=c(-30,20),ylim=c(-20,20),
                  title="",rev_pc1=1,rev_pc2=1,plot_leg="TRUE",cex_leg=1,leg_loc="topright")
{
  eigenvects_to_plot=pca_obj$li
  variance_explained=(pca_obj$eig)/sum(pca_obj$eig)
  
  col_vect=NULL
  for(i in 1:length(colors))
  {
    rep_vec=rep(colors[i],pop_samp_sizes[i])
    col_vect=c(col_vect,rep_vec)
  }
  xlabel=paste("PC",pcs_to_plot[1]," ",round(variance_explained[pcs_to_plot[1]],3)*100,"% of variance",sep="")
  ylabel=paste("PC",pcs_to_plot[2]," ",round(variance_explained[pcs_to_plot[2]],3)*100,"% of variance",sep="")
  
  plot(rev_pc1*eigenvects_to_plot[,pcs_to_plot[1]],rev_pc2*eigenvects_to_plot[,pcs_to_plot[2]],
       bg=col_vect, pch=21, col="black",xlab=xlabel,ylab=ylabel,
       xlim=xlim, ylim=ylim)
  abline(h=0)
  abline(v=0)
  title(title)
  if(plot_leg=="TRUE"){legend(leg_loc,pop_names,cex=cex_leg, fill=colors,bty="n")}
  return(eigenvects_to_plot)
}

###Plots

pca_All_pops <- plot_pca(pca_obj=All_Pops_PCA, xlim=c(-75,150), ylim=c(-50,75),
                         pop_names=c("PAMUR10","PAMUR11","PSUSIT13","PHAYLY09",
                                     "PHAYLY10","PKOPE91","PKOPE96","PKUSHI06","PKUSHI07","PLAKEL06",
                                     "PLAKEL07","PNOME91","PNOME94","PSNOH03","PSNOH96","PSUSIT14","PTAUY09","PTAUY12"),
                         pop_samp_sizes=c(32,32,24,32,32,24,24,32,32,24,24,24,24,20,24,24,24,32,32),
                         colors=c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", 
                                  "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", 
                                  "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", 
                                  "#771122", "#AA4455", "#DD7788"),
                         plot_leg="TRUE",cex_leg=0.7,leg_loc="bottomright",pcs_to_plot=c(1,2))

pca_Even_pops <- plot_pca(pca_obj=Even_Pops_PCA, xlim=c(-75,150), ylim=c(-50,75),
                      pop_names=c("PAMUR10","PHAYLY10","PKOPE96","PKUSHI06","PLAKEL06",
                                  "PNOME94","PSNOH96","PSUSIT14","PTAUY12"),
                      pop_samp_sizes=c(32,32,24,32,24,20,24,24,32),
                      colors=c("blue","grey","purple","yellow","green","light green","red","orange","magenta"),
                      plot_leg="TRUE",cex_leg=0.7,leg_loc="bottomright",pcs_to_plot=c(1,2))

pca_Odd_pops <- plot_pca(pca_obj=Odd_Pops_PCA, xlim=c(-75,150), ylim=c(-50,75),
                          pop_names=c("PAMUR10","PHAYLY10","PKOPE96","PKUSHI06","PLAKEL06",
                                      "PNOME94","PSNOH96","PSUSIT14","PTAUY12"),
                          pop_samp_sizes=c(32,32,24,32,24,20,24,24,32),
                          colors=c("blue","grey","purple","yellow","green","light green","red","orange","magenta"),
                          plot_leg="TRUE",cex_leg=0.7,leg_loc="bottomright",pcs_to_plot=c(1,2))

pca_NA_pops <- plot_pca(pca_obj=NA_Pops_PCA, xlim=c(-75,150), ylim=c(-50,75),
                          pop_names=c("PAMUR10","PHAYLY10","PKOPE96","PKUSHI06","PLAKEL06",
                                      "PNOME94","PSNOH96","PSUSIT14","PTAUY12"),
                          pop_samp_sizes=c(32,32,24,32,24,20,24,24,32),
                          colors=c("blue","grey","purple","yellow","green","light green","red","orange","magenta"),
                          plot_leg="TRUE",cex_leg=0.7,leg_loc="bottomright",pcs_to_plot=c(1,2))

pca_Ber_pops <- plot_pca(pca_obj=Ber_Pops_PCA, xlim=c(-75,150), ylim=c(-50,60),
                          pop_names=c("PAMUR10","PHAYLY10","PKOPE96","PKUSHI06","PLAKEL06",
                                      "PNOME94","PSNOH96","PSUSIT14","PTAUY12"),
                          pop_samp_sizes=c(32,32,24,32,24,20,24,24,32),
                          colors=c("blue","grey","purple","yellow","green","light green","red","orange","magenta"),
                          plot_leg="TRUE",cex_leg=0.7,leg_loc="bottomright",pcs_to_plot=c(1,2))


pca_Old_pops <- plot_pca(pca_obj=Old_Pops_PCA, xlim=c(-75,100), ylim=c(-50,40),
                          pop_names=c("PAMUR10","PAMUR11","PHAYLY09",
                                      "PHAYLY10","PKOPE91","PKOPE96","PKUSHI06","PKUSHI07",
                                      "PNOME91","PNOME94","PSNOH03","PSNOH96","PTAUY09","PTAUY12"),
                          pop_samp_sizes=c(32,32,32,32,24,24,32,32,24,24,20,24,24,32,32),
                          colors=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", 
                                   "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C"),
                        plot_leg="TRUE",cex_leg=0.7,leg_loc="bottomright",pcs_to_plot=c(2,3))

pca_Old_383_pops <- plot_pca(pca_obj=Old_383_Pops_PCA, xlim=c(-100,100), ylim=c(-100,140),
                         pop_names=c("PAMUR10","PAMUR11","PHAYLY09",
                                     "PHAYLY10","PKOPE91","PKOPE96","PKUSHI06","PKUSHI07",
                                     "PNOME91","PNOME94","PSNOH03","PSNOH96","PTAUY09","PTAUY12"),
                         pop_samp_sizes=c(32,32,32,32,24,24,32,32,24,24,20,24,24,32,32),
                         colors=c("#882E72", "#B178A6", "#D6C1DE", "#1965B0", "#5289C7", "#7BAFDE", 
                                  "#4EB265", "#90C987", "#CAE0AB", "#F7EE55", "#F6C141", "#F1932D", "#E8601C", "#DC050C"),
                         plot_leg="TRUE",cex_leg=0.7,leg_loc="bottomright",pcs_to_plot=c(2,3))

#Scatter plots

scatter(Even_Pops_PCA, scree.da=FALSE, bg="white", pch=14:22,  cell=0, cstar=0, solid=.4, clab=0,
        cex=3,  leg=TRUE, txt.leg=c("PAMUR10","PHAYLY10","PKOPE96","PKUSHI06","PLAKEL06",
                                        "PNOME94","PSNOH96","PSUSIT14","PTAUY12"))

scatter(Odd_Pops_PCA, scree.da=FALSE, bg="white", pch=14:22,  cell=0, cstar=0, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=c("PAMUR11","PSUSIT13","PHAYLY09","PKOPE91","PKUSHI07",
                                          "PLAKEL07","PNOME91","PSNOH03","PTAUY09"))

scatter(NA_Pops_PCA, scree.da=FALSE, bg="white", pch=15:22,  cell=0, cstar=0, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=c("PSUSIT13","PKOPE91","PKOPE96","PLAKEL06",
                                          "PLAKEL07","PSNOH03","PSNOH96","PSUSIT14","PTAUY09","PTAUY12"))

scatter(Ber_Pops_PCA, scree.da=FALSE, bg="white", pch=13:22,  cell=0, cstar=0, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=("PAMUR10","PAMUR11","PHAYLY09","PHAYLY10","PKUSHI06",
                                         "PKUSHI07","PNOME91","PNOME94","PTAUY09","PTAUY12"))

scatter(All_Pops_PCA, scree.da=FALSE, bg="white", pch=7:22,  cell=0, cstar=0, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=c("PAMUR10","PAMUR11","PSUSIT13","PHAYLY09",
                                          "PHAYLY10","PKOPE91","PKOPE96","PKUSHI06","PKUSHI07","PLAKEL06",
                                          "PLAKEL07","PNOME91","PNOME94","PSNOH03","PSNOH96","PSUSIT14","PTAUY09","PTAUY12"))

scatter(Old_Pops_PCA, scree.da=FALSE, bg="white", pch=9:22,  cell=0, cstar=0, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=c("PAMUR10","PAMUR11","PHAYLY09",
                                          "PHAYLY10","PKOPE91","PKOPE96","PKUSHI06","PKUSHI07",
                                          "PNOME91","PNOME94","PSNOH03","PSNOH96","PTAUY09","PTAUY12"))

scatter(New_Pops_PCA, scree.da=FALSE, bg="white", pch=8:22,  cell=0, cstar=0, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=c("PSUSIT13","PLAKEL06", "PLAKEL07","PSUSIT14"))

##################################

### Plotting the PCA from PLINK
###   Final versions 
### Carolyn Tarpey | April 2016 
### ---------------------------------------

library("RColorBrewer")
library(ggplot2)
library(colorspace)
library(plyr)
library(colorRamps)
#install.packages("colorspace")
#install.packages("colorRamps")
###############DATA

pop_key <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/POPINFO_original.txt", header = TRUE, sep = '\t')

ALL_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Original_populations/16681_80_pca.eigenvec")
ALL_tt = merge(ALL_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')
ALL_geo <- ALL_tt[order(ALL_tt$Order_geo),] #sort by a geographical order, then odd then even 

Beringia_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Original_populations/Beringia_pops_pca.eigenvec")
Beringia_tt = merge(Beringia_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')
Beringia_geo <- Beringia_tt[order(Beringia_tt$Order_geo),] #sort by a geographical order, then odd then even 

Asia_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Original_populations/Asia_pops_pca.eigenvec")
Asia_tt = merge(Asia_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')
Asia_geo <- Asia_tt[order(Asia_tt$Order_geo),] #sort by a geographical order, then odd then even 

NorthAmerica_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Original_populations/NorthAmerica_pops_pca.eigenvec")
NorthAmerica_tt = merge(NorthAmerica_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')
NorthAmerica_geo <- NorthAmerica_tt[order(NorthAmerica_tt$Order_geo),] #sort by a geographical order, then odd then even 

Beringia_Odd_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Original_populations/Beringia_Odd_pops_pca.eigenvec")
Beringia_Odd_tt = merge(Beringia_Odd_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')
#ALL_geo <- ALL_tt[order(ALL_tt$Order_geo),] #sort by a geographical order, then odd then even 

Beringia_Even_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Original_populations/Beringia_Even_pops_pca.eigenvec")
Beringia_Even_tt = merge(Beringia_Even_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')
#ALL_geo <- ALL_tt[order(ALL_tt$Order_geo),] #sort by a geographical order, then odd then even 

Odd_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Original_populations/Odd_pops_pca.eigenvec")
Odd_tt = merge(Odd_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')
#ALL_geo <- ALL_tt[order(ALL_tt$Order_geo),] #sort by a geographical order, then odd then even 

Even_eigenvec_table <- read.table("C:/Users/Carolyn/Desktop/16681_80_PLINK/Original_populations/Even_pops_pca.eigenvec")
Even_tt = merge(Even_eigenvec_table,pop_key, by.x = 'V1', by.y = 'CLUSTER')
#ALL_geo <- ALL_tt[order(ALL_tt$Order_geo),] #sort by a geographical order, then odd then even 

names(ALL_tt)

#####COLORS
cols_OddEven <-c(AMUR10 = "#EE547D", AMUR11 = "#4E70C8", HAYLY09 = "#4C7DE7", HAYLY10 = "#ED2DBE", 
                 KOPPE96 = "#D186CF", KOPPE91 = "#B0B0D8", KUSHI06 = "#AD2868", KUSHI07 = "#354A7C",
                 NOME91 = "#95E0CA", NOME94 = "#E63D25", SNOH03 = "#708FEC", SNOH96 ="#a12787" , TAUY09 = "#3C92A8", TAUY12 = "#DF9E39")

pal<- choose_palette()


## These are pretty , but the mid colors are too light 
#cols_ord7 <- c(AMUR10 = "#fc8d59", AMUR11 ="#fc8d59",  HAYLY09 = "#ffffbf", HAYLY10 = "#ffffbf", 
#                KOPPE96 = "#91bfdb", KOPPE91 = "#91bfdb", KUSHI06 = "#d73027", KUSHI07 = "#d73027",
#                NOME91 = "#e0f3f8", NOME94 = "#e0f3f8", SNOH03 = "#4575b4", SNOH96 ="#4575b4" , TAUY09 = "#fee090", TAUY12 = "#fee090")
# #d73027  #fc8d59  #fee090  #ffffbf  #e0f3f8  #91bfdb  #4575b4


## These are harsh and still not good colorRamps::blue2red(7) 
#cols_ord7 <- c(AMUR10 = "#0055FF", AMUR11 ="#0055FF",  HAYLY09 = "#00FFFF", HAYLY10 = "#00FFFF", 
#               KOPPE96 = "#FF5500", KOPPE91 = "#FF5500", KUSHI06 = "#0000FF", KUSHI07 = "#0000FF",
#               NOME91 = "#FFAA00", NOME94 = "#FFAA00", SNOH03 = "#FF0000", SNOH96 ="#FF0000" , TAUY09 = "#00AAFF", TAUY12 = "#00AAFF")

#"#0000FF" "#0055FF" "#00AAFF" "#00FFFF" "#FFAA00" "#FF5500" "#FF0000"

## These are too close together colorRampPalette(c("red", "blue"))(7)
#cols_ord7 <- c(AMUR10 = "#D4002A", AMUR11 ="#D4002A",  HAYLY09 = "#7F007F", HAYLY10 = "#7F007F", 
#               KOPPE96 = "#2A00D4", KOPPE91 = "#2A00D4", KUSHI06 = "#FF0000", KUSHI07 = "#FF0000",
#               NOME91 = "#5500AA", NOME94 = "#5500AA", SNOH03 = "#0000FF", SNOH96 ="#0000FF" , TAUY09 = "#AA0055", TAUY12 = "#AA0055")
#"#FF0000" "#D4002A" "#AA0055" "#7F007F" "#5500AA" "#2A00D4" "#0000FF"

#shapes_OddEven <- c(6, 1, 1, 6, 6, 1, 6, 1, 1, 6, 1, 6, 1, 6)

#scale_shape_discrete(name ="Population", breaks=c("AMUR10", "AMUR11", "HAYLY09", "HAYLY10", 
#  "KOPPE96", "KOPPE91", "KUSHI06", "KUSHI07","NOME91", "NOME94", "SNOH03", "SNOH96", "TAUY09", "TAUY12"),
#labels=c("Amur even", "Amur odd", "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", 
#  "Prince William Sound even", "Hokkaido even", "Hokkaido odd", "Norton Sound odd", "Norton Sound even", 
#  "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")

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




###@@@@@@@@@@@@@@@@@@@@@@Plots for the Manuscript, have labels as Principal components

###ALL 
pdf("G:/Analysis/Pop_analysis/Populations_b3_may/PCA_R_images_plink/Updated PCA plots_4manuscript.pdf", width = 9, height = 7)

ggplot(data = ALL_geo) + geom_point(aes(x = V3, y = V4,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = ALL_col, name ="Population", labels = c("Amur even", "Amur odd",
                                                                       "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
                                                                       "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Principal Components 1 & 2", x = "Principal Component 1 (30.24%)", y = "Principal Component 2 (13.83%)", size = 30))

##@@@@@@@@@@@@@@@@@@@@@@  one and three

###ALL 

ggplot(data = ALL_geo) + geom_point(aes(x = V3, y = V5,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = ALL_col, name ="Population", labels = c("Amur even", "Amur odd",
                                                                       "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
                                                                       "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Principal Components 1 & 3", x = "Principal Component 1 (30.24%)", y = "Principal Component 3 (8.73%)", size = 30))

##@@@@@@@@@@@@@@@@@@@@@@ two and three

###ALL 

ggplot(data = ALL_geo) + geom_point(aes(x = V4, y = V5,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = ALL_col, name ="Population", labels = c("Amur even", "Amur odd",
                                                                       "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
                                                                       "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Principal Components 2 & 3", x = "Principal Component 2 (13.83%)", y = "Principal Component 3 (8.73%)", size = 30))

dev.off()



###################################################These use the term dimension instead of Principal Component

###@@@@@@@@@@@@@@@@@@@@@@dimmension one and two

###ALL 
pdf("G:/Analysis/Pop_analysis/Populations_b3_may/PCA_R_images_plink/PDF of updated PCA plots_4.pdf", width = 9, height = 7)

ggplot(data = ALL_geo) + geom_point(aes(x = V3, y = V4,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = ALL_col, name ="Population", labels = c("Amur even", "Amur odd",
                                                                       "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
                                                                       "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Dimensions 1 & 2", x = "Dimension 1 (30.24%)", y = "Dimension 2 (13.83%)", size = 30))


###Beringia_tt

ggplot(data = Beringia_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = Ber_col, name ="Population", labels = c("Amur even", "Amur odd",
                                                                      "Kamchatka odd", "Kamchatka even", "Hokkaido even",
                                                                      "Hokkaido odd", "Norton Sound odd", "Norton Sound even",  "Magadan odd","Magadan even")) +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( title= "Beringia Dimensions 1 & 2", x = "Dimension 1 (33.75%)", y = "Dimension 2 (7.16%)", size = 30))


###Asia_tt

ggplot(data = Asia_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_color_manual(values = Asia_col, name ="Population", labels = c("Amur even", "Amur odd", "Kamchatka odd", "Kamchatka even", "Hokkaido even", "Hokkaido odd","Magadan odd","Magadan even"))+ 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Asia Dimensions 1 & 2", x = "Dimension 1 (32.10%)", y = "Dimension 2 (6.83%)", size = 30))

###NorthAmerica_tt

ggplot(data = NorthAmerica_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_color_manual(values = NA_col, name ="Population", labels = c("Prince William Sound odd", "Prince William Sound even", 
                                                                     "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even")) +  scale_x_continuous(trans="reverse")+
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "North America Dimensions 1 & 2", x = "Dimension 1 (22.76%)", y = "Dimension 2 (10.36%)", size = 30)) 


###@@@@@@@@@@@@@@@@@@@@@@dimmension one and three

###ALL 

ggplot(data = ALL_geo) + geom_point(aes(x = V3, y = V5,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = ALL_col, name ="Population", labels = c("Amur even", "Amur odd",
                                                                       "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
                                                                       "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Dimensions 1 & 3", x = "Dimension 1 (30.24%)", y = "Dimension 3 (8.73%)", size = 30))


###Beringia_tt

ggplot(data = Beringia_geo) + geom_point(aes(x = V3, y = V5, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = Ber_col, name ="Population", labels = c("Amur even", "Amur odd",
                                                                      "Kamchatka odd", "Kamchatka even", "Hokkaido even", "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Magadan odd","Magadan even")) +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(), legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "Beringia Dimensions 1 & 3", x = "Dimension 1 (33.75%)", y = "Dimension 3 (6.18%)", size = 30))


###Asia_tt

ggplot(data = Asia_geo) + geom_point(aes(x = V3, y = V5, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_color_manual(values = Asia_col, name ="Population", labels = c("Amur even", "Amur odd",
                                                                       "Kamchatka odd", "Kamchatka even", "Hokkaido even", "Hokkaido odd",  "Magadan odd","Magadan even")) + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "Asia Dimensions 1 & 3", x = "Dimension 1 (32.10%)", y = "Dimension 3 (6.04%)", size = 30))

###NorthAmerica_tt

ggplot(data = NorthAmerica_geo) + geom_point(aes(x = V3, y = V5, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_color_manual(values = NA_col, name ="Population", labels = c("Prince William Sound odd", "Prince William Sound even", 
                                                                     "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even")) +  scale_x_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "North America Dimensions 1 & 3", x = "Dimension 1 (22.76%)", y = "Dimension 3 (6.73%)", size = 30)) 


###@@@@@@@@@@@@@@@@@@@@@@dimmension two and three

###ALL 

ggplot(data = ALL_geo) + geom_point(aes(x = V4, y = V5,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = ALL_col, name ="Population", labels = c("Amur even", "Amur odd",
                                                                       "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
                                                                       "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Dimensions 2 & 3", x = "Dimension 2 (13.83%)", y = "Dimension 3 (8.73%)", size = 30))


###Beringia_tt

ggplot(data = Beringia_geo) + geom_point(aes(x = V4, y = V5, color = POPNAME), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = Ber_col, name ="Population", labels = c("Amur even", "Amur odd",
                                                                      "Kamchatka odd", "Kamchatka even", "Hokkaido even", "Hokkaido odd", "Norton Sound odd", "Norton Sound even", 
                                                                      "Magadan odd","Magadan even")) +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "Beringia Dimensions 2 & 3", x = "Dimension 2 (7.16%)", y = "Dimension 3 (6.18%)", size = 30))


###Asia_tt

ggplot(data = Asia_geo) + geom_point(aes(x = V4, y = V5, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_color_manual(values = Asia_col, name ="Population", labels = c("Amur even", "Amur odd",
                                                                       "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
                                                                       "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "Asia Dimensions 2 & 3", x = "Dimension 2 (6.83%)", y = "Dimension 3 (6.04%)", size = 30))

###NorthAmerica_tt

ggplot(data = NorthAmerica_geo) + geom_point(aes(x = V4, y = V5, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_color_manual(values = NA_col, name ="Population", labels = c("Amur even", "Amur odd",
                                                                     "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
                                                                     "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) +  scale_x_continuous(trans="reverse")+
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "North America Dimensions 2 & 3", x = "Dimension 2 (10.36%)", y = "Dimension 3 (6.73%)", size = 30)) 



###@@@@@@@@@@@@@@@@@@@@@@dimmension  two and four

###ALL 

ggplot(data = ALL_tt) + geom_point(aes(x = V4, y = V6,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = cols_OddEven, name ="Population", labels = c("Amur even", "Amur odd",
                                                                            "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
                                                                            "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Dimensions 2 & 4", x = "Dimension 2 (13.83%)", y = "Dimension 4 (5.13%)", size = 30))


###Beringia_tt

ggplot(data = Beringia_tt) + geom_point(aes(x = V4, y = V6, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = cols_OddEven, name ="Population", labels = c("Amur even", "Amur odd",
                                                                           "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
                                                                           "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "Beringia Dimensions 2 & 4", x = "Dimension 2 (7.16%)", y = "Dimension 4 (3.54%)", size = 30))


###Asia_tt

ggplot(data = Asia_tt) + geom_point(aes(x = V4, y = V6, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population", labels = c("Amur even", "Amur odd",
                                                                           "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
                                                                           "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "Asia Dimensions 2 & 4", x = "Dimension 2 (6.83%)", y = "Dimension 4 (3.40%)", size = 30))

###NorthAmerica_tt

ggplot(data = NorthAmerica_tt) + geom_point(aes(x = V4, y = V6, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population", labels = c("Amur even", "Amur odd",
                                                                           "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
                                                                           "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) +  scale_x_continuous(trans="reverse")+
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title = "North America Dimensions 2 & 4", x = "Dimension 2 (10.36%)", y = "Dimension 4 (5.88%)", size = 30)) 


###@@@@@@@@@@@@@@@@@@@@@@  dimmension three and four

###ALL 


ggplot(data = ALL_tt) + geom_point(aes(x = V5, y = V6,  color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_colour_manual(values = cols_OddEven, name ="Population", labels = c("Amur even", "Amur odd",
                                                                            "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
                                                                            "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "All Populations Dimensions 3 & 4", x = "Dimension 3 (8.73%)", y = "Dimension 4 (5.13%)", size = 30))


###Beringia_tt

ggplot(data = Beringia_tt) + geom_point(aes(x = V5, y = V6, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  scale_color_manual(values = cols_OddEven, name ="Population", labels = c("Amur even", "Amur odd",
                                                                           "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
                                                                           "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) +  
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Beringia Dimensions 3 & 4", x = "Dimension 3 (6.18%)", y = "Dimension 4 (3.54%)", size = 30))


###Asia_tt

ggplot(data = Asia_tt) + geom_point(aes(x = V5, y = V6, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population", labels = c("Amur even", "Amur odd",
                                                                           "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
                                                                           "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) + 
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "Asia Dimensions 3 & 4", x = "Dimension 3 (6.04%)", y = "Dimension 4 (3.40%)", size = 30))

###NorthAmerica_tt

ggplot(data = NorthAmerica_tt) + geom_point(aes(x = V5, y = V6, color = POPNAME, shape = LINEAGE), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population", labels = c("Amur even", "Amur odd",
                                                                           "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", "Prince William Sound even", "Hokkaido even",
                                                                           "Hokkaido odd", "Norton Sound odd", "Norton Sound even", "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) +  scale_x_continuous(trans="reverse")+
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),legend.text = element_text(size= 17),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list(title= "North America Dimensions 3 & 4", x = "Dimension 3 (6.73%)", y = "Dimension 4 (5.88%)", size = 30)) 

dev.off()









###############################  Lineages

###Beringia_Even_tt

ggplot(data = Beringia_Even_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))


###Beringia_ODD_tt

ggplot(data = Beringia_Odd_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))

###Even_tt

ggplot(data = Even_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .9, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))

###Odd_tt

ggplot(data = Odd_tt) + geom_point(aes(x = V3, y = V4, color = POPNAME), alpha = .8, size = 4) + 
  scale_color_manual(values = cols_OddEven, name ="Population") +  
  scale_x_continuous(trans="reverse") + scale_y_continuous(trans="reverse") +
  theme_classic() + theme(text = element_text(size= 20), legend.title = element_blank(),
                          axis.line.x = element_line(), axis.line.y = element_line()) +
  labs(list( x = "Dimension 1", y = "Dimension 2", size = 30))


# 
# scale_shape_discrete(name ="Population", breaks=c("AMUR10", "AMUR11", "HAYLY09", "HAYLY10", 
#                                                   "KOPPE96", "KOPPE91", "KUSHI06", "KUSHI07","NOME91", "NOME94", "SNOH03", "SNOH96", "TAUY09", "TAUY12"),
#                      labels=c("Amur even", "Amur odd", "Kamchatka odd", "Kamchatka even", "Prince William Sound odd", 
#                               "Prince William Sound even", "Hokkaido even", "Hokkaido odd", "Norton Sound odd", "Norton Sound even", 
#                               "Puget Sound odd", "Puget Sound even", "Magadan odd","Magadan even")) +


#____________________________________________________________________________________________________________
## PCA subsetting the genind file

# gen16662<-import2genind("Z:/WORK/TARPEY/Exp_Pink_Pops/Adegenet/batch_4_16662.genepop.gen")
table(pop(gen16662))
#barplot(table(pop(gen16662)), col=funky(18), xpd=FALSE, xlab="Population", ylab="Sample size")

pop_samplesize<- table(pop(gen16662)
temp_summary<- summary(gen16662)

#DAPC of different groups 
split_pops<- seppop(gen16662)
names(split_pops)

# "PAMUR10_0032"             "PAMUR11_0033_comb"        "PDISAP13_0095.1.combined"
# "PHAYLY09_0043_comb"       "PHAYLY10_0046_comb"       "PKOPE91T_0037"           
# "PKOPE96T_0024"            "PKUSHI06_0032_comb"       "PKUSHI07_0032_comb"      
# "PLAKEL06_0101.1.combined" "PLAKEL07_0200.1.combined" "PNOME91_0033"            
# "PNOME94_0020"             "PSNOH03_0095"             "PSNOH96_0028"            
# "PSPINK14_0052.1.combined" "PTAUY09_0046"             "PTAUY12_0042"  

## dAPC of All
dapc1<- dapc(gen16662)
scatter(dapc1)

scatter(dapc1, scree.da=FALSE, bg="white", pch=20,  cell=0, cstar=0, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=c("PAMUR10","PAMUR11","PDISAP13","PHAYLY09",
                                          "PHAYLY10","PKOPE91","PKOPE96","PKUSHI06","PKUSHI07","PLAKEL06",
                                          "PLAKEL07","PNOME91","PNOME94","PSNOH03","PSNOH96","PSPINK14","PTAUY09","PTAUY12"))
## dAPC of Old Pops 

old_pops <- repool(split_pops[c( "PAMUR10_0032", "PAMUR11_0033_comb","PHAYLY09_0043_comb","PHAYLY10_0046_comb",
                                 "PKOPE91T_0037", "PKOPE96T_0024","PKUSHI06_0032_comb","PKUSHI07_0032_comb","PNOME91_0033",            
                                 "PNOME94_0020", "PSNOH03_0095","PSNOH96_0028","PTAUY09_0046","PTAUY12_0042")])
dapc_old_pops<- dapc(old_pops)
scatter(dapc_old_pops)

scatter(dapc_old_pops, scree.da=FALSE, bg="white", pch=20,  cell=0, cstar=0, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=c("PAMUR10","PAMUR11","PHAYLY09",
                                          "PHAYLY10","PKOPE91","PKOPE96","PKUSHI06","PKUSHI07",
                                          "PNOME91","PNOME94","PSNOH03","PSNOH96","PTAUY09","PTAUY12"))


