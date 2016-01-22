### Identify centromere location with phased haploid data: A more recent edit
###    Written by Morten Limborg edited for my data by me
### Carolyn Tarpey | January 2016 
### ---------------------------------------


setwd("G:/Analysis/Mapping/AllHaps/Centromeres" )


######################################################################
#######Place centromeres and count y distribution with haploids#######
######################################################################


rm(list=ls())


plot.pseudo.y <- function(LG, plot){
  
  ##Read in map with phased genotypes:

  #Map <- as.matrix(read.csv(file = 'Cents_05.txt', sep = '\t')) #Read map w. haploid genotypes
  #Map <- as.matrix(read.csv(file = 'Cents_01.txt', sep = '\t')) #Read map w. haploid genotypes
  #Map <- as.matrix(read.csv(file = 'Cents_110.txt', sep = '\t')) #Read map w. haploid genotypes
  Map <- as.matrix(read.csv(file = 'Cents_108.txt', sep = '\t'))  #Read map w. haploid genotypes

is.matrix(Map)
  rownames(Map) <- Map[,1]                #Change column 1 into rownames
  Map <- Map[,2:ncol(Map)]                
  
  sample.size <- nrow(Map)-5
  
  #Making subset matrix for LG_X:
  #LG <- 17 #For Testing
  
  LG_X <- NULL
  for (i in 1:length(Map[1,])){
    if (Map[1,i] == LG){  
      LG_X <- cbind(LG_X, Map[,i])
    }
  }
  
  ##Make new matrix where a & b phases are recoded to 1 and 0:
  XO.vec.kw <- matrix(nrow=nrow(Map), ncol=ncol(LG_X))  #Make Cross-over (XO) vector (vec) using kernel-window (kw) for focal LG
  XO.vec.kw[1:5,] <- LG_X[1:5,]
  rownames(XO.vec.kw) <- rownames(LG_X)
  
  for (j in 6:nrow(XO.vec.kw)){
    for (i in 1:ncol(XO.vec.kw)){
      if (as.character(LG_X[j,i]) == "a"){
        XO.vec.kw[j,i] <- 0  
      }
      if (as.character(LG_X[j,i]) == "b"){
        XO.vec.kw[j,i] <- 1  
      }
    }
  }
  
  #Impute first and last genotype if missing:
  for (j in 6:nrow(XO.vec.kw)){
    first.ten.phases <- as.numeric(XO.vec.kw[j,1:10])
    start.phase <- round(mean(first.ten.phases, na.rm=T), digits=0)
    if (as.character(LG_X[j,1]) != start.phase){
      XO.vec.kw[j,1] <- start.phase  
    }
    
    last.ten.phases <-  as.numeric(XO.vec.kw[j,(ncol(XO.vec.kw)-9):ncol(XO.vec.kw)])
    end.phase <- round(mean(last.ten.phases, na.rm=T), digits=0)
    if (as.character(LG_X[j,ncol(XO.vec.kw)]) != end.phase){
      XO.vec.kw[j,ncol(XO.vec.kw)] <- end.phase  
    }
  } 
      
  #Make sub-matrix with only genotypes and change NA to "-":
  XO.phased.genot <- XO.vec.kw[6:nrow(Map),]
  
  for (j in 1:nrow(XO.phased.genot)){
    for (i in 1:ncol(XO.phased.genot)){
      if (is.na(XO.phased.genot[j,i])){
        XO.phased.genot[j,i] <- "-"  
      }
    }
  }
  
  #Sort matrix by first marker and count number of "0"
  XO.phased.genot <- as.data.frame(XO.phased.genot)
  colnames(XO.phased.genot) <- paste("marker",c(1:ncol(XO.phased.genot)), sep="") 
  XO.phased.genot <- XO.phased.genot[order(XO.phased.genot$marker1),] #Sorted by first marker
  XO.phased.genot.right <- XO.phased.genot[order(XO.phased.genot[,ncol(XO.phased.genot)]),] #Sorted by last marker
  
  #Find number of individuals having 0 at first marker
  no.zero <- sum(XO.phased.genot$marker1 == 0)
  no.zero.right <- sum(XO.phased.genot[,ncol(XO.phased.genot)] == 0)
  
  #Convert XO.phased.genot back to matrix in order to edit numeric contents
  XO.phased.genot <- as.matrix(XO.phased.genot)
  XO.phased.genot.right <- as.matrix(XO.phased.genot.right)
  #is.matrix(XO.phased.genot.right)
  
  
  ##Make new temporary matrix for the ~50% offspring that needs to get phases switched (i.e. having 0 at marker1)
  
  #From left:
  Matrix.convert <- NULL
  converted.genot <- NULL #Temporary vector for each individual to be used within loop
  
  if (no.zero < sample.size){
    for (j in (no.zero+1):nrow(XO.phased.genot)){  #Loop over individuals having the genotype 1 at the first marker
      for (i in 1:ncol(XO.phased.genot)){ #For individual j: Loop over all markers in LG_X
        converted.genot <- cbind(converted.genot, abs(as.numeric(XO.phased.genot[j,i])-1))
        
      }
      Matrix.convert <- rbind(Matrix.convert, converted.genot)
      converted.genot <- NULL
    }
  }
  #Else assign entire XO.phased.genot.right matrix -> Matrix.convert.right
  if (no.zero >= sample.size){
    Matrix.convert <- NULL
  }


  #From right:
  Matrix.convert.right <- NULL
  converted.genot.right <- NULL #Temporary vector for each individual to be used within loop
  
  if (no.zero.right < sample.size){
    for (j in (no.zero.right+1):nrow(XO.phased.genot.right)){  #Loop over individuals having the genotype 1 at the last marker
      for (i in 1:ncol(XO.phased.genot.right)){ #For individual j: Loop over all markers in LG_X
        converted.genot.right <- cbind(converted.genot.right, abs(as.numeric(XO.phased.genot.right[j,i])-1))
        
      }
      Matrix.convert.right <- rbind(Matrix.convert.right, converted.genot.right)
      converted.genot.right <- NULL
    }
  }
  #Else assign empty matrix -> Matrix.convert.right
  if (no.zero.right >= sample.size){
    Matrix.convert.right <- NULL
  }
    

  ##Merge synced phased offspring into final matrix for counting cross-over frequencies
  #From left:
  phased.matrix <- rbind(XO.phased.genot[1:no.zero,], Matrix.convert)
  rownames(phased.matrix) <- rownames(XO.phased.genot)
  
  #Make vector for pseudo y values for LG_X from left
  pseudo.y.LG_X <- NULL
  
  for (i in 1:ncol(phased.matrix)){
    pseudo.y <- sum(as.numeric(phased.matrix[,i]), na.rm = TRUE)/sample.size
    pseudo.y.LG_X <- cbind(pseudo.y.LG_X, pseudo.y)
  }
  
  #From right:
  phased.matrix.right <- rbind(XO.phased.genot.right[1:no.zero.right,], Matrix.convert.right)
  rownames(phased.matrix.right) <- rownames(XO.phased.genot.right)
  
  #Make vector for pseudo y values for LG_X from right
  pseudo.y.LG_X.right <- NULL
  
  for (i in 1:ncol(phased.matrix.right)){
    pseudo.y.right <- sum(as.numeric(phased.matrix.right[,i]), na.rm = TRUE)/sample.size
    pseudo.y.LG_X.right <- cbind(pseudo.y.LG_X.right, pseudo.y.right)
  }
    
  #Plot pseudo-y values from left and right:
  if (plot == "Y"){
  plot(NULL, pch=19, col="red", xlab="cM", ylab="RF (phase-shift)", xlim=c(0,150), ylim=c(0,1.2), las=1, main=paste("LG",LG))
  lines(x=c(-10,150), y=c(0.5, 0.5), lty=2)  
  lines(x=c(-10,150), y=c(1, 1), lty=2, col="grey")  
#  points(x=as.numeric(LG_X[2,]), y=1-pseudo.y.LG_X*2, pch=19, col="red")
#  points(x=as.numeric(LG_X[2,]), y=1-pseudo.y.LG_X.right*2, pch=19, col="blue")
  
  points(x=as.numeric(LG_X[2,]), y=pseudo.y.LG_X, pch=19, col="red")
  points(x=as.numeric(LG_X[2,]), y=pseudo.y.LG_X.right, pch=19, col="blue")
  
#  points(x=as.numeric(LG_X[2,]), y=LG_X[5,], pch=19, col="black") # y values
  }
  
  #Making output matrix with pseudo y values
  output <- rbind(as.numeric(LG_X[2,]), as.character(LG_X[4,]), as.numeric(LG_X[5,]), 
                  pseudo.y.LG_X, pseudo.y.LG_X.right)
  
  #Rounding values above 0.5 to 0.500 to facilitate comparison with y:
  output[4,][which(output[4,] > 0.5)] <- 0.5
  output[5,][which(output[5,] > 0.5)] <- 0.5
  
  
  colnames(output) <- colnames(XO.phased.genot)
  rownames(output) <- c(paste("LG", LG, "_cM", sep=""), paste("LG_arm", sep=""), paste("y (Ho)", sep=""),
                        paste("LG", LG, "_Left", sep=""), paste("LG", LG, "_Right", sep=""))
  return(output)
}

#############
## TESTING ##
#############

x <- plot.pseudo.y(LG=1, plot="Y")


#PLOT LINKAGE GROUPS BY NAME, MUST HAVE SAME NAME AS INPUT FILE

pdf(file = "G:/Analysis/Mapping/AllHaps/Centromeres/plots/morten_edit_Bird13_108")

plot.pseudo.y(LG="1", plot="Y")
plot.pseudo.y(LG="2", plot="Y")
plot.pseudo.y(LG="3", plot="Y")
plot.pseudo.y(LG="4", plot="Y")
plot.pseudo.y(LG="5", plot="Y")
plot.pseudo.y(LG="6", plot="Y")
plot.pseudo.y(LG="7", plot="Y")
plot.pseudo.y(LG="8", plot="Y")
plot.pseudo.y(LG="9", plot="Y")
plot.pseudo.y(LG="10", plot="Y")
plot.pseudo.y(LG="11", plot="Y")
plot.pseudo.y(LG="12", plot="Y")
plot.pseudo.y(LG="13", plot="Y")
plot.pseudo.y(LG="14", plot="Y")
plot.pseudo.y(LG="15", plot="Y")
plot.pseudo.y(LG="16", plot="Y")
plot.pseudo.y(LG="17", plot="Y")
plot.pseudo.y(LG="18", plot="Y")
plot.pseudo.y(LG="19", plot="Y")
plot.pseudo.y(LG="20", plot="Y")
plot.pseudo.y(LG="21", plot="Y")
plot.pseudo.y(LG="22", plot="Y")
plot.pseudo.y(LG="23", plot="Y")
plot.pseudo.y(LG="24", plot="Y")
plot.pseudo.y(LG="25", plot="Y")
plot.pseudo.y(LG="26", plot="Y")

dev.off()

##############################################################################################
#######END HERE
##############################################################################################


# #######################################################
# ### Prepping plot for all LG's with RFm for kokanee ###
# #######################################################
# 
# #Define kokanee LG names:
# LGs <- as.character(c(1:26)) 
# LGs_Y_axis <- LGs[c(1,7,13,19,25)]
# LGs_X_axis <- LGs[c(25:30)]
# 
# xmax <- 150
# 
# LG.types <- as.matrix(read.csv(file = 'LG types.csv', sep = ',')) #Read map w. haploid genotypes
# 
# pdf(file="RFm_plots_kokaneePSV_map_inclCentro.pdf")
# par(mfrow=c(5,6),oma=c(4,3,3,2), mar=c(1.5,1,1,0))
# 
# #Figure 1 with LGs 1 & 3
# LGs <- c(1,3)
# LGs_Y_axis <- LGs[1]
# LGs_X_axis <- LGs
# #xmax <- 150
# xmax <- c(100,100,135)
# 
# pdf(file="RFm_plots_kokanee LGs 1 & 3_inclCentro.pdf")
# par(mfrow=c(1,2),oma=c(6,3,6,2), mar=c(1.5,1,1,0))


# ########################################################
# ### Prepping plot for all LG's with RFm for pink X05 ###
# ########################################################
# 
# #Define pink LG names:
# LGs <- as.character(c(1:26))
# LGs_Y_axis <- LGs[c(1,7,13,19,25)]
# LGs_X_axis <- c(21:26)
# 
# xmax <- 155
# 
# pdf(file="RFm_plots_PinkX05_map.pdf")
# par(mfrow=c(5,6),oma=c(4,3,3,2), mar=c(1.5,1,1,0))
# 
# ######################################################
# ### Prepping plot for all LG's with RFm for Barley ###
# ######################################################
# 
# #Define Barley LG names:
# LGs <- as.character(c(1:7))
# LGs_Y_axis <- LGs[c(1,3,5,7)]
# LGs_X_axis <- LGs[c(6:7)]
# 
# xmax <- 250
# 
# pdf(file="RFm_plots_Barley_map2.pdf")
# par(mfrow=c(4,2),oma=c(4,3,3,2), mar=c(3,2,2,1))
# 
# 
# ###############################################################
# ### Prepping plot for all LG's with RFm for mykiss (Miller) ###
# ###############################################################
# 
# #Define mykiss LG names:
# LGs <- as.character(c(1:29))
# LGs_Y_axis <- LGs[c(1,7,13,19,25)]
# LGs_X_axis <- c(24:29)
# 
# xmax <- 130
# 
# pdf(file="RFm_plots_mykiss_map.pdf")
# par(mfrow=c(5,6),oma=c(4,3,3,2), mar=c(1.5,1.5,1,0))
# 
# 
# ##################################################################
# ### Prepping plot for all LG's with RFm for Chinook (McKinney) ###
# ##################################################################
# 
# #Define mykiss LG names:
# LGs <- as.character(c(1:34))
# LGs_Y_axis <- LGs[c(1,7,13,19,25,31)]
# LGs_X_axis <- c(29:34)
# 
# xmax <- 140
# 
# pdf(file="RFm_plots_Chinook_KUWDip_F13_M7_Sire_MSTMapPhased_LGpos.pdf")
# par(mfrow=c(6,6),oma=c(4,3,3,2), mar=c(1.5,1.5,1,0))
# 


##################################################################
### Prepping plot for all LG's with RFm for pink ###
##################################################################

#Define mykiss LG names:
LGs <- as.character(c(1:26))
LGs_Y_axis <- LGs[c(1,7,13,19,25)]
LGs_X_axis <- c(29:34)

xmax <- 140

pdf(file="RFm_plots_pink_psv_wcentro.pdf")
par(mfrow=c(6,6),oma=c(4,3,3,2), mar=c(1.5,1.5,1,0))

#######################
### PLOTTING SCRIPT ###
#######################
par(mfrow=c(5,6),oma=c(4,3,3,2), mar=c(1.5,1,1,0))
pdf(file="RFm_plots_pink_psv_wcentro.pdf")

#######

RFm_matrix <- NULL
for (j in LGs){
  #Running RFm function
  pseudo.y.matrix <- plot.pseudo.y(LG=j, plot="N")
  RFm_matrix <- cbind(RFm_matrix, pseudo.y.matrix)
  
  #Ploting box, axes, labels etc.
  #plot(x=NULL, y=NULL, xlim=c(0,xmax[j]), axes=F,ann=F,
  plot(x=NULL, y=NULL, xlim=c(0,xmax), axes=F,ann=F,
       ylim=c(0,0.6))
  #axis(1, las=1, at=c(0,50,100,150,200,250), labels=rep("",6), col="gray80")
  lines(x=c(-10,260), y=c(0.5, 0.5), lty=2) #Maximum RFm value expected with compelte interference
  #lines(x=c(-10,260), y=c(0.45, 0.45), lty=2) #RFm threshold for defining centromeres according to y < 0.1
  box(col="gray80")

  if (j %in% LGs_Y_axis) {     #true when i is 1, 5, 9, 13 (%% is remainder when i divided by 4)
    axis(2, las=1, at=seq(0,1.4,0.2), col="gray80")
    mtext(text=expression(paste("RF"[m])) , side=2, line=2.5, cex=0.6, las=0)
  }

  if (j %in% LGs_X_axis){
    axis(1, las=1, at=seq(0,150,25),col="gray80")  #reduce the number of axis labels
    mtext(text="cM", side=1, line=2, cex=0.6)
  }

#   #Plotting centromeric region for kokanee map:
#   Centromeric_pos <- NULL #Making vector of all LG positions in centromere
#   
#   for (i in 1:ncol(pseudo.y.matrix)){ 
#     if (pseudo.y.matrix[2,i] == "c"){
#       Centromeric_pos <- cbind(Centromeric_pos, as.numeric(as.character(pseudo.y.matrix[1,i])))    
#     }  
#   }
#   
#   centro_start <- as.numeric(min(Centromeric_pos))
#   centro_end <- as.numeric(max(Centromeric_pos))
#   
#   if (centro_end - centro_start < 2){
#     centro_start = centro_start-1
#     centro_end = centro_end+1
#   }
#   
#   centro_x <- c(centro_start, centro_end, centro_end, centro_start)
#   centro_y <- c(-0.1, -0.1, 2, 2)
#   polygon(x=centro_x, y=centro_y, col="lightgrey", border= NA)  
#   
  #Plotting RFm
  #points(x=pseudo.y.matrix[1,], y=as.numeric(pseudo.y.matrix[3,]), pch=19, col="#00000033", cex=0.6) #Plotting y for kokanee
  points(x=pseudo.y.matrix[1,], y=as.numeric(pseudo.y.matrix[4,]), pch=19, col="#0000FF33", cex=0.6) #Plotting RFm from left
  points(x=pseudo.y.matrix[1,], y=as.numeric(pseudo.y.matrix[5,]), pch=19, col="#FF000033", cex=0.6) #Plotting RFm from right
  
  #Plotting LG title
  mtext(text=paste("LG", j), side=3, line=0.1, cex=0.6)
  #mtext(text=paste("LG", j, " (", LG.types[which(LG.types[,1]==j),2],")", sep=""), side=3, line=0.1, cex=0.6)
}  

dev.off()

write.csv(t(RFm_matrix), file="pink_RFm.csv")

####################################################
### Plotting LG's 1 & 3 with y & RFm for kokanee ###
####################################################

#Read in map with renamed LGs according to sockeye map & Gyno Ho values:
mapped_loci_PSV <- read.csv(file = 'PSV_map_nerka2.csv', sep = ',')

#Remove LGs 31 & 32:
#LGs_Keep <- c(1,3)
#mapped_loci_PSV <- subset(mapped_loci_PSV, mapped_loci_PSV$LG_nerka %in% LGs_Keep)

###PSV map:
#LGs <- as.character(LGs_Keep)
#LGs_Y_axis <- LGs[c(1)]
#LGs_X_axis <- LGs
#col <- as.character(mapped_loci_PSV$color)
#alpha <- as.character(mapped_loci_PSV$alpha)
#pairs <- read.csv(file = 'PSV_map_nerka2B_pairs.csv', sep = ',')

pdf(file="y_&_Hy_plots_LGs_1&5.pdf")
par(mfrow=c(1,2),oma=c(8,3,8,1), mar=c(1.5,1,1,0))

for (j in LGs){
  plot(x=NULL, y=NULL, xlim=c(0,145), axes=F,ann=F,
       ylim=c(0,1.3))
  mtext(text=paste("LG", j), side=3, line=0.5, cex=1)
  box(col="gray80")
  lines(x=c(-10,160), y=c(0.67, 0.67), lty=2) #Maximum y value expected with no interference
  lines(x=c(-10,160), y=c(1, 1), lty=2, col="grey")  #Y value expected on opposite arm of the one synced (i.e. random a's and b's)
  lines(x=c(-10,160), y=c(1.2, 1.2), lty=1, col="grey")  # Lower edge for LG box
  if (j %in% LGs_Y_axis) {     #true when i is 1 (%% is remainder when i divided by 4)
    axis(2, las=1, at=seq(0,1.2,0.2), col="gray80")
    mtext(text="y / Hy", side=2, line=2.5, cex=1, adj=0.45)
  }
  if (j %in% LGs_X_axis){
    axis(1, las=1, col="gray80")  #reduce the number of axis labels
    mtext(text="cM", side=1, line=2, cex=1)
  }
  
  for (i in 1:length(mapped_loci_PSV$LG_nerka)){  
    
    if (mapped_loci_PSV$LG_nerka[i]==j){
      #Plotting true y from gyno-dips
      points(mapped_loci_PSV$cM[i], as.vector(mapped_loci_PSV$Y_New_Filt10[i]), 
             type="p", pch=19, cex=0.6, col="#00000033")
  
      #Plotting map with NonPSVs and PSVs
      points(mapped_loci_PSV$cM[i],  1.28, type="p", pch=19, col=col[i], cex=0.6)
      
    }}
  
    #Plotting pseudo-y
    pseudo.y.matrix <- plot.pseudo.y(LG=j, plot="N")
    points(x=pseudo.y.matrix[1,], y=1-as.numeric(pseudo.y.matrix[4,])*2, pch=19, col="#0000FF33", cex=0.6)
    points(x=pseudo.y.matrix[1,], y=1-as.numeric(pseudo.y.matrix[5,])*2, pch=19, col="#FF000033", cex=0.6)
  
  
}

dev.off()


###############################################################
###  Figure 1 with phase plot (a) and RFm illustration (b)  ###
### This is a draft figure that was polished in Illustrator ###
###############################################################

##Read in map with phased genotypes:
Map <- as.matrix(read.csv(file = 'Map_phased_genotypes_haploids_for_R.csv', sep = ',')) #Haploid map

LG <- 1

##Making subset matrix for LG_X
LG_X <- NULL
LG_X_arms <- NULL
for (i in 1:length(Map[1,])){
  if (Map[1,i] == LG){  
    LG_X <- cbind(LG_X, Map[,i])
    LG_X_arms <- unique(LG_X[4,])
  }
}


##Make new matrix where a & b are changed to 0 and 1:
XO.vec.kw <- matrix(nrow=nrow(Map), ncol=ncol(LG_X))  #Make Cross-over (XO) vector (vec) using kernel-window (kw) for focal LG
colnames(XO.vec.kw) <- colnames(LG_X)

for (j in 6:nrow(XO.vec.kw)){
  for (i in 1:ncol(XO.vec.kw)){
    if (as.character(LG_X[j,i]) == "a"){
      XO.vec.kw[j,i] <- 0  
    }
    if (as.character(LG_X[j,i]) == "b"){
      XO.vec.kw[j,i] <- 1  
    }  
  }
}


#############################################################
#### Plot haplotype (0 or 1) across LG for each offspring ###
#############################################################

#Defining centromeric region:
Centromeric_pos <- NULL
for (i in 1:ncol(LG_X)){ #Making vector of all LG positions in centromere
  if (LG_X[4,i] == "c"){
    Centromeric_pos <- cbind(Centromeric_pos, LG_X[2,i])    
  }  
}
centro_start <- as.numeric(min(Centromeric_pos))
centro_end <- as.numeric(max(Centromeric_pos))
if (centro_end-centro_start < 1){ #If centromeric region is less than 1 cM wide it is extended 1 cM in both directions
  centro_start <- centro_start-1
  centro_end <- centro_end+1
}
centro_x <- c(centro_start, centro_end, centro_end, centro_start)
centro_y <- c(-0.5, -0.5, 1.9, 1.9)


###################
### Making plot ###
###################

#Plotting Figure 1a:
cM <- LG_X[2,]
offspring <- 6
title <- paste("LG",LG,"offspring", offspring)


pdf(file="Figure 1b.pdf")
par(mfrow=c(2,1), mar=c(3,2,3,2), oma=c(4,6,4,6))

plot(NULL, xlab="", ylab="", xlim=c(-5,80), ylim=c(-0.5,1.5), xaxs="i",yaxs="i", xaxt="n", yaxt="n")
#Plotting y-axes
axis(side=2, at=c(0,1), labels=c("a", "b"), las=1)
mtext(text="Phase", side=2, line=2, adj=0.5, cex=0.8)  
  
#Plotting x-axes 
axis(side=1, at=seq(0,120,20), labels=seq(0,120,20))
mtext(text="cM", side=1, line=2, cex=0.8)  
  
#Plotting phases
polygon(x=centro_x, y=centro_y, col="lightgrey", border= NA)
points(x=cM, y=XO.vec.kw[offspring,], type="p", pch=19, col=c(rep("red", 47), rep("blue", 45)))  
mtext(side=3, line=0.2, text=title, cex=0.8)
box(col="black")


#Plotting Figure 1b:
plot(NULL, xlab="", ylab="", xlim=c(-10,100), ylim=c(-30,220), xaxs="i",yaxs="i", xaxt="n", yaxt="n")
#1st bar
polygon(x=c(10,90,90,10), y=c(180,180,190,190), col="red", border= NA)
#2nd bar
polygon(x=c(10,15,15,10), y=c(160,160,170,170), col="red", border= NA)
polygon(x=c(15,90,90,15), y=c(160,160,170,170), col="blue", border= NA)
#3rd bar
polygon(x=c(10,90,90,10), y=c(140,140,150,150), col="red", border= NA)
#4th bar
polygon(x=c(10,30,30,10), y=c(120,120,130,130), col="red", border= NA)
polygon(x=c(30,90,90,30), y=c(120,120,130,130), col="blue", border= NA)
#5th bar
polygon(x=c(10,90,90,10), y=c(110,110,100,100), col="red", border= NA)
#6th bar
polygon(x=c(10,45,45,10), y=c(80,80,90,90), col="red", border= NA)
polygon(x=c(45,90,90,45), y=c(80,80,90,90), col="blue", border= NA)
#7th bar
polygon(x=c(10,90,90,10), y=c(70,70,60,60), col="red", border= NA)
#8th bar
polygon(x=c(10,60,60,10), y=c(40,40,50,50), col="red", border= NA)
polygon(x=c(60,90,90,60), y=c(40,40,50,50), col="blue", border= NA)
#9th bar
polygon(x=c(10,90,90,10), y=c(30,30,20,20), col="red", border= NA)
#10th bar
polygon(x=c(10,75,75,10), y=c(0,0,10,10), col="red", border= NA)
polygon(x=c(75,90,90,75), y=c(0,0,10,10), col="blue", border= NA)
#Centromere:
#polygon(x=c(90,95,95,90), y=c(0,0,190,190), col="lightgrey", border= NA)

#Plotting marker position transects:
points(x=c(12.5, 12.5), y=c(-5,195), type="l", lty=2)
points(x=c(22.5, 22.5), y=c(-5,195), type="l", lty=2)
points(x=c(37.5, 37.5), y=c(-5,195), type="l", lty=2)
points(x=c(52.5, 52.5), y=c(-5,195), type="l", lty=2)
points(x=c(67.5, 67.5), y=c(-5,195), type="l", lty=2)
points(x=c(82.5, 82.5), y=c(-5,195), type="l", lty=2)

#Plotting marker names above phase data:
text(x=c(12.5, 22.5, 37.5, 52.5, 67.5, 82.5), y=(200), labels=c(expression(paste("m"[1])), 
                                                                expression(paste("m"[2])), 
                                                                expression(paste("m"[3])),
                                                                expression(paste("m"[4])), 
                                                                expression(paste("m"[5])), 
                                                                expression(paste("m"[6]))))

#Plotting marker names above phase data:
text(x=-5, y=200, labels="marker")
text(x=-5, y=-12, labels=expression(paste("RF"[m])))
text(x=c(12.5, 22.5, 37.5, 52.5, 67.5, 82.5), y=(-12), labels=as.character(seq(0.0,0.5,0.1)))


dev.off()



#######################
### END FOR PAPER 3 ###
#######################

#####################################################
### Plotting y and RFm distribution barplots ###
#####################################################

###Making vectors with y and pseudo.y values for all LGs excluding 18B
LGs <- as.character(c(1:8,"9A_(X2)", "9B_(X1)", 10:17, "18A", 19:28)) #Ignoring 18B

y.All.LG <- NULL
y.All.LG.arms <- NULL
pseudo.y.All.LG <- NULL
pseudo.y.All.LG.arms <- NULL

for (j in LGs){
  LG_X <- plot.pseudo.y(LG=j, plot="N")
  pseudo.y.LG_x.arms <- NULL
  pseudo.y.LG_x.cen <- NULL
  paste("LG",j, ", No. markers:", ncol(LG_X), sep="")  
  for (i in 1:ncol(LG_X)){
    y.All.LG <- cbind(y.All.LG, LG_X[3,i])
    if (LG_X[2,i] == "a1"){
      pseudo.y.LG_x.arms <- cbind(pseudo.y.LG_x.arms, (1-2*as.numeric(LG_X[4,i])))
      y.All.LG.arms <- cbind(y.All.LG.arms, LG_X[3,i])
    }
    if (LG_X[2,i] == "a2"){
      pseudo.y.LG_x.arms <- cbind(pseudo.y.LG_x.arms, (1-2*as.numeric(LG_X[5,i])))
      y.All.LG.arms <- cbind(y.All.LG.arms, LG_X[3,i])
    }
    if (LG_X[2,i] == "c"){
      pseudo.y.LG_x.cen <- cbind(pseudo.y.LG_x.cen, (1-2*as.numeric(LG_X[4,i]))) #Centromeres are scored as the same phase as a1 to cover whole centromeres in acrocentrics as well
    }
  }
  pseudo.y.All.LG.arms <- cbind(pseudo.y.All.LG.arms, pseudo.y.LG_x.arms)
  pseudo.y.All.LG <- cbind(pseudo.y.All.LG, pseudo.y.LG_x.arms, pseudo.y.LG_x.cen)
}


##Making histogram counts to be plotted:

#Individual histograms:
bins <- seq(0,1,0.1)

H.y.All.LG.arms <- hist(as.numeric(y.All.LG.arms), breaks=bins)
H.y.All.LG <- hist(as.numeric(y.All.LG), breaks=bins)
H.pseudo.y.arms <- hist(pseudo.y.All.LG.arms, breaks=bins)
H.pseudo.y.All.LG <- hist(pseudo.y.All.LG, breaks=bins)

# #Number of loci for counts with and without centromeres used to estimate proportions
# no.loci.y.All.LG.arms <- length(y.All.LG.arms)
# no.loci.y.All.LG <- length(y.All.LG)
# no.loci.pseudo.y.All.LG.arms <- length(pseudo.y.All.LG.arms)
# no.loci.pseudo.y.All.LG <- length(pseudo.y.All.LG)

sum(H.y.All.LG.arms$counts)

##Making barplots
#Setting universal parameters:
las = 1
width = 0.1
cex=0.7
xlab <- seq(0.1,1,0.1)
xlim=NULL #c(0,1)  
xlab1="y"
xlab2="Hy"
cex.names=0.7
cex.axis=0.7
ylim1=c(0,0.50)  
ylim2=c(0,0.50)  

#Individual Barplots
dev.off()
pdf(file="y_and_Hy_distributions.pdf")
par(mfrow=c(2,1),oma=c(3,8,2,8))

# y #
par(mar=c(0.5,4,2,0))
barplot(H.y.All.LG$counts/sum(H.y.All.LG$counts), col="grey", names.arg=NULL, xlab=xlab1,
        ylab="Proportion of loci", ylim=ylim1, las=las, width=width, xlim=xlim,
        cex.names=cex.names, cex.axis=cex.axis)

box(col="gray80")
mtext(text=paste("y", "\n", "Non-duplicated loci"), side=3, line=-2, cex=cex, las=1)

# Hy #
par(mar=c(2,4,0.5,0))
barplot(H.pseudo.y.All.LG$counts/sum(H.pseudo.y.All.LG$counts), col="grey", names.arg=xlab,
        xlab=xlab2, ylab="Proportion of loci", ylim=ylim2, las=las, width=width, xlim=xlim, 
        cex.names=cex.names, cex.axis=cex.axis)

box(col="gray80")
mtext(text=paste("Hy", "\n", "All loci"), side=3, line=-2, cex=cex, las=1)
mtext(text="y / Hy", side=1, line=1.7, cex=cex, las=1)

dev.off()





### OLD PLOTS FOR EXCL. CENTROMERES ###
par(mar=c(2,4,2,0))
barplot(H.y.All.LG.arms$counts/sum(H.y.All.LG.arms$counts), col="pink", names.arg=NULL,
        ylab="Proportion of loci", xlab=xlab1, ylim=ylim1, las=las, width=width, xlim=xlim,
        cex.axis=cex.axis)
box(col="gray80")
mtext(text="y \n (excl. centromeres)", side=3, line=-2.2, cex=cex, las=1)


par(mar=c(4,4,0,0))
barplot(H.pseudo.y.arms$counts/sum(H.pseudo.y.arms$counts), col="lightblue", names.arg=xlab,
        ylab="Proportion of loci", xlab=xlab2, ylim=ylim2, las=las, width=width, xlim=xlim,
        cex.names=cex.names, cex.axis=cex.axis)
box(col="gray80")
mtext(text="haplo y \n (excl. centromeres)", side=3, line=-2.2, cex=cex, las=1)