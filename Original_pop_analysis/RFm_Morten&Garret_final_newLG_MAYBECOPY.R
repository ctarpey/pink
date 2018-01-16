### Identify centromere location with phased haploid data
###  uses a Loess regression to report the centromere distances
###    Written by Morten Limborg and edited by Garrett McKinney 
###    This version rounds RFm down to 0.5
###  This is to place the centormeres on the 2 new LGs that garrett re-ordered
### Carolyn Tarpey | October 2016
### ---------------------------------------

##Nov 2016 changed the colors to green, purple and yellow
##so as not to interfere with the red/blue color scheme of the duplicates

###########################################################################
####### There is an issue with getting LG 1 to plot in most families#######
###########################################################################


#MAKE THESE CHANGES TO PLOT LG1. The rest of the LGs dont need these modifications: 

# for Hood 11 _05 it needs the right cent flank choose 2
# for Hood 11_01 it needs both the right and left cent flank choose 2
# for Bird13_108 it nees only the the right cent flank choose 2
# for bird13_110 it seems to work with the right and left choose left off but failed for a lot of the other lgs. 
# 
# #
# LeftCentValues<-mapply("-",CentLimit,leftValues)
# LeftCentFlank<-which(diff(sign(LeftCentValues))!=0)
# 
# # RW 1-12                             #edit to make it work with LG 1 (comment these three lines out for the other LGs)
# print(LeftCentFlank)
# LeftCentFlank <- LeftCentFlank[[2]]
# print(LeftCentFlank)
# 
# LeftCentStart<-LeftCentFlank+1
# 
# 
# Or else you will get this error: 
# 
# [1] 175
# Error in RightCentFlank[[2]] : subscript out of bounds
# In addition: There were 50 or more warnings (use warnings() to see the first 50)
# 
#I ended up plotting LG 1 separate for each family and merging the PDFs to the other LGs. 


setwd("Y:/WORK/TARPEY/Pink Linkage Map/Centromeres/NewLGs" )

######################################################################
#######Place centromeres and count y distribution with haploids#######
######################################################################

rm(list=ls())

plot.pseudo.y <- function(LG, plot){
  
  ##Read in map with phased genotypes:
  #Map <- as.matrix(read.csv(file = 'Cents_05.txt', sep = '\t')) #Read map with imputated genotypes for missing values at first and last markers on each LG
  Map <- as.matrix(read.csv(file = 'Cents_01.txt', sep = '\t')) #Read map with imputated genotypes for missing values at first and last markers on each LG
  #Map <- as.matrix(read.csv(file = 'Cents_110.txt', sep = '\t')) #Read map with imputated genotypes for missing values at first and last markers on each LG
  #Map <- as.matrix(read.csv(file = 'Cents_108.txt', sep = '\t')) #Read map with imputated genotypes for missing values at first and last markers on each LG
  
  #Map <- as.matrix(read.csv(file = 'New_cents_05.txt', sep = '\t')) #Read map with imputated genotypes for missing values at first and last markers on each LG
  #Map <- as.matrix(read.csv(file = 'New_cents_01.txt', sep = '\t')) #Read map with imputated genotypes for missing values at first and last markers on each LG
  #Map <- as.matrix(read.csv(file = 'New_cents_110.txt', sep = '\t')) #Read map with imputated genotypes for missing values at first and last markers on each LG
  #Map <- as.matrix(read.csv(file = 'New_cents_108.txt', sep = '\t')) #Read map with imputated genotypes for missing values at first and last markers on each LG
  
  is.matrix(Map)
  
  rownames(Map) <- Map[,1]                #Change column 1 into rownames
  Map <- Map[,2:ncol(Map)] 
  
  RowNum<-dim(Map)[1]
  FamSize<-RowNum-5
  
  
  #Making subset matrix for LG_X:
  LG <- LG #Could be looped over in the end  
  
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
  XO.phased.genot <- XO.vec.kw[6:RowNum,]
  
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
  
  for (j in (no.zero+1):nrow(XO.phased.genot)){  #Loop over individuals having the genotype 1 at the first marker
    for (i in 1:ncol(XO.phased.genot)){ #For individual j: Loop over all markers in LG_X
      converted.genot <- cbind(converted.genot, abs(as.numeric(XO.phased.genot[j,i])-1))
      
    }
    Matrix.convert <- rbind(Matrix.convert, converted.genot)
    converted.genot <- NULL
  }
  
  #From right:
  Matrix.convert.right <- NULL
  converted.genot.right <- NULL #Temporary vector for each individual to be used within loop
  
  for (j in (no.zero.right+1):nrow(XO.phased.genot.right)){  #Loop over individuals having the genotype 1 at the last marker
    for (i in 1:ncol(XO.phased.genot.right)){ #For individual j: Loop over all markers in LG_X
      converted.genot.right <- cbind(converted.genot.right, abs(as.numeric(XO.phased.genot.right[j,i])-1))
      
    }
    Matrix.convert.right <- rbind(Matrix.convert.right, converted.genot.right)
    converted.genot.right <- NULL
  }
  
  
  ##Merge synced phased offspring into final matrix for counting cross-over frequencies
  #From left:
  phased.matrix <- rbind(XO.phased.genot[1:no.zero,], Matrix.convert)
  rownames(phased.matrix) <- rownames(XO.phased.genot)
  
  #Make vector for pseudo y values for LG_X from right
  pseudo.y.LG_X <- NULL
  
  #######################################################################################
  #Test Modifications Below
  #######################################################################################
  for (i in 1:ncol(phased.matrix)){
    #Pseudo.Y is calculated by effectively taking the average based on family size, but excluding
    #NAs from the sum in the numerator
    pseudo.y <- mean(as.numeric(phased.matrix[,i]),na.rm=TRUE)
    #pseudo.y <- sum(as.numeric(phased.matrix[,i]), na.rm = TRUE)/FamSize
    pseudo.y.LG_X <- cbind(pseudo.y.LG_X, pseudo.y)
  }
  
  #From right:
  phased.matrix.right <- rbind(XO.phased.genot.right[1:no.zero.right,], Matrix.convert.right)
  rownames(phased.matrix.right) <- rownames(XO.phased.genot.right)
  
  #Make vector for pseudo y valyes for LG_X from right
  pseudo.y.LG_X.right <- NULL
  
  for (i in 1:ncol(phased.matrix.right)){
    #Make same modification to pseudo.y calculation below by taking mean instead
    pseudo.y.right <- mean(as.numeric(phased.matrix.right[,i]),na.rm=TRUE)
    #pseudo.y.right <- sum(as.numeric(phased.matrix.right[,i]), na.rm = TRUE)/FamSize
    pseudo.y.LG_X.right <- cbind(pseudo.y.LG_X.right, pseudo.y.right)
  }
  
  ##########################################################################
  #Try adding trendline to plotting below, choose intersect of line as max 
  #value for recombination threshold
  ##########################################################################
  
  
  #Plot pseudo-y values from left and right:
  if (plot == "Y"){
    Xvalues<-as.numeric(LG_X[2,])
    Yvalues<-as.numeric(pseudo.y.LG_X)
    Yvalues_right<-as.numeric(pseudo.y.LG_X.right)
    loess_fit_left <- loess(Yvalues ~ Xvalues)
    loess_fit_right <- loess(Yvalues_right ~ Xvalues)
    
    ########################################################################
    #add intersect to plot
    ########################################################################
    
    leftValues<-predict(loess_fit_left)
    #output<-leftValues
    rightValues<-predict(loess_fit_right)
    subTest<-mapply("-",leftValues,rightValues)
    which(diff(sign(subTest))!=0)
    Flank1<-which(diff(sign(subTest))!=0)
    Flank2<-Flank1+1
    #
    leftYflank1<-leftValues[[Flank1]]
    leftYflank2<-leftValues[[Flank2]]
    rightYflank1<-rightValues[[Flank1]]
    rightYflank2<-rightValues[[Flank2]]
    flank1X<-Xvalues[[Flank1]]
    flank2X<-Xvalues[[Flank2]]
    leftSlope<-(leftYflank2-leftYflank1)/(flank2X-flank1X)
    rightSlope<-(rightYflank2-rightYflank1)/(flank2X-flank1X)
    leftIntercept<-leftYflank1-(leftSlope*flank1X)
    rightIntercept<-rightYflank1-(rightSlope*flank1X)
    Xintersect<-(rightIntercept-leftIntercept)/(leftSlope-rightSlope)
    Yintersect<-leftSlope*Xintersect+leftIntercept
    Max<-Yintersect
    AcroMax<-0.1
    CentLimit<-Max*0.9
    #
    LeftCentValues<-mapply("-",CentLimit,leftValues)
    LeftCentFlank<-which(diff(sign(LeftCentValues))!=0)
    
    # RW 1-12                             #edit to make it work with LG 1
    #print(LeftCentFlank)
    #LeftCentFlank <- LeftCentFlank[[2]]
    #print(LeftCentFlank)
    
    LeftCentStart<-LeftCentFlank+1
    LeftCentLimit<-Xvalues[[LeftCentStart]]
    RightCentValues<-mapply("-",CentLimit,rightValues)
    RightCentFlank<-which(diff(sign(RightCentValues))!=0)
    
    # RW 1-12                             #edit to make it work with LG 1
    #print(RightCentFlank)
    #RightCentFlank <- RightCentFlank[[2]]
    #print(RightCentFlank)
  
    RightCentLimit<-Xvalues[[RightCentFlank]]
#     LeftAcroValues<-mapply("-",AcroMax,leftValues)
#     LeftAcroFlank<-which(diff(sign(LeftAcroValues))!=0)
#     LeftAcroLimit<-Xvalues[[LeftAcroFlank]]  #May need to adjust this value, double check
#     RightAcroValues<-mapply("-",AcroMax,rightValues)
#     RightAcroFlank<-which(diff(sign(RightAcroValues))!=0)
#     RightAcroStart<-RightAcroFlank+1  #This will start at Y values<0.1 but not including 0.1
#     RightAcroLimit<-Xvalues[[RightAcroStart]]
    #
    plot(NULL, pch=19, col="#00b05f", xlab="cM", ylab="RF (phase-shift)", xlim=c(0,160), ylim=c(0,1), las=1, main=paste("LG",LG))
    lines(x=c(-10,160), y=c(0.5, 0.5), lty=2)  
    lines(x=c(-10,160), y=c(1, 1), lty=2, col="grey")  
    points(x=as.numeric(LG_X[2,]), y=pseudo.y.LG_X, pch=19, col="#00b05f")
    points(x=as.numeric(LG_X[2,]), y=pseudo.y.LG_X.right, pch=19, col="#6965dd")
    lines(Xvalues,predict(loess_fit_left),col="black",lwd=2)
    lines(Xvalues,predict(loess_fit_right),col="black",lwd=2)
    points(x=Xintersect,y=Yintersect,pch= 24, col="orange", bg="orange")
    
    text(140,0.4,paste("X: ",round(Xintersect,2),"\nY: ",round(Yintersect,2)),pos=4)
    #text(140,0.2,paste("Limit: ", round(CentLimit,3)),pos=4)
    text(140,0.2,paste("Lx: ", round(LeftCentLimit,2)),pos=4)
    text(140,0.15,paste("Rx: ", round(RightCentLimit,2)),pos=4)
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

#PLOT LINKAGE GROUPS BY NAME, MUST HAVE SAME NAME AS INPUT FILE


pdf(file = "Y:/WORK/TARPEY/Pink Linkage Map/Centromeres/NewLGs/Plots/110_1.pdf")

plot.pseudo.y(LG="1", plot="Y")

dev.off()

################################## Change the settings

pdf(file = "Y:/WORK/TARPEY/Pink Linkage Map/Centromeres/NewLGs/Plots/110_2-16.pdf")

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

dev.off()

pdf(file = "Y:/WORK/TARPEY/Pink Linkage Map/Centromeres/NewLGs/Plots/110_18-23.pdf")

plot.pseudo.y(LG="18", plot="Y")
plot.pseudo.y(LG="19", plot="Y")
plot.pseudo.y(LG="20", plot="Y")
plot.pseudo.y(LG="21", plot="Y")
plot.pseudo.y(LG="22", plot="Y")
plot.pseudo.y(LG="23", plot="Y")
dev.off()

pdf(file = "Y:/WORK/TARPEY/Pink Linkage Map/Centromeres/NewLGs/Plots/110_25-26.pdf")
plot.pseudo.y(LG="25", plot="Y")
plot.pseudo.y(LG="26", plot="Y")

dev.off()

################################### change to new input file

pdf(file = "Y:/WORK/TARPEY/Pink Linkage Map/Centromeres/NewLGs/Plots/110_24.1.pdf")
plot.pseudo.y(LG="24", plot="Y")

dev.off()

pdf(file = "Y:/WORK/TARPEY/Pink Linkage Map/Centromeres/NewLGs/Plots/110_17.pdf")

plot.pseudo.y(LG="17", plot="Y")

dev.off()



##############################################################################################
#######END HERE
##############################################################################################


