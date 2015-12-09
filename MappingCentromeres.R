### Identify centromere location with phased haploid data
###    Written by Morten Limborg and edited by Garrett McKinney 
### Carolyn Tarpey | December  2015 
### ---------------------------------------


setwd("G:/Analysis/Mapping/AllHaps/Centromeres" )

######################################################################
#######Place centromeres and count y distribution with haploids#######
######################################################################

rm(list=ls())

##Read in map with phased genotypes:
Map <- as.matrix(read.csv(file = 'centmapping_01.txt', sep = '\t')) #Read map with imputated genotypes for missing values at first and last markers on each LG
is.matrix(Map)
rownames(Map) <- Map[,1]                #Change column 1 into rownames
Map <- Map[,2:ncol(Map)] #Should go in function                
RowNum<-dim(Map)[1]
FamSize<-RowNum-5

plot.pseudo.y <- function(LG, plot){
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
  
  #Make vector for pseudo y valyes for LG_X from right
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
   # output <-LG_X[2,]
    #dim(pseudo.y.LG_X*2)
    #dim(as.numeric(LG_X[2,]))
    Xvalues<-as.numeric(LG_X[2,])
    Yvalues<-as.numeric(pseudo.y.LG_X*2)
    Yvalues_right<-as.numeric(pseudo.y.LG_X.right*2)
    loess_fit_left <- loess(Yvalues ~ Xvalues)
    loess_fit_right <- loess(Yvalues_right ~ Xvalues)
    #output<-Xvalues
    #output<-Yvalues
    #output<-YvaluesRight
    #output<-predict(loess_fit_left)
#  plot(x=as.numeric(LG_X[2,]), y=pseudo.y.LG_X*2, pch=19, col="red", xlab="cM", ylab="Pseudo y", 
#       xlim=c(0,1.05*max(x)), ylim=c(0,1.2), las=1)
  plot(x=as.numeric(LG_X[2,]), y=pseudo.y.LG_X*2, pch=19, col="red", xlab="cM", ylab="Pseudo y", main=LG, 
     xlim=c(0,160), ylim=c(0,1.2), las=1)
       #xlim=c(0,160), ylim=c(0,1.2), las=1)
  lines(x=c(-10,160), y=c(0.67, 0.67), lty=2)  
  lines(x=c(-10,160), y=c(1, 1), lty=2, col="grey")  
  points(x=as.numeric(LG_X[2,]), y=pseudo.y.LG_X.right*2, pch=19, col="blue")
  points(x=as.numeric(LG_X[2,]), y=LG_X[5,], pch=19, col="black")
  lines(Xvalues,predict(loess_fit_left),col="black",lwd=2)
  lines(Xvalues,predict(loess_fit_right),col="black",lwd=2)

  }

  
  
  #Making output matrix with pseudo y values
  output <- rbind(as.numeric(LG_X[2,]), as.character(LG_X[4,]), as.numeric(LG_X[5,]), 
                  pseudo.y.LG_X, pseudo.y.LG_X.right)
  colnames(output) <- colnames(XO.phased.genot)
  rownames(output) <- c(paste("LG", LG, "_cM", sep=""), paste("LG_arm", sep=""), paste("y (Ho)", sep=""),
                        paste("LG", LG, "_Left", sep=""), paste("LG", LG, "_Right", sep=""))
  
  
  return(output)
}


#PLOT LINKAGE GROUPS BY NAME, MUST HAVE SAME NAME AS INPUT FILE
plot.pseudo.y(LG="1", plot="Y")
plot.pseudo.y(LG="19", plot="Y")

##############################################################################################
#######END HERE
##############################################################################################
