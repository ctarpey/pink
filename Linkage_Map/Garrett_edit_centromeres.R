######################################################################
#######Place centromeres and count y distribution with haploids#######
######################################################################

setwd("G:/Analysis/Mapping/AllHaps/Centromeres" )

rm(list=ls())

##Read in map with phased genotypes:

Map <- as.matrix(read.csv(file = 'Cents_05.txt', sep = '\t')) #Read map with imputated genotypes for missing values at first and last markers on each LG
#Map <- as.matrix(read.csv(file = 'Cents_01.txt', sep = '\t')) #Read map with imputated genotypes for missing values at first and last markers on each LG
#Map <- as.matrix(read.csv(file = 'Cents_110.txt', sep = '\t')) #Read map with imputated genotypes for missing values at first and last markers on each LG
#Map <- as.matrix(read.csv(file = 'Cents_108.txt', sep = '\t')) #Read map with imputated genotypes for missing values at first and last markers on each LG


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
    #output<-loess_fit_left
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
    print(LeftCentFlank)
    LeftCentFlank <- LeftCentFlank[[1]]
    print(LeftCentFlank)
    
    LeftCentStart<-LeftCentFlank+1
    LeftCentLimit<-Xvalues[[LeftCentStart]]
    RightCentValues<-mapply("-",CentLimit,rightValues)
    RightCentFlank<-which(diff(sign(RightCentValues))!=0)
    
    # RW 1-12                             #edit to make it work with LG 1
    #print(RightCentFlank)
   # RightCentFlank <- RightCentFlank[[2]]
    #print(RightCentFlank)
  
    
    RightCentLimit<-Xvalues[[RightCentFlank]]
    LeftAcroValues<-mapply("-",AcroMax,leftValues)
    LeftAcroFlank<-which(diff(sign(LeftAcroValues))!=0)
    LeftAcroLimit<-Xvalues[[LeftAcroFlank]]  #May need to adjust this value, double check
    RightAcroValues<-mapply("-",AcroMax,rightValues)
    RightAcroFlank<-which(diff(sign(RightAcroValues))!=0)
    RightAcroStart<-RightAcroFlank+1  #This will start at Y values<0.1 but not including 0.1
    RightAcroLimit<-Xvalues[[RightAcroStart]]
    
    
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
  points(x=Xintersect,y=Yintersect,col="green")
  #
  text(140,0.4,paste("X: ",round(Xintersect,2),"\nY: ",round(Yintersect,2)),pos=4)
  #text(140,0.2,paste("Limit: ", round(CentLimit,3)),pos=4)
  text(140,0.2,paste("Lx: ", round(LeftCentLimit,2)),pos=4)
  text(140,0.15,paste("Rx: ", round(RightCentLimit,2)),pos=4)
  text(140,0.1,paste("ALx: ", round(LeftAcroLimit,2)),pos=4)
  text(140,0.05,paste("ARx: ", round(RightAcroLimit,2)),pos=4)

  }

  
  
  #Making output matrix with pseudo y values
  output <- rbind(as.numeric(LG_X[2,]), as.character(LG_X[4,]), as.numeric(LG_X[5,]), 
                  pseudo.y.LG_X, pseudo.y.LG_X.right)
  colnames(output) <- colnames(XO.phased.genot)
  rownames(output) <- c(paste("LG", LG, "_cM", sep=""), paste("LG_arm", sep=""), paste("y (Ho)", sep=""),
                        paste("LG", LG, "_Left", sep=""), paste("LG", LG, "_Right", sep=""))
  
  
  #return(output)
}

pdf(file = "G:/Analysis/Mapping/AllHaps/Centromeres/plots/Hood11_05_gar")

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



locator()


plot.pseudo.y(LG="lg81", plot="Y")
LG107 <- plot.pseudo.y(LG="lg107", plot="N")
LG107
write.table(LG107,file='C:/Users/Garrett_2/Documents/Research/Marblemount/Centromere_Mapping/LG107.txt',sep="\t",eol="\n",row.names=TRUE,col.names=TRUE,qmethod = c("escape", "double"))

Map[4,]

Map

#Use Par to set number of graphs per image, look into using lattice instead....
par(mfrow=c(2,2))

#KUWHap14LG results, along with mergeFile1, mergeFile2, mergeFile3
for (y in 74:77){
  prefix<-c("lg")
  currentLG<-paste(prefix,y,sep="")
  print(currentLG)
  plot.pseudo.y(LG=currentLG, plot="Y")
}

#Hap10LG results
Hap10LG<-c("lg14","lg10","lg23","lg32","lg13","lg2","lg4","lg17","lg31","lg5","lg18","lg16","lg8","lg36","lg21","lg9","lg15","lg6","lg28","lg33","lg1","lg40","lg35","lg24","lg0","lg20","lg22","lg30","lg44","lg39","lg11","lg38","lg27","lg12","lg41","lg25","lg37","lg7","lg3","lg25")
for (y in Hap10LG){
    print (y)
    plot.pseudo.y(LG=y, plot="Y")
}


########################################################
###Plotting All LG's with Gyno y AND Haploid pseudo y###
########################################################
#Read in map with renamed LGs according to sockeye map & Gyno Ho values:
mapped_loci_PSV <- read.csv(file = 'PSV_map_nerka2.csv', sep = ',')

#Remove LGs 31 & 32:
LGs_Keep <- c(0:30)
mapped_loci_PSV <- subset(mapped_loci_PSV, mapped_loci_PSV$LG %in% LGs_Keep)

###PSV map:
LGs <- as.character(c(1:8,"9A_(X1)", "9B_(X2)", 10:17, "18A", "18B", 19:28)) 
LGs_Y_axis <- LGs[c(1,7,13,19,25)]
LGs_X_axis <- LGs[c(25:30)]
col <- as.character(mapped_loci_PSV$color)
#alpha <- as.character(mapped_loci_PSV$alpha)
#pairs <- read.csv(file = 'PSV_map_nerka2B_pairs.csv', sep = ',')

pdf(file="Y_plots_PSV_map_gynos_and_haplos_large.pdf")
par(mfrow=c(5,6),oma=c(4,4,3,1), mar=c(1.5,1,1,0))
for (j in LGs){
  plot(x=NULL, y=NULL, xlim=c(0,145), axes=F,ann=F,
       ylim=c(0,1.5))
  mtext(text=paste("LG", j), side=3, line=0.5, cex=0.6)
  box(col="gray80")
  lines(x=c(-10,160), y=c(0.67, 0.67), lty=2) #Maximum Y value expected with no interference
  lines(x=c(-10,160), y=c(1, 1), lty=2, col="grey")  #Y value expected on opposite arm of the one synced (i.e. random a's and b's)
  lines(x=c(-10,160), y=c(1.44, 1.44), lty=1, col="grey")  #Y value expected on opposite arm of the one synced (i.e. random a's and b's)
#  if (j %in% LGs_Y_axis) {     #true when i is 1, 5, 9, 13 (%% is remainder when i divided by 4)
    axis(2, las=1, col="gray80")
    mtext(text="Y", side=2, line=3, cex=0.6, las=2)
#  }
#  if (j %in% LGs_X_axis){
    axis(1, las=1, col="gray80")  #reduce the number of axis labels
    mtext(text="cM", side=1, line=2, cex=0.6)
#  }
  
  for (i in 1:length(mapped_loci_PSV$LG_nerka)){  
    
    if (mapped_loci_PSV$LG_nerka[i]==j){
      #Plotting true y from gyno-dips
      points(mapped_loci_PSV$cM[i], as.vector(mapped_loci_PSV$Ho_Gyn_Filt10.[i]), 
             type="p", pch=19, cex=0.6)
      #Plotting map with NonPSVs and PSVs
      points(mapped_loci_PSV$cM[i],  1.5, type="p", pch=19,  
             col=col[i], cex=0.6)
      
    }}

  #Plotting pseudo-y
  pseudo.y.matrix <- plot.pseudo.y(LG=j, plot="N")
  points(x=pseudo.y.matrix[1,], y=pseudo.y.matrix[2,]*2, pch=19, col="blue", cex=0.6)
  points(x=pseudo.y.matrix[1,], y=pseudo.y.matrix[3,]*2, pch=19, col="red", cex=0.6)
  

}

dev.off()


###########################################################
###Plotting pseudo LG's with Gyno y AND Haploid pseudo y###
###########################################################

###Making vectors with y and pseudo.y values for all LGs excluding 18B
LGs <- as.character(c(1:8,"9A_(X1)", "9B_(X2)", 10:17, "18A", 19:28)) #Ignoring 18B

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
      pseudo.y.LG_x.arms <- cbind(pseudo.y.LG_x.arms, (1-as.numeric(LG_X[4,i])))
      y.All.LG.arms <- cbind(y.All.LG.arms, LG_X[3,i])
    }
    if (LG_X[2,i] == "a2"){
      pseudo.y.LG_x.arms <- cbind(pseudo.y.LG_x.arms, (1-as.numeric(LG_X[5,i])))
      y.All.LG.arms <- cbind(y.All.LG.arms, LG_X[3,i])
    }
    if (LG_X[2,i] == "c"){
      pseudo.y.LG_x.cen <- cbind(pseudo.y.LG_x.cen, (1-as.numeric(LG_X[4,i]))) #Centromeres are scored as the same phase as a1 to cover whole centromere in acrocentrics as well
    }
  }
  pseudo.y.All.LG.arms <- cbind(pseudo.y.All.LG.arms, pseudo.y.LG_x.arms)
  pseudo.y.All.LG <- cbind(pseudo.y.All.LG, pseudo.y.LG_x.arms, pseudo.y.LG_x.cen)
}

###Plotting y and pseudo y distribution barplots

##Making histogram counts to be plotted:

#Number of loci for counts with and without centromeres
no.loci.y.All.LG.arms <- length(y.All.LG.arms)
no.loci.y.All.LG <- length(y.All.LG)
no.loci.pseudo.y.All.LG.arms <- length(pseudo.y.All.LG.arms)
no.loci.pseudo.y.All.LG <- length(pseudo.y.All.LG)

#Individual historgrams:
bins <- seq(0,1,0.1)

H.y.All.LG.arms <- hist(as.numeric(y.All.LG.arms), breaks=bins)
H.y.All.LG <- hist(as.numeric(y.All.LG), breaks=bins)
H.pseudo.y.arms <- hist(pseudo.y.All.LG.arms, breaks=bins)
H.pseudo.y.All.LG <- hist(pseudo.y.All.LG, breaks=bins)

##Making barplots
#Setting universal parameters:
las = 1
width = 0.1
cex=0.7
xlab <- seq(0.1,1,0.1)
xlim=NULL #c(0,1)  
xlab1=NULL #"y"
xlab2=NULL #"Pseudo y"
cex.names=0.7
cex.axis=0.7
ylim1=c(0,0.4)  
ylim2=c(0,0.6)  

#Individual Barplots
dev.off()
pdf(file="y_and_pseudo_y_distributions.pdf")
par(mfrow=c(2,2),oma=c(2,1,2,0))

par(mar=c(2,4,2,0))
barplot(H.y.All.LG.arms$counts/no.loci.y.All.LG.arms, col="pink", names.arg=NULL,
        ylab="Proportion of loci", xlab=xlab1, ylim=ylim1, las=las, width=width, xlim=xlim,
        cex.axis=cex.axis)
box(col="gray80")
mtext(text="y \n (excl. centromeres)", side=3, line=-2.2, cex=cex, las=1)

par(mar=c(2,1,2,3))
barplot(H.y.All.LG$counts/no.loci.y.All.LG, col="darkred", names.arg=NULL, xlab=xlab1,
        ylim=ylim1, las=las, width=width, axes=F, xlim=xlim)
box(col="gray80")
mtext(text="y", side=3, line=-1.2, cex=cex, las=1)

par(mar=c(4,4,0,0))
barplot(H.pseudo.y.arms$counts/no.loci.pseudo.y.All.LG.arms, col="lightblue", names.arg=xlab,
        ylab="Proportion of loci", xlab=xlab2, ylim=ylim2, las=las, width=width, xlim=xlim,
        cex.names=cex.names, cex.axis=cex.axis)
box(col="gray80")
mtext(text="haplo y \n (excl. centromeres)", side=3, line=-2.2, cex=cex, las=1)

par(mar=c(4,1,0,3))
barplot(H.pseudo.y.All.LG$counts/no.loci.pseudo.y.All.LG, col="darkblue", names.arg=xlab,
        xlab=xlab2, ylim=ylim2, las=las, width=width, axes=F, xlim=xlim, cex.names=cex.names)
box(col="gray80")
mtext(text="haplo y", side=3, line=-1.2, cex=cex, las=1)

dev.off()

