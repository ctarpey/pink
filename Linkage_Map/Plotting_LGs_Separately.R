#######################################
### Running plot code for single LG ###
#######################################

rm(list=ls())
library(KernSmooth)

############################
## DEFINING DATA TO PLOT  ##
############################
#Kokanee Final
setwd("C:/Users/moli/Google Drev/Morten arbejde/Sockeye vs Kokanee/Paper 1-2/Pops/Data/Kokanee")
Fst_map <- read.csv(file = 'Kokanee_Map_w_data.csv', sep = ',')
Boot_results <- read.csv(file = 'Kokanee_Bootstrap_results.csv', sep = ',')

#PSV map LG names:
LGs <- as.character(unique(Fst_map$LG_nerka))


###########################
## Defining data to plot ##
###########################

pop.comparisons <- colnames(Fst_map)[7:21]
data_to_plot <- pop.comparisons[15] #Defining pop comparison to plot:
file = paste(data_to_plot, "_Boot_kernel.pdf", sep="") #Name of output file

###############################
## FIXED PLOTTING PARAMETERS ##
###############################

boot.mean.to.plot <- 'BSM.mean'
boot.p.to.plot <- 'p.value.4N.corr'

label1 <- paste(data_to_plot, "_2n", sep="")
label2 <- paste(data_to_plot, "_4n", sep="")

LGs_Y_axis <- LGs[c(1,7,13,19,25)] #LGs for which y-axes will be drawn:
LGs_X_axis <- LGs[c(24:30)]        #LGs for which x-axes will be drawn:

ylim=c(-0.14,0.65)
ylab = expression(paste("1-BSM"))
#ylab = expression(paste("1-D"[BSM]))
#ylab = expression(paste(italic("F")[ST]," / 1-D"[BSM]))
#ylab = expression(paste(italic("F")[ST]))
xlab = "cM"

#col=c("#00FF0080", "#0000FF80") #GREEN & BLUE
col=c("#0000FF40", "#FD3C0640") #Ryan's colors

BSM.limit <- 0.000
bandwidth = 3.5


##################################
## START LOOP over ALL LGs HERE ##
##################################
#Writing pdf file:
pdf(file=file)

par(mfrow=c(5,6),oma=c(4,4,3,1), mar=c(1.5,1,1,0))
for (j in LGs){
    
  ##################################
  ## Subset All Map data for LG j ##
  ##################################
  #j <- c("1")
  LG_j_loci <- which(Fst_map$LG_nerka == j)
  LG_j_Map <- Fst_map[LG_j_loci,]
  
  # Define LG j data for plotting linkage map below:
  LG_j_cM <- LG_j_Map$cM
  LG_j_color <- LG_j_Map$color
  
  # Define LG j data for plotting BSM values:
  LG_j_data <- cbind(LG_j_Map$cM, LG_j_Map[,data_to_plot], LG_j_Map$color) #DATA SELECTOR#
  LG_j_data_excl_NA <- na.omit(LG_j_data)

  #Make plotting frame for LG j:
  plot(x=NULL, y=NULL, xlim=c(-5,230), axes=F,ann=F, ylim=ylim, xaxs="i",yaxs="i")
  mtext(text=paste("So", j, sep=""), side=3, line=0.3, cex=0.65)
  box(col="gray80")
  abline(h=-0.05, lty=1, col="gray80")
  axis(1, las=1, cex.axis=0.5, col="gray80", labels=FALSE)  #reduce the number of axis labels  
  axis(2, las=1, cex.axis=0.5, col="gray80", at=seq(0,ylim[2], 0.2), labels=FALSE)
  
  if (j %in% LGs_Y_axis) {     #true when i is 1, 5, 9, 13 (%% is remainder when i divided by 4)
    axis(2, las=1, cex.axis=0.5, col="gray80", at=seq(0,ylim[2], 0.2))
    mtext(text=ylab, side=2, line=1.9, cex=0.3, las=0)
  }
  
  if (j %in% LGs_X_axis){
    axis(1, las=1, cex.axis=0.5, col="gray80")
    mtext(text=xlab, side=1, line=1.9, cex=0.3)
  }


  #Plot mapped loci
  points(x=LG_j_cM, y=rep(-0.09, length(LG_j_cM)), type="p", pch="|", col=col[LG_j_color], cex=0.3)
  
  #Defining centromeric region:
  Centromeric_pos <- NULL #Making vector of all LG positions in centromere
  for (i in 1:nrow(LG_j_Map)){ 
    if (LG_j_Map$LG_nerka[i]==j){
      if (LG_j_Map$Chr_Arm[i] == "c"){
        Centromeric_pos <- cbind(Centromeric_pos, as.numeric(as.character(LG_j_Map$cM[i])))    
      }  
    }
  }
  
  centro_start <- as.numeric(min(Centromeric_pos))
  centro_end <- as.numeric(max(Centromeric_pos))
  if (centro_end - centro_start < 2){
    centro_start = centro_start-1
    centro_end = centro_end+1
  }
  centro_x <- c(centro_start, centro_end, centro_end, centro_start)
  centro_y <- c(-0.05, -0.05, 1.2, 1.2)
  polygon(x=centro_x, y=centro_y, col="#D3D3D330", border= NA)
  
  ##Fitting a smoothed kernel using locpoly()
  filtered.loci <- which(LG_j_data_excl_NA[,2] <= 1-BSM.limit)
  fit <- locpoly(x=LG_j_data_excl_NA[,1][filtered.loci], y=1-LG_j_data_excl_NA[,2][filtered.loci], bandwidth = bandwidth)
  #lines(fit, lwd=1.5, col="#00000080")
  
  ##Loop over markers in LG j to Plot data  
  for (i in 1:nrow(LG_j_data_excl_NA)){  
    if ((1-LG_j_data_excl_NA[i,2]) > BSM.limit){  
      points(x=LG_j_data_excl_NA[i,1], y=1-LG_j_data_excl_NA[i,2], type="p", pch=19, cex=0.4, col=col[LG_j_data_excl_NA[i,3]])    
    }
  }
}

##Plot legend in lower right corner
plot(x=NULL, y=NULL, xlim=c(-5,230), axes=F,ann=F, ylim=ylim)
text(x=220, y=ylim[2]*0.50, labels=label1, col="#0000FF90", cex=0.7, adj=1)
text(x=220, y=ylim[2]*0.35, labels=label2, col="#FD3C0690", cex=0.7, adj=1)

dev.off()

###############
### THE END ###
###############
