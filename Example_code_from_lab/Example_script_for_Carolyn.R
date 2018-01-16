# Charlie's example code for Carolyn
# April 6, 2016

par()  ##view default par settings
opar<-par()   ##save a copy of default par settings so that I can revert back to them while plotting


# Import Fst values and save values for each pop (example of full code)
Fst_mapped <-read.csv("Fst/Fst_mapped_loci.csv",header=TRUE)
Fst_1998_2002INT_mapped <- Fst_mapped$Fst_1998_2002INT


# Plot Fst across the genome
qvals_summary <- read.csv("Qvalues_bayescan_ftemp_manhattenplot.csv",header=TRUE)  ##This file has cumulative map positions

box_starts<-c(190.03,495.16,778.51,1103.49,1411.39,1687.68,2014.15,2340.95,2542.96,2745.38,2928,3079.19,3233.45,3377.65,3526.45,3682.85,3908.19) #Colored boxes for chromosomes
box_ends<-c(318.91,652.67,945.25,1232.97,1552.29,1835.08,2186.92,2468.58,2623.57,2859.71,3009.7,3136.97,3299.6,3453.31,3615.91,3768.25,3970.65)

par(opar) #option to restore default settings

par(mfrow=c(3,2),mar=c(6,6,5,6),oma=c(7,1,2,2),mai=c(0.5,1,0.5,0.5),new=FALSE) #number of plots per window and space of margins

outlier_sum<- read.csv('C:/Users/Waters/Outlier_summary.csv', header=TRUE)

Bayescan_outs<-outlier_sum[outlier_sum$Outlier=="Bayescan",]  #2 loci for seg line; 11 for INT line (combine bayescan and Baye_Ftemp)
FtempSEG_outs<-outlier_sum[outlier_sum$Outlier=="FtempSEG",]  ##53 loci for SEG line
DAPC_structure<-outlier_sum[outlier_sum$Outlier=="DAPC",]     ##33 loci for SEG line; 37 for INT line (combine DAPC with Ftemp_DAPC)
Baye_Ftemp<-outlier_sum[outlier_sum$Outlier=="Bayescan_Ftemp",]  ##Only for SEG line, 9 loci
Ftemp_DAPC<-outlier_sum[outlier_sum$Outlier=="Ftemp_DAPC",] #4 loci for SEG line
All_methods<-outlier_sum[outlier_sum$Outlier=="All",]  #2 loci for SEG line; these 2 loci are Bayes_DAPC for INT line

sliding_window_outliers<-read.csv('C:/Users/Waters/sliding_window_outlier_regions_populations.csv', header=TRUE)
sw_2002int <- sliding_window_outliers[sliding_window_outliers$Population=="2002INT",] #give me rows where Population column equals 2002INT but give me all columns
sw_2002seg <- sliding_window_outliers[sliding_window_outliers$Population=="2002SEG",]
sw_2006int <- sliding_window_outliers[sliding_window_outliers$Population=="2006INT",]
sw_2006seg <- sliding_window_outliers[sliding_window_outliers$Population=="2006SEG",]
sw_2010int <- sliding_window_outliers[sliding_window_outliers$Population=="2010INT",]
sw_2010seg <- sliding_window_outliers[sliding_window_outliers$Population=="2010SEG",]


par(xpd=FALSE)
#2002SEG
plot(qvals_summary$Cumulative_position,qvals_summary$fst2002seg,ylim=c(-0.02,0.25),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,-0.025,box_ends,0.25,col="gray91",border=NA)
points(qvals_summary$Cumulative_position,qvals_summary$fst2002seg,ylim=c(-0.02,0.25),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(Bayescan_outs$Cumulative_position,Bayescan_outs$fst2002seg,type='p',pch=15,cex=1.2,col="blue")
points(FtempSEG_outs$Cumulative_position,FtempSEG_outs$fst2002seg,type='p',pch=16,cex=1.2,col="red")
points(DAPC_structure$Cumulative_position,DAPC_structure$fst2002seg,type='p',pch=17,cex=1.2,col="purple")
points(Baye_Ftemp$Cumulative_position,Baye_Ftemp$fst2002seg,type='p',pch=15,cex=1.5,col="orange")
points(Ftemp_DAPC$Cumulative_position,Ftemp_DAPC$fst2002seg,type='p',pch=16,cex=1.5,col="green")
points(All_methods$Cumulative_position,All_methods$fst2002seg,type='p',pch=17,cex=1.8,col="cyan")
rect(sw_2002seg$Cum_start,0.245,sw_2002seg$Cum_end,0.25,col='black',lwd=3)
mtext("2002SEG",outer=FALSE,line=0.5,cex=2,at=2000)

## 2002INT
plot(qvals_summary$Cumulative_position,qvals_summary$fst2002int,ylim=c(-0.02,0.25),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,-0.025,box_ends,0.25,col="gray91",border=NA)
points(qvals_summary$Cumulative_position,qvals_summary$fst2002int,ylim=c(-0.02,0.25),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(Bayescan_outs$Cumulative_position,Bayescan_outs$fst2002int,type='p',pch=15,cex=1.2,col="blue")  #Bayescan
points(Baye_Ftemp$Cumulative_position,Baye_Ftemp$fst2002int,type='p',pch=15,cex=1.2,col="blue")       #Also Bayescan only for INT line
points(DAPC_structure$Cumulative_position,DAPC_structure$fst2002int,type='p',pch=17,cex=1.2,col="purple")  #DAPC
points(Ftemp_DAPC$Cumulative_position,Ftemp_DAPC$fst2002int,type='p',pch=17,cex=1.2,col="purple")         #Also DAPC only for INT line
points(All_methods$Cumulative_position,All_methods$fst2002int,type='p',pch=17,cex=1.5,col="brown4") ##Just Bayescan and DAPC for INT line
rect(sw_2002int$Cum_start,0.245,sw_2002int$Cum_end,0.25,col='black',lwd=3)
mtext("2002INT",outer=FALSE,line=0.5,cex=2,at=2000)

#2006 SEG
plot(qvals_summary$Cumulative_position,qvals_summary$fst2006seg,ylim=c(-0.02,0.25),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,-0.025,box_ends,0.25,col="gray91",border=NA)
points(qvals_summary$Cumulative_position,qvals_summary$fst2006seg,ylim=c(-0.02,0.25),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(Bayescan_outs$Cumulative_position,Bayescan_outs$fst2006seg,type='p',pch=15,cex=1.2,col="blue")
points(FtempSEG_outs$Cumulative_position,FtempSEG_outs$fst2006seg,type='p',pch=16,cex=1.2,col="red")
points(DAPC_structure$Cumulative_position,DAPC_structure$fst2006seg,type='p',pch=17,cex=1.2,col="purple")
points(Baye_Ftemp$Cumulative_position,Baye_Ftemp$fst2006seg,type='p',pch=15,cex=1.5,col="orange")
points(Ftemp_DAPC$Cumulative_position,Ftemp_DAPC$fst2006seg,type='p',pch=16,cex=1.5,col="green")
points(All_methods$Cumulative_position,All_methods$fst2006seg,type='p',pch=17,cex=1.8,col="cyan")
rect(sw_2006seg$Cum_start,0.245,sw_2006seg$Cum_end,0.25,col='black',lwd=3)
mtext("2006SEG",outer=FALSE,line=0.5,cex=2,at=2000)

#2006INT
plot(qvals_summary$Cumulative_position,qvals_summary$fst2006int,ylim=c(-0.02,0.25),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,-0.025,box_ends,0.25,col="gray91",border=NA)
points(qvals_summary$Cumulative_position,qvals_summary$fst2006int,ylim=c(-0.02,0.25),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(Bayescan_outs$Cumulative_position,Bayescan_outs$fst2006int,type='p',pch=15,cex=1.2,col="blue")  #Bayescan
points(Baye_Ftemp$Cumulative_position,Baye_Ftemp$fst2006int,type='p',pch=15,cex=1.2,col="blue")       #Also Bayescan only for INT line
points(DAPC_structure$Cumulative_position,DAPC_structure$fst2006int,type='p',pch=17,cex=1.2,col="purple")  #DAPC
points(Ftemp_DAPC$Cumulative_position,Ftemp_DAPC$fst2006int,type='p',pch=17,cex=1.2,col="purple")         #Also DAPC only for INT line
points(All_methods$Cumulative_position,All_methods$fst2006int,type='p',pch=17,cex=1.5,col="brown4") ##Just Bayescan and DAPC for INT line
rect(sw_2006int$Cum_start,0.245,sw_2006int$Cum_end,0.25,col='black',lwd=3)
mtext("2006INT",outer=FALSE,line=0.5,cex=2,at=2000)

#2010SEG
plot(qvals_summary$Cumulative_position,qvals_summary$fst2010seg,ylim=c(-0.02,0.25),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,-0.025,box_ends,0.25,col="gray91",border=NA)
points(qvals_summary$Cumulative_position,qvals_summary$fst2010seg,ylim=c(-0.02,0.25),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(Bayescan_outs$Cumulative_position,Bayescan_outs$fst2010seg,type='p',pch=15,cex=1.2,col="blue")
points(FtempSEG_outs$Cumulative_position,FtempSEG_outs$fst2010seg,type='p',pch=16,cex=1.2,col="red")
points(DAPC_structure$Cumulative_position,DAPC_structure$fst2010seg,type='p',pch=17,cex=1.2,col="purple")
rect(sw_2010seg$Cum_start,0.245,sw_2010seg$Cum_end,0.25,col='black',lwd=3)
points(Baye_Ftemp$Cumulative_position,Baye_Ftemp$fst2010seg,type='p',pch=15,cex=1.5,col="orange")
points(Ftemp_DAPC$Cumulative_position,Ftemp_DAPC$fst2010seg,type='p',pch=16,cex=1.5,col="green")
points(All_methods$Cumulative_position,All_methods$fst2010seg,type='p',pch=17,cex=1.8,col="cyan")
mtext("2010SEG",outer=FALSE,line=0.5,cex=2,at=2000)

#2010INT
plot(qvals_summary$Cumulative_position,qvals_summary$fst2010int,ylim=c(-0.02,0.25),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,-0.025,box_ends,0.25,col="gray91",border=NA)
points(qvals_summary$Cumulative_position,qvals_summary$fst2010int,ylim=c(-0.02,0.25),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(Bayescan_outs$Cumulative_position,Bayescan_outs$fst2010int,type='p',pch=15,cex=1.2,col="blue")  #Bayescan
points(Baye_Ftemp$Cumulative_position,Baye_Ftemp$fst2010int,type='p',pch=15,cex=1.2,col="blue")       #Also Bayescan only for INT line
points(DAPC_structure$Cumulative_position,DAPC_structure$fst2010int,type='p',pch=17,cex=1.2,col="purple")  #DAPC
points(Ftemp_DAPC$Cumulative_position,Ftemp_DAPC$fst2010int,type='p',pch=17,cex=1.2,col="purple")         #Also DAPC only for INT line
points(All_methods$Cumulative_position,All_methods$fst2010int,type='p',pch=17,cex=1.5,col="brown4") ##Just Bayescan and DAPC for INT line
rect(sw_2010int$Cum_start,0.245,sw_2010int$Cum_end,0.25,col='black',lwd=3)
mtext("2010INT",outer=FALSE,line=0.5,cex=2,at=2000)

mtext("Cumulative Map Position (cM)",outer=TRUE,side=1,line=1.5,cex=2,at=0.525)
mtext(expression(paste("F"[ST])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Bayescan", "Ftemp", "DAPC", "Baye&Ftemp","Baye&DAPC","DAPC&Ftemp", "All Methods","Sliding Window"), xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n", pch = c(15,16,17,17,15,16,17,NA), lty=c(NA,NA,NA,NA,NA,NA,NA,1),lwd=c(NA,NA,NA,NA,NA,NA,NA,3),col =c("blue","red","purple","brown4","orange","green","cyan","black"), cex = 2)
