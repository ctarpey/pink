### Manhattan plots of Outliers from LFMM and Arlequin
###   Using code that Charlie Waters wrote
### Carolyn Tarpey | April 2016
### ---------------------------------------

setwd('G:/Analysis/Mapping/Outliers_Manhattan_Plot/Arlequin_LFMM')

par()  ##view default par settings
opar<-par()   ##save a copy of default par settings so that I can revert back to them while plotting

# Import Fst values and save values for each pop (example of full code)
Fst_mapped <-read.table("Outliers_Master.txt", header=TRUE) ##This file has cumulative map positions
head(Fst_mapped)

#Colored boxes for even chromosomes
box_starts<-c(125.18, 374.86, 610.3, 836.07, 1074.26, 1303.02, 1522.59, 1765.54, 1996.05, 2247.89, 2463.55, 2692.96, 2922.46)
box_ends<-c(251.86, 495.05, 728.59, 952.17, 1189.4, 1411.97, 1642.5, 1893.93, 2123.51, 2352.32, 2560.11, 2807.45, 3037.61)


#these make breakout groups for each of the outlier types to overplot on FST
LFMM_odd<-Fst_mapped[Fst_mapped$odd_any=="yes",]  #outlier in any LFMM odd test
LFMM_o_PT1<-Fst_mapped[Fst_mapped$odd_PT1=="yes",] 
LFMM_o_PT2<-Fst_mapped[Fst_mapped$odd_PT2=="yes",]
LFMM_o_LO<-Fst_mapped[Fst_mapped$odd_LO=="yes",]
LFMM_o_LA<-Fst_mapped[Fst_mapped$odd_LA=="yes",]

LFMM_even<-Fst_mapped[Fst_mapped$even_any=="yes",]  ##outlier in any LFMM even test
LFMM_e_PT1<-Fst_mapped[Fst_mapped$even_PT1=="yes",]  
LFMM_e_PT2<-Fst_mapped[Fst_mapped$even_PT2=="yes",]
LFMM_e_LO<-Fst_mapped[Fst_mapped$even_LO=="yes",]
LFMM_e_LA<-Fst_mapped[Fst_mapped$even_LA=="yes",]

LFmm_overlap<-Fst_mapped[Fst_mapped$LFMM_oe_overlap=="yes",] #outlier in any LFMM odd and any LFMM even
LFMM_PT1<-Fst_mapped[Fst_mapped$overlap_PT1=="yes",]  # PT1 outlier in both odd and even 
LFMM_PT2<-Fst_mapped[Fst_mapped$overlap_PT2=="yes",]
LFMM_LO<-Fst_mapped[Fst_mapped$overlap_LO=="yes",]
LFMM_LA<-Fst_mapped[Fst_mapped$overlap_LA=="yes",]

Arlequin_odd<-Fst_mapped[Fst_mapped$odd_0.01=="yes",]  #outlier in Arlequin at sig level 0.01 odd
Arlequin_even<-Fst_mapped[Fst_mapped$even_0.01 =="yes",] 
Arlequin_overlap <- Fst_mapped[Fst_mapped$overlap_0.01 =="yes",]  #outlier in both arlequin 

Overlap_odd<-Fst_mapped[Fst_mapped$LFMM_arl_odd=="yes",]  #outlier in any LFMM odd and odd Arlequin
Overlap_even<-Fst_mapped[Fst_mapped$LFMM_arl_even=="yes",] #outlier in any LFMM even and even Arlequin

Any_outlier<-Fst_mapped[Fst_mapped$any_outlier=="yes",]   #outlier in any test

onlyOverlap_odd<-Fst_mapped[Fst_mapped$onlyLFMM_arl_odd=="yes",]  #outlier only in LFMM odd and odd Arlequin
onlyOverlap_even<-Fst_mapped[Fst_mapped$onlyLFMM_arl_even=="yes",] #outlier only in LFMM even and even Arlequin
Arlequin_LFMM_overlap<-Fst_mapped[Fst_mapped$LFMM_arl_overlap=="yes",]   #outlier overlapping between arl/lfmm odd and arl/lfmm even


par(opar) #option to restore default settings

#cols_OddEven <-c(AMUR10 = "#EE547D", AMUR11 = "#4E70C8", HAYLY09 = "#4C7DE7", HAYLY10 = "#ED2DBE", 
#                 KOPPE96 = "#D186CF", KOPPE91 = "#B0B0D8", KUSHI06 = "#AD2868", KUSHI07 = "#354A7C",
#                 NOME91 = "#95E0CA", NOME94 = "#E63D25", SNOH03 = "#708FEC", SNOH96 ="#a12787" , TAUY09 = "#3C92A8", TAUY12 = "#DF9E39")




######################### Overlap in Arlequin

pdf("Plots/Arlequin Overlap.pdf", width = 12, height = 8)

plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab= "Global Fst", xlab="Cumulative Map Position (cM)",cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab= "Fst", xlab="Cumulative Map Position (cM)", pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(Arlequin_overlap$BP_seq,Arlequin_overlap$global_Fst,type='p',pch=19,cex=.8,col="magenta")
mtext("Arlequin outliers in both odd and even-year lineage",outer=FALSE,line=0.5,cex=1.5,at=1500)

dev.off()

######################### Overlap in any LFMM odd and even

pdf("Plots/Overlap LFMM any odd any even.pdf", width = 12, height = 8)

plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab= "Global Fst", xlab="Cumulative Map Position (cM)",cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab= "Fst", xlab="Cumulative Map Position (cM)", pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFmm_overlap$BP_seq,LFmm_overlap$global_Fst,type='p',pch=19,cex=.8,col="coral")
mtext("LFMM ouliers overlap any test odd and even-year lineage",outer=FALSE,line=0.5,cex=1.5,at=1500)

dev.off()

######################### Outlier overlap between LFMM and Arl 1x2 three colors

pdf("Plots/Outlier overlap Arl and LFMM_4.pdf", width = 12, height = 7)

plot(onlyOverlap_even$BP_seq,onlyOverlap_even$global_Fst,ylim=c(-0.02,0.7),
     type="n",col="gray",ylab= "Global Fst", xlab="Cumulative Map Position (cM)",cex.axis=1.5,cex.lab=1.5)
rect(box_starts,0,box_ends,0.7,col="gray91",border=NA)
abline(h=0)
points(onlyOverlap_even$BP_seq,onlyOverlap_even$global_Fst,type='p',pch=19,cex=.8,col="#AD2868") #any even LFMM test outlier
points(onlyOverlap_odd$BP_seq,onlyOverlap_odd$global_Fst,type='p',pch=19,cex=.8,col="navy") #any odd LFMM test outlier
points(Arlequin_LFMM_overlap$BP_seq,Arlequin_LFMM_overlap$global_Fst,type='p',pch=19,cex=.8,col="green") #any odd LFMM test outlier
mtext("Overlapping LFMM & Arlequin Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)

legend("bottom", c("Odd-year Outliers", "Even-year Outliers", "Outliers in both lineages"), xpd = TRUE, horiz = TRUE, inset = c(0,0),
       bty = "n", pch = 19,cex=1.2,lty=c(NA,NA,NA),lwd=c(NA,NA,NA),col = c("navy" ,"#AD2868",  "green"))

dev.off()


######################### Outlier overlap between LFMM and Arl 1x2

pdf("Plots/Outlier overlap by lineage Arl and LFMM.pdf", width = 12, height = 8)

#number of plots per window and space of margins
par(mfrow=c(1,2),mar=c(6,6,5,6),oma=c(7,1,2,2),mai=c(0.5,1,0.5,0.5),new=FALSE)

par(xpd=FALSE)

#LFMM_odd any
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',xlab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',xlab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(Overlap_odd$BP_seq,Overlap_odd$global_Fst,type='p',pch=19,cex=.8,col="navy") #any odd LFMM test outlier
mtext("Overlapping LFMM & Arlequin Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)

#LFMM_even any
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',xlab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',xlab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(Overlap_even$BP_seq,Overlap_even$global_Fst,type='p',pch=19,cex=.8,col="#AD2868") #any even LFMM test outlier
mtext("Overlapping LFMM & Arlequin Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)

mtext("Cumulative Map Position (cM)",outer=TRUE,side=1,line=1.5,cex=2,at=0.525)
mtext(expression(paste("Global F"[ST])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Odd-year Outliers", "Even-year Outliers"), xpd = TRUE, horiz = TRUE, inset = c(0,0),
       bty = "n", pch = 19, lty=c(NA,NA),lwd=c(NA,NA),col =c("navy" ,"#AD2868"), cex = 2)

dev.off()

######################### Raw LFMM, 4x2 

pdf("Plots/LFMM by test by lineage.pdf", width = 8, height = 12)

#number of plots per window and space of margins
par(mfrow=c(4,2),mar=c(6,6,5,6),oma=c(7,1,2,2),mai=c(0.5,1,0.5,0.5),new=FALSE)

par(xpd=FALSE)

#LFMM_odd PT1
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',xlab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',xlab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_o_PT1$BP_seq,LFMM_o_PT1$global_Fst,type='p',ylab='',pch=19,cex=.8,col="red")
mtext("Odd-year PT1 Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)

#LFMM_even PT1
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',xlab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',xlab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_e_PT1$BP_seq,LFMM_e_PT1$global_Fst,type='p',ylab='',pch=19,cex=.8,col="red")
mtext("Even-year PT1 Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)

#LFMM_odd PT2 
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',xlab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',xlab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_o_PT2$BP_seq,LFMM_o_PT2$global_Fst,type='p',ylab='',pch=19,cex=.8,col="blue")
mtext("Odd-year PT2 Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)

#LFMM_even PT2
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',xlab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',xlab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_e_PT2$BP_seq,LFMM_e_PT2$global_Fst,type='p',ylab='',pch=19,cex=.8,col="blue")
mtext("Even-year PT2 Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)


#LFMM_odd LO 
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',xlab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',xlab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_o_LO$BP_seq,LFMM_o_LO$global_Fst,type='p',ylab='',pch=19,cex=.8,col="orange")
mtext("Odd-year LO Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)

#LFMM_even LO
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',xlab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',xlab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_e_LO$BP_seq,LFMM_e_LO$global_Fst,type='p',ylab='',pch=19,cex=.8,col="orange")
mtext("Even-year LO Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)


#LFMM_odd LA
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',xlab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',xlab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_o_LA$BP_seq,LFMM_o_LA$global_Fst,type='p',ylab='',pch=19,cex=.8,col="darkgreen")
mtext("Odd-year LA Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)

#LFMM_even LA
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',xlab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',xlab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_e_LA$BP_seq,LFMM_e_LA$global_Fst,type='p',ylab='',pch=19,cex=.8,col="darkgreen")
mtext("Even-year LA Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)


mtext("Cumulative Map Position (cM)",outer=TRUE,side=1,line=1.5,cex=2,at=0.525)
mtext(expression(paste("Global F"[ST])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("PT1 Outliers", "PT2 Outliers", "Longitude Outliers", "Latitude Outliers"), xpd = TRUE, horiz = TRUE, inset = c(0,0),
       bty = "n", pch = 19, lty=c(NA,NA,NA,NA),lwd=c(NA,NA,NA,NA),col =c("red","blue","orange","darkgreen"),  cex = 1.2)

dev.off()



######################### ANY, 2x2 

pdf("Plots/Outlier in any test by lineage.pdf", width = 12, height = 8)

#number of plots per window and space of margins
par(mfrow=c(2,2),mar=c(6,6,5,6),oma=c(7,1,2,2),mai=c(0.5,1,0.5,0.5),new=FALSE)

par(xpd=FALSE)

#LFMM_odd any
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',xlab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',xlab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_odd$BP_seq,LFMM_odd$global_Fst,type='p',pch=19,cex=.8,col="tomato") #any odd LFMM test outlier
mtext("LFMM Odd-year Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)

#LFMM_even any
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',xlab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',xlab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_even$BP_seq,LFMM_even$global_Fst,type='p',pch=19,cex=.8,col="tomato") #any even LFMM test outlier
mtext("LFMM Even-year Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)

#Arlequin_odd
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',xlab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',xlab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(Arlequin_odd$BP_seq,Arlequin_odd$global_Fst,type='p',pch=19,cex=.8,col="darkviolet") 
mtext("Arlequin Odd-year Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)

#Arlequin_even
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',xlab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',xlab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(Arlequin_even$BP_seq,Arlequin_even$global_Fst,type='p',pch=19,cex=.8,col="darkviolet") 
mtext("Arlequin Even-year Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)

mtext("Cumulative Map Position (cM)",outer=TRUE,side=1,line=1.5,cex=2,at=0.525)
mtext(expression(paste("Global F"[ST])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Outliers Any LFmm Test", "Outliers in Arlequin"), xpd = TRUE, horiz = TRUE, inset = c(0,0),
       bty = "n", pch = 19, lty=c(NA,NA),lwd=c(NA,NA),col =c("tomato","darkviolet"), cex = 2)

dev.off()

######################### LFMM overlap, 2x2

pdf("Plots/LFMM overlap by env test.pdf", width = 12, height = 8)

#number of plots per window and space of margins
par(mfrow=c(2,2),mar=c(6,6,5,6),oma=c(7,1,2,2),mai=c(0.5,1,0.5,0.5),new=FALSE)

par(xpd=FALSE)

#LFMM_odd overlap
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',xlab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',xlab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_PT1$BP_seq,LFMM_PT1$global_Fst,type='p',pch=19,cex=.75,col="red")
mtext("LFMM PT1 Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)

#LFMM_odd overlap
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_PT2$BP_seq,LFMM_PT2$global_Fst,type='p',pch=19,cex=.75,col="blue")
mtext("LFMM PT2 Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)

#LFMM_odd overlap
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',xlab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',xlab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_LO$BP_seq,LFMM_LO$global_Fst,type='p',pch=19,cex=.75,col="orange")
mtext("LFMM LO Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)

#LFMM_odd overlap
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',xlab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',xlab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_LA$BP_seq,LFMM_LA$global_Fst,type='p',pch=19,cex=.75,col="darkgreen")
mtext("LFMM LA Outliers",outer=FALSE,line=0.5,cex=1.5,at=1500)

mtext("Cumulative Map Position (cM)",outer=TRUE,side=1,line=1.5,cex=2,at=0.525)
mtext(expression(paste("Global F"[ST])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("PT1 Outliers", "PT2 Outliers", "Longitude Outliers", "Latitude Outliers"), xpd = TRUE, horiz = TRUE, inset = c(0,0),
       bty = "n", pch = 19, lty=c(NA,NA,NA,NA),lwd=c(NA,NA,NA,NA),col =c("red","blue","orange","darkgreen"), cex = 1)

dev.off()





######################### LFMM overlap, 2x2

#Arlequin_odd
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(Arlequin_odd$BP_seq,Arlequin_odd$global_Fst,type='p',pch=19,cex=.8,col="cyan") 
rect(sw_2002seg$Cum_start,0.245,sw_2002seg$Cum_end,0.25,col='black',lwd=3)
mtext("Arlequin odd year outliers",outer=FALSE,line=0.5,cex=2,at=2000)

#Arlequin_even
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(Arlequin_even$BP_seq,Arlequin_even$global_Fst,type='p',pch=19,cex=.8,col="cyan") 
rect(sw_2002seg$Cum_start,0.245,sw_2002seg$Cum_end,0.25,col='black',lwd=3)
mtext("Arlequin even year outliers",outer=FALSE,line=0.5,cex=2,at=2000)

#Arlequin_overlap
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="gray",ylab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="gray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(Arlequin_overlap$BP_seq,Arlequin_overlap$global_Fst,type='p',pch=19,cex=.8,col="cyan") 
rect(sw_2002seg$Cum_start,0.245,sw_2002seg$Cum_end,0.25,col='black',lwd=3)
mtext("Arlequin outliers in both lineages",outer=FALSE,line=0.5,cex=2,at=2000)


mtext("Cumulative Map Position (cM)",outer=TRUE,side=1,line=1.5,cex=2,at=0.525)
mtext(expression(paste("Global F"[ST])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Bayescan", "Ftemp", "DAPC", "Baye&Ftemp","Baye&DAPC","DAPC&Ftemp", "All Methods","Sliding Window"), xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n", pch = c(15,16,17,17,15,16,17,NA), lty=c(NA,NA,NA,NA,NA,NA,NA,1),lwd=c(NA,NA,NA,NA,NA,NA,NA,3),col =c("blue","red","purple","brown4","orange","green","cyan","black"), cex = 2)









#number of plots per window and space of margins
par(mfrow=c(4,2),mar=c(6,6,5,6),oma=c(7,1,2,2),mai=c(0.5,1,0.5,0.5),new=FALSE)

# Plot global Fst across the genome

par(xpd=FALSE)

#LFMM_odd
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_o_PT1$BP_seq,LFMM_o_PT1$global_Fst,type='p',pch=19,cex=.5,col="red")
points(LFMM_o_PT2$BP_seq,LFMM_o_PT2$global_Fst,type='p',pch=9,cex=.5,col="purple")
points(LFMM_o_LO$BP_seq,LFMM_o_LO$global_Fst,type='p',pch=19,cex=.5,col="orange")
points(LFMM_o_LA$BP_seq,LFMM_o_LA$global_Fst,type='p',pch=19,cex=.5,col="green")
rect(sw_2002seg$Cum_start,0.245,sw_2002seg$Cum_end,0.25,col='black',lwd=3)
mtext("LFMM odd year outliers",outer=FALSE,line=0.5,cex=2,at=2000)

#LFMM_odd
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_odd$BP_seq,LFMM_odd$global_Fst,type='p',pch=19,cex=.5,col="cyan") #any odd LFMM test outlier
rect(sw_2002seg$Cum_start,0.245,sw_2002seg$Cum_end,0.25,col='black',lwd=3)
mtext("LFMM odd year outliers",outer=FALSE,line=0.5,cex=2,at=2000)

#LFMM_even
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_e_PT1$BP_seq,LFMM_e_PT1$global_Fst,type='p',pch=19,cex=.5,col="red")
points(LFMM_e_PT2$BP_seq,LFMM_e_PT2$global_Fst,type='p',pch=9,cex=.5,col="purple")
points(LFMM_e_LO$BP_seq,LFMM_e_LO$global_Fst,type='p',pch=19,cex=.5,col="orange")
points(LFMM_e_LA$BP_seq,LFMM_e_LA$global_Fst,type='p',pch=19,cex=.5,col="green")
#points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,type='p',pch=17,cex=1.8,col="cyan")
rect(sw_2002seg$Cum_start,0.245,sw_2002seg$Cum_end,0.25,col='black',lwd=3)
mtext("LFMM even year outliers",outer=FALSE,line=0.5,cex=2,at=2000)

#LFMM_even
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(LFMM_even$BP_seq,LFMM_even$global_Fst,type='p',pch=19,cex=.5,col="cyan") #any even LFMM test outlier
rect(sw_2002seg$Cum_start,0.245,sw_2002seg$Cum_end,0.25,col='black',lwd=3)
mtext("LFMM even year outliers",outer=FALSE,line=0.5,cex=2,at=2000)


#Arlequin_odd
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(Arlequin_odd$BP_seq,Arlequin_odd$global_Fst,type='p',pch=19,cex=.5,col="cyan") 
rect(sw_2002seg$Cum_start,0.245,sw_2002seg$Cum_end,0.25,col='black',lwd=3)
mtext("Arlequin odd year outliers",outer=FALSE,line=0.5,cex=2,at=2000)

#Arlequin_even
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(Arlequin_even$BP_seq,Arlequin_even$global_Fst,type='p',pch=19,cex=.5,col="cyan") 
rect(sw_2002seg$Cum_start,0.245,sw_2002seg$Cum_end,0.25,col='black',lwd=3)
mtext("Arlequin even year outliers",outer=FALSE,line=0.5,cex=2,at=2000)

#Arlequin_overlap
plot(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
abline(h=0)
rect(box_starts,0,box_ends,0.9,col="gray91",border=NA)
points(Fst_mapped$BP_seq,Fst_mapped$global_Fst,ylim=c(-0.02,0.9),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
points(Arlequin_overlap$BP_seq,Arlequin_overlap$global_Fst,type='p',pch=19,cex=.5,col="cyan") 
rect(sw_2002seg$Cum_start,0.245,sw_2002seg$Cum_end,0.25,col='black',lwd=3)
mtext("Arlequin outliers in both lineages",outer=FALSE,line=0.5,cex=2,at=2000)


mtext("Cumulative Map Position (cM)",outer=TRUE,side=1,line=1.5,cex=2,at=0.525)
mtext(expression(paste("F"[ST])),outer=TRUE,side=2,line=-2,cex=2,at=0.525)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("bottom", c("Bayescan", "Ftemp", "DAPC", "Baye&Ftemp","Baye&DAPC","DAPC&Ftemp", "All Methods","Sliding Window"), xpd = TRUE, horiz = TRUE, inset = c(0,0), bty = "n", pch = c(15,16,17,17,15,16,17,NA), lty=c(NA,NA,NA,NA,NA,NA,NA,1),lwd=c(NA,NA,NA,NA,NA,NA,NA,3),col =c("blue","red","purple","brown4","orange","green","cyan","black"), cex = 2)




#################################################original code from charlie# 
# ## 2002INT
# plot(qvals_summary$Cumulative_position,qvals_summary$fst2002int,ylim=c(-0.02,0.25),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
# abline(h=0)
# rect(box_starts,-0.025,box_ends,0.25,col="gray91",border=NA)
# points(qvals_summary$Cumulative_position,qvals_summary$fst2002int,ylim=c(-0.02,0.25),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
# points(Bayescan_outs$Cumulative_position,Bayescan_outs$fst2002int,type='p',pch=15,cex=1.2,col="blue")  #Bayescan
# points(Baye_Ftemp$Cumulative_position,Baye_Ftemp$fst2002int,type='p',pch=15,cex=1.2,col="blue")       #Also Bayescan only for INT line
# points(DAPC_structure$Cumulative_position,DAPC_structure$fst2002int,type='p',pch=17,cex=1.2,col="purple")  #DAPC
# points(Ftemp_DAPC$Cumulative_position,Ftemp_DAPC$fst2002int,type='p',pch=17,cex=1.2,col="purple")         #Also DAPC only for INT line
# points(All_methods$Cumulative_position,All_methods$fst2002int,type='p',pch=17,cex=1.5,col="brown4") ##Just Bayescan and DAPC for INT line
# rect(sw_2002int$Cum_start,0.245,sw_2002int$Cum_end,0.25,col='black',lwd=3)
# mtext("2002INT",outer=FALSE,line=0.5,cex=2,at=2000)
# 
# #2006 SEG
# plot(qvals_summary$Cumulative_position,qvals_summary$fst2006seg,ylim=c(-0.02,0.25),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
# abline(h=0)
# rect(box_starts,-0.025,box_ends,0.25,col="gray91",border=NA)
# points(qvals_summary$Cumulative_position,qvals_summary$fst2006seg,ylim=c(-0.02,0.25),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
# points(Bayescan_outs$Cumulative_position,Bayescan_outs$fst2006seg,type='p',pch=15,cex=1.2,col="blue")
# points(FtempSEG_outs$Cumulative_position,FtempSEG_outs$fst2006seg,type='p',pch=16,cex=1.2,col="red")
# points(DAPC_structure$Cumulative_position,DAPC_structure$fst2006seg,type='p',pch=17,cex=1.2,col="purple")
# points(Baye_Ftemp$Cumulative_position,Baye_Ftemp$fst2006seg,type='p',pch=15,cex=1.5,col="orange")
# points(Ftemp_DAPC$Cumulative_position,Ftemp_DAPC$fst2006seg,type='p',pch=16,cex=1.5,col="green")
# points(All_methods$Cumulative_position,All_methods$fst2006seg,type='p',pch=17,cex=1.8,col="cyan")
# rect(sw_2006seg$Cum_start,0.245,sw_2006seg$Cum_end,0.25,col='black',lwd=3)
# mtext("2006SEG",outer=FALSE,line=0.5,cex=2,at=2000)
# 
# #2006INT
# plot(qvals_summary$Cumulative_position,qvals_summary$fst2006int,ylim=c(-0.02,0.25),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
# abline(h=0)
# rect(box_starts,-0.025,box_ends,0.25,col="gray91",border=NA)
# points(qvals_summary$Cumulative_position,qvals_summary$fst2006int,ylim=c(-0.02,0.25),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
# points(Bayescan_outs$Cumulative_position,Bayescan_outs$fst2006int,type='p',pch=15,cex=1.2,col="blue")  #Bayescan
# points(Baye_Ftemp$Cumulative_position,Baye_Ftemp$fst2006int,type='p',pch=15,cex=1.2,col="blue")       #Also Bayescan only for INT line
# points(DAPC_structure$Cumulative_position,DAPC_structure$fst2006int,type='p',pch=17,cex=1.2,col="purple")  #DAPC
# points(Ftemp_DAPC$Cumulative_position,Ftemp_DAPC$fst2006int,type='p',pch=17,cex=1.2,col="purple")         #Also DAPC only for INT line
# points(All_methods$Cumulative_position,All_methods$fst2006int,type='p',pch=17,cex=1.5,col="brown4") ##Just Bayescan and DAPC for INT line
# rect(sw_2006int$Cum_start,0.245,sw_2006int$Cum_end,0.25,col='black',lwd=3)
# mtext("2006INT",outer=FALSE,line=0.5,cex=2,at=2000)
# 
# #2010SEG
# plot(qvals_summary$Cumulative_position,qvals_summary$fst2010seg,ylim=c(-0.02,0.25),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
# abline(h=0)
# rect(box_starts,-0.025,box_ends,0.25,col="gray91",border=NA)
# points(qvals_summary$Cumulative_position,qvals_summary$fst2010seg,ylim=c(-0.02,0.25),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
# points(Bayescan_outs$Cumulative_position,Bayescan_outs$fst2010seg,type='p',pch=15,cex=1.2,col="blue")
# points(FtempSEG_outs$Cumulative_position,FtempSEG_outs$fst2010seg,type='p',pch=16,cex=1.2,col="red")
# points(DAPC_structure$Cumulative_position,DAPC_structure$fst2010seg,type='p',pch=17,cex=1.2,col="purple")
# rect(sw_2010seg$Cum_start,0.245,sw_2010seg$Cum_end,0.25,col='black',lwd=3)
# points(Baye_Ftemp$Cumulative_position,Baye_Ftemp$fst2010seg,type='p',pch=15,cex=1.5,col="orange")
# points(Ftemp_DAPC$Cumulative_position,Ftemp_DAPC$fst2010seg,type='p',pch=16,cex=1.5,col="green")
# points(All_methods$Cumulative_position,All_methods$fst2010seg,type='p',pch=17,cex=1.8,col="cyan")
# mtext("2010SEG",outer=FALSE,line=0.5,cex=2,at=2000)
# 
# #2010INT
# plot(qvals_summary$Cumulative_position,qvals_summary$fst2010int,ylim=c(-0.02,0.25),type="n",col="dimgray",ylab='',cex.axis=1.5,cex.lab=1.5)
# abline(h=0)
# rect(box_starts,-0.025,box_ends,0.25,col="gray91",border=NA)
# points(qvals_summary$Cumulative_position,qvals_summary$fst2010int,ylim=c(-0.02,0.25),type="p",col="dimgray",ylab='',pch=19,cex=0.8,cex.axis=1.5,cex.lab=1.5)
# points(Bayescan_outs$Cumulative_position,Bayescan_outs$fst2010int,type='p',pch=15,cex=1.2,col="blue")  #Bayescan
# points(Baye_Ftemp$Cumulative_position,Baye_Ftemp$fst2010int,type='p',pch=15,cex=1.2,col="blue")       #Also Bayescan only for INT line
# points(DAPC_structure$Cumulative_position,DAPC_structure$fst2010int,type='p',pch=17,cex=1.2,col="purple")  #DAPC
# points(Ftemp_DAPC$Cumulative_position,Ftemp_DAPC$fst2010int,type='p',pch=17,cex=1.2,col="purple")         #Also DAPC only for INT line
# points(All_methods$Cumulative_position,All_methods$fst2010int,type='p',pch=17,cex=1.5,col="brown4") ##Just Bayescan and DAPC for INT line
# rect(sw_2010int$Cum_start,0.245,sw_2010int$Cum_end,0.25,col='black',lwd=3)
# mtext("2010INT",outer=FALSE,line=0.5,cex=2,at=2000)



# sw_2002int <- sliding_window_outliers[sliding_window_outliers$Population=="2002INT",] #give me rows where Population column equals 2002INT but give me all columns
# sw_2002seg <- sliding_window_outliers[sliding_window_outliers$Population=="2002SEG",]
# sw_2006int <- sliding_window_outliers[sliding_window_outliers$Population=="2006INT",]
# sw_2006seg <- sliding_window_outliers[sliding_window_outliers$Population=="2006SEG",]
# sw_2010int <- sliding_window_outliers[sliding_window_outliers$Population=="2010INT",]
# sw_2010seg <- sliding_window_outliers[sliding_window_outliers$Population=="2010SEG",]
