### Manhattan plots of Outliers from LFMM and Arlequin
###Carolyn Tarpey | April 2016
### ---------------------------------------


setwd('G:/Analysis/Mapping/Outliers_Manhattan_Plot/Arlequin_LFMM')


outliers = read.table('Outliers_Map_self.txt', header = TRUE)
as.data.frame(outliers)

head(outliers)
table(outliers$CHR)

#SNP CHR   BP BP_seq Paralog    o_p_PT1   o_p_PT2    o_p_LO    o_p_LA    e_p_PT1    e_p_PT2    e_p_LO   
#e_p_LA o_p_Fst_arl e_p_Fst_arl odd_PT1 odd_PT2 odd_LO odd_LA
#odd_any even_PT1 even_PT2 even_LO even_LA even_any odd_0.01 even_0.01 overlap_0.01 overlap_PT1 overlap_PT2 overlap_LO overlap_LA


#to make the tick marks in the middle only has to be done once
chrNum=26
bpMidVec <- vector(length=chrNum)
for (i in 1:chrNum){ndx <- which(outliers[, 2]==i)
posSub <- outliers[ndx, 4]
bpMidVec[i] <- ((max(posSub) - min(posSub))/2) + min(posSub)
}
bpMidVec


#This is super basic, no outliers colored,outliers determine the alpha.  requires you figure out where the -log of p values makes outliers 
ggplot(outliers) + geom_point(aes(x=BP_seq, y=o_p_PT1, colour=as.factor(CHR), alpha=odd_PT1), size = 2) + 
  scale_color_manual(values=rep(c('black', 'dark grey'), 13)) + theme_bw(base_size=24) + theme(legend.position='none') + 
  scale_x_continuous(labels=as.character(1:26), breaks=bpMidVec) +  scale_y_continuous(limits = c(0, 20)) +
  ggtitle('Odd Lineage Outliers PT1') + xlab('Linkage Group') +  ylab('-log10(P)')  
  geom_hline(aes(yintercept= 2.40945563289108), linetype=1, col='red', lwd=1)


#pdf("plots/LFMMresults_Manhattan_plots.pdf", width = 8, height = 5)


#	odd_0.01	even_0.01	overlap_0.01	overlap_PT1	overlap_PT2	overlap_LO	overlap_LA

###############################################
#LFMM
# odd pvalues: o_p_PT1	o_p_PT2	o_p_LO	o_p_LA	
# outlisers: odd_PT1	odd_PT2	odd_LO	odd_LA	odd_any	

ggplot(outliers) + geom_point(aes(x=BP_seq, y=o_p_PT1, colour=as.factor(CHR), alpha=odd_PT1), size = 2) + 
  scale_color_manual(values=rep(c('black', 'dark grey'), 13)) + theme_bw(base_size=24) + theme(legend.position='none') + 
  scale_x_continuous(labels=as.character(1:26), breaks=bpMidVec) +  scale_y_continuous(limits = c(0, 20)) +
  ggtitle('Odd Lineage Outliers PT1') + xlab('Linkage Group') +  ylab('-log10(P)')  +
  geom_hline(aes(yintercept= 2.40945563289108), linetype=1, col='red', lwd=1)

ggplot(outliers) + geom_point(aes(x=BP_seq, y=o_p_PT2, colour=as.factor(CHR), alpha=odd_PT2), size = 2) + 
  scale_color_manual(values=rep(c('black', 'dark grey'), 13)) + theme_bw(base_size=24) + theme(legend.position='none') + 
  scale_x_continuous(labels=as.character(1:26), breaks=bpMidVec) +  scale_y_continuous(limits = c(0, 23)) +
  ggtitle('Odd Lineage Outliers PT2') + xlab('Linkage Group') +  ylab('-log10(P)')  +
  geom_hline(aes(yintercept= 2.47715752), linetype=1, col='red', lwd=1)

ggplot(outliers) + geom_point(aes(x=BP_seq, y=o_p_LO, colour=as.factor(CHR), alpha=odd_LO), size = 2) + 
  scale_color_manual(values=rep(c('black', 'dark grey'), 13)) + theme_bw(base_size=24) + theme(legend.position='none') + 
  scale_x_continuous(labels=as.character(1:26), breaks=bpMidVec) +  scale_y_continuous(limits = c(0, 21)) +
  ggtitle('Odd Lineage Outliers LONG') + xlab('Linkage Group') +  ylab('-log10(P)')  +
  geom_hline(aes(yintercept= 2.372077804), linetype=1, col='red', lwd=1)

ggplot(outliers) + geom_point(aes(x=BP_seq, y=o_p_LA, colour=as.factor(CHR), alpha=odd_LA), size = 2) + 
  scale_color_manual(values=rep(c('black', 'dark grey'), 13)) + theme_bw(base_size=24) + theme(legend.position='none') + 
  scale_x_continuous(labels=as.character(1:26), breaks=bpMidVec) +  scale_y_continuous(limits = c(0, 20)) +
  ggtitle('Odd Lineage Outliers LAT') + xlab('Linkage Group') +  ylab('-log10(P)')  +
  geom_hline(aes(yintercept= 2.446908887), linetype=1, col='red', lwd=1)

#even pvalues: e_p_PT1	e_p_PT2	e_p_Le	e_p_LA	
# outliers: even_PT1	even_PT2	even_LO	even_LA	even_any


ggplot(outliers) + geom_point(aes(x=BP_seq, y=e_p_PT1, colour=as.factor(CHR), alpha=even_PT1), size = 2) + 
  scale_color_manual(values=rep(c('black', 'dark grey'), 13)) + theme_bw(base_size=24) + theme(legend.position='none') + 
  scale_x_continuous(labels=as.character(1:26), breaks=bpMidVec) +  scale_y_continuous(limits = c(0, 17)) +
  ggtitle('Even Lineage Outliers PT1') + xlab('Linkage Group') +  ylab('-log10(P)')  +
  geom_hline(aes(yintercept= 2.553768935), linetype=1, col='red', lwd=1)

ggplot(outliers) + geom_point(aes(x=BP_seq, y=e_p_PT2, colour=as.factor(CHR), alpha=even_PT2), size = 2) + 
  scale_color_manual(values=rep(c('black', 'dark grey'), 13)) + theme_bw(base_size=24) + theme(legend.position='none') + 
  scale_x_continuous(labels=as.character(1:26), breaks=bpMidVec) +  scale_y_continuous(limits = c(0, 15)) +
  ggtitle('Even Lineage Outliers PT2') + xlab('Linkage Group') +  ylab('-log10(P)')  +
  geom_hline(aes(yintercept= 2.704105379), linetype=1, col='red', lwd=1)

ggplot(outliers) + geom_point(aes(x=BP_seq, y=e_p_LO, colour=as.factor(CHR), alpha=even_LO), size = 2) + 
  scale_color_manual(values=rep(c('black', 'dark grey'), 13)) + theme_bw(base_size=24) + theme(legend.position='none') + 
  scale_x_continuous(labels=as.character(1:26), breaks=bpMidVec) +  scale_y_continuous(limits = c(0, 22)) +
  ggtitle('Even Lineage Outliers LONG') + xlab('Linkage Group') +  ylab('-log10(P)')  +
  geom_hline(aes(yintercept= 2.642001537), linetype=1, col='red', lwd=1)

ggplot(outliers) + geom_point(aes(x=BP_seq, y=e_p_LA, colour=as.factor(CHR), alpha=even_LA), size = 2) + 
  scale_color_manual(values=rep(c('black', 'dark grey'), 13)) + theme_bw(base_size=24) + theme(legend.position='none') + 
  scale_x_continuous(labels=as.character(1:26), breaks=bpMidVec) +  scale_y_continuous(limits = c(0, 15)) +
  ggtitle('Even Lineage Outliers LAT') + xlab('Linkage Group') +  ylab('-log10(P)')  +
  geom_hline(aes(yintercept= 2.590279851), linetype=1, col='red', lwd=1)


#overlap
#pvalues:
# outliers: overlap_PT1	overlap_PT2	overlap_LO	overlap_LA

#plot = manhattan(outliers, p = "LFMM_ASE_K2_LO", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASE_K2_LO", col = c("blue4", "orange3"))

#odd vs even
#pvalues:
# outliers: 
#plot = manhattan(outliers, p = "LFMM_ASE_K2_LO", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASE_K2_LO", col = c("blue4", "orange3"))


###############################################
#Arlequin

#pvalues: o_p_Fst_arl	e_p_Fst_arl	
#ouliers: odd_0.01	even_0.01

ggplot(outliers) + geom_point(aes(x=BP_seq, y=o_p_Fst_arl, colour=as.factor(CHR), alpha=odd_0.01), size = 2) + 
  scale_color_manual(values=rep(c('black', 'dark grey'), 13)) + theme_bw(base_size=24) + theme(legend.position='none') + 
  scale_x_continuous(labels=as.character(1:26), breaks=bpMidVec) +  scale_y_continuous(limits = c(0.3, 10)) +
  ggtitle('Arlequin Odd Lineage Outliers') + xlab('Linkage Group') +  ylab('-log10(P)')  +
  geom_hline(aes(yintercept= 2.000284122), linetype=1, col='red', lwd=1)

ggplot(outliers) + geom_point(aes(x=BP_seq, y=e_p_Fst_arl, colour=as.factor(CHR), alpha=even_0.01), size = 2) + 
  scale_color_manual(values=rep(c('black', 'dark grey'), 13)) + theme_bw(base_size=24) + theme(legend.position='none') + 
  scale_x_continuous(labels=as.character(1:26), breaks=bpMidVec) +  scale_y_continuous(limits = c(0.3, 10)) +
  ggtitle('Arlequin Even Lineage Outliers') + xlab('Linkage Group') +  ylab('-log10(P)')  +
  geom_hline(aes(yintercept= 2.002051232), linetype=1, col='red', lwd=1)

#overlap
#pvalues:
# outliers: 
#plot = manhattan(outliers, p = "LFMM_ASE_K2_LO", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASE_K2_LO", col = c("blue4", "orange3"))

#odd vs even
#pvalues:
# outliers: 

#plot = manhattan(outliers, p = "LFMM_ASE_K2_LA", chr = "CHR", bp = "BP", snp = "SNP", main= "LFMM_ASE_K2_LA", col = c("blue4", "orange3"))


#dev.off()


################
#This is super basic, no outliers colored, requires you figure out where the -log of p values makes outliers 
# ggplot(outliers) + geom_point(aes(x=BP_seq, y=o_p_PT1, colour=as.factor(CHR)), size = 2) + 
#   scale_color_manual(values=rep(c('black', 'dark grey'), 13)) + theme_bw(base_size=24) + theme(legend.position='none') + 
#   scale_x_continuous(labels=as.character(1:26), breaks=bpMidVec) +  scale_y_continuous(limits = c(0, 20)) +
#   ggtitle('Odd Lineage Outliers PT1') + xlab('Linkage Group') +  ylab('-log10(P)') + 
#   geom_hline(aes(yintercept= 2.40945563289108), linetype=1, col='red', lwd=1)

# 
# ggplot(outliers) + geom_point(aes(x=BP_seq, y=o_p_PT1, size = 0.53, colour=as.factor(o_p_PT1)) +
# scale_color_manual(values=rep(c('black', "grey"), 13)) + theme_bw(base_size=24) + theme(legend.position='none') + 
# scale_x_continuous(labels=as.character(1:26), breaks=bpMidVec) +  scale_y_continuous(limits = c(0, 20)) +
# ggtitle('Odd Lineage Outliers PT1') + xlab('Linkage Group') +  ylab('-log10(P)') 
#                               
# ggplot(outliers) + geom_point(aes(x=BP_seq, y=o_p_PT1, colour=as.factor(odd_PT1))) + 
# theme_bw(base_size=24) + theme(legend.position='none') + 
# scale_x_continuous(labels=as.character(1:26), breaks=bpMidVec) +  scale_y_continuous(limits = c(0, 20)) +
# ggtitle('Odd Lineage Outliers PT1') + xlab('Linkage Group') +  ylab('-log10(P)') 
#                               
# ggplot(outliers) + geom_point(aes(x=BP_seq, y=o_p_PT1, colour=as.factor(CHR))) +  
# scale_color_manual(values=rep(c('black', 'dark grey'), 13)) + geom_point(aes(x=BP_seq, y=o_p_PT1, colour=as.factor(outliers$odd_PT1), colour= c("black","pink"))) + 
# theme_bw(base_size=24) + theme(legend.position='none') + 
# scale_x_continuous(labels=as.character(1:26), breaks=bpMidVec) +  scale_y_continuous(limits = c(0, 20)) +
# ggtitle('Odd Lineage Outliers PT1') + xlab('Linkage Group') +  ylab('-log10(P)')
#                               
# 
