### This is not a working R code
###    it is a set of examples of things I'm always looking up
###
###   Carolyn Tarpey | January 2018
### ---------------------------------------


#install.packages("vcfR")

library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)
library(hierfstat)
library(adegenet)
library(argparse)
library(stringi)


####<------------------------------------------------START HERE,
#####################################################################################################################################

####Write a table to a file
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Even_filtered_noMultAlleles_Haplo.txt", "wb")
write.table(Even_filteredHaplotypes,outputFile,quote=FALSE,row.names=TRUE,col.names=TRUE,eol="\n")
close(outputFile)

####Read a data file in
FST_file<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Genepop/HWE_Genepop_FST_R.txt", header=TRUE, 
                     stringsAsFactors = FALSE, na.strings = "-" )





####Data manipulation
Even_tag_pos <- loci_table_t[loci_table_t$Tag%in%just_singletons_tags,]

genotype_file_t<-gsub("X","",colnames(genotype_file))

allGenos_oneSNP_temp_tags<-data.frame(str_split_fixed(allGenos_oneSNP_temp$Locus,"_",2))

FST_file_ASIAeven_sort <- FST_ASIA_even[order(FST_ASIA_even$ASIA_even, decreasing=TRUE),]

### Set opperations

#get the unique tags 
singletons_tags<- unique(singletons$Tag)

#get the tags for singletons that are not in the paralog list
just_singletons_tags<- setdiff(singletons_tags, paralog_tags)




###Structure of for if/else loops
for (s in 1:length(Even_tags)) {
  tested_tag <- Even_tags[s]
  trial_set <- Even_SNPS[which(Even_SNPS$Tag %in% tested_tag),]
  if (any(trial_set$Btw_17_73 == FALSE)){
    Even_SNPs_InRange[[s]]<- "FALSE"
  } else {
    Even_SNPs_InRange[[s]]<- "TRUE"
  }
}




## convert the vector (or list of atomic vectors) to a Data Frame
Even_filt_snps <-data.frame(stri_list2matrix(Keep_even, byrow=TRUE))
colnames(Even_filt_snps)<-"Even_SNP_Index"

##make a dataframe with the individuals of one lineage
Even_inds<-popMapLineages[popMapLineages$Lineage=="Even",]

###Combine the output into a matrix that has the tag and the flag
Odd_matrix <-cbind(Odd_tags, Odd_SNPs_InRange)
colnames(Odd_matrix) <- c("Tags","SNPs_InRange")


##apply and Functions
locusGenoRate<-apply(allGenos,2,function(x) 1-(sum(x=="0000")/dim(allGenos)[1]))



###@@@@@@@@@@@@@@@@@@@@@@Plots for the Manuscript, have labels as Principal components
ALL_col<- c("#cd5490","#cd5490","#c9cf4a","#c9cf4a","#7297ee","#7297ee","#61005e","#61005e", "#677f3e","#677f3e","#00158a","#00158a","#e0957e","#e0957e")

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

dev.off()
