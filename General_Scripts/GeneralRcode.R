### This is not a working R code
###    it is a set of examples of things I'm always looking up
###
###   Carolyn Tarpey | Ongoing
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
library(reshape2)


####<------------------------------------------------START HERE,
#####################################################################################################################################

####Write a table to a file
outputFile <- file("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Even_filtered_noMultAlleles_Haplo.txt", "wb")
write.table(Even_filteredHaplotypes,outputFile,quote=FALSE,row.names=TRUE,col.names=TRUE,eol="\n")
close(outputFile)

####Read a data file in
FST_file<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/Genepop/HWE_Genepop_FST_R.txt", header=TRUE, 
                     stringsAsFactors = FALSE, na.strings = "-" )

#load genepop files Keeping the four digit format, 0303, 0000 etc
genos1<-read.table("Z:/WORK/TARPEY/Exp_Pink_Pops/STACKS/Ustacks/OriginalWhitelists/Whitelist_Genepops/whitelist_1.genepop",colClasses="factor")




####Data manipulation
Even_tag_pos <- loci_table_t[loci_table_t$Tag%in%just_singletons_tags,]

genotype_file_t<-gsub("X","",colnames(genotype_file))

allGenos_oneSNP_temp_tags<-data.frame(str_split_fixed(allGenos_oneSNP_temp$Locus,"_",2))

FST_file_ASIAeven_sort <- FST_ASIA_even[order(FST_ASIA_even$ASIA_even, decreasing=TRUE),]

#change the variable to the number you want 
NAeven_var <- 500 #<--------------------------------------------dictates the number of loci
NAeven_set <- FST_file_NAeven_sort_inRange[1:NAeven_var,]
NAeven_set_avgFST <- mean(NAeven_set$NA_even)
NAeven_set_maxFST <- max(NAeven_set$NA_even)
NAeven_set_minFST <- min(NAeven_set$NA_even)
cat("NAeven FST max, min, average: ",NAeven_set_maxFST, NAeven_set_minFST, NAeven_set_avgFST )

## Add these list of failed to the list of passed to get a list to pull the sequences to make a FASTA file 
Even_failed_passed_list <- append(NA_ASIA_even_union_FST_sort$Locus, Even_failed_union)

# merge the SNP_inRange designations
Even_FST_file <- merge(Even_FST_file, Even_Tag_Flags, by = "Tag")

#create new data frame from columns of another
FST_NA_even<- Even_FST_file[,c("Locus","NA_even","Even_SNPs_InRange")]

#add new empty column to the end of each data frame for the flags
Even_SNPS$Btw_17_73 <- NA

#delete a column
#add new empty column to the end of each data frame for the flags
Even_SNPS$Btw_17_73 <- NULL

#add comma to name if not already there!
EVEN_pops$Indv_name <- paste0(EVEN_pops$Indv_name, ",")

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





#<----------------------------------------PLOTTING

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



#histograms with bars next to eachother 

Histo_7_colors<- c("#6e016b","#88419d","#8c6bb1","#8c96c6","#7bccc4",'#4eb3d3','#2b8cbe')
Histo_14_colors<- c("#6e016b","#6e016b","#88419d","#88419d","#8c6bb1","#8c6bb1","#8c96c6","#8c96c6","#7bccc4","#7bccc4",'#4eb3d3','#4eb3d3','#2b8cbe','#2b8cbe')
E_O_Colors <- c("#00158a", "#c9cf4a")

even_snps_counts <- as.data.frame(table(Even_SNPS$Tag))
even_snps_counts <- as.data.frame(table(even_snps_counts$Freq))
colnames(even_snps_counts) <-c("SNPs","Even_Count")
even_snps_counts

odd_snps_counts <- as.data.frame(table(Odd_SNPS$Tag))
odd_snps_counts <- as.data.frame(table(odd_snps_counts$Freq))
colnames(odd_snps_counts) <-c("SNPs","Odd_Count")
odd_snps_counts

snp_counts<-merge(odd_snps_counts,even_snps_counts, by = "SNPs")
snp_counts_melt <- melt(snp_counts)

#pdf("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLOTs_Snps_per_Haplotype_lineage.pdf", width = 9, height = 7)

# histogram with the even and odd next to eachother in one plot 
ggplot(data=snp_counts_melt, aes(snp_counts_melt$SNPs, snp_counts_melt$value)) + theme_bw() +
  geom_bar(aes(fill= snp_counts_melt$variable), alpha = 0.5, position = "dodge", stat="identity") +
  xlab("Number of SNPs per Haplotype") + ylab("Count") + ggtitle("Variable SNPs per Haplotype ") +
  scale_fill_manual(values = E_O_Colors, labels = c("Odd", "Even")) + labs(fill="Lineage")

#dev.off()  

