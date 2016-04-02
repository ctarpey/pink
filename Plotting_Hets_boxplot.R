### Boxplots of the heterozygosity of loci
###    Het from genepop file, outliers by arlequin, sig at 0.01%
### Carolyn Tarpey | April 2016
### ---------------------------------------


library(adegenet)
library(ggplot2)
library(ape)
library(reshape2)

library(devtools)
library(ggplot2)

give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}


setwd('G:/Analysis/Pop_analysis/Populations_b3_may/OutlierObservedHet')

#making dataframes of hets

odd_out <- read.table ("odd_het.txt", sep = "\t", header = TRUE)
odd_neut <- read.table ("odd_het_neut.txt", sep = "\t", header = TRUE)
odd_neut_nodup <- read.table ("odd_het_neut_nodup.txt", sep = "\t", header = TRUE)
even_out<- read.table ("even_het.txt", sep = "\t", header = TRUE)
even_neut <- read.table ("even_het_neut.txt", sep = "\t", header = TRUE)
even_neut_nodup <- read.table ("even_het_neut_nodup.txt", sep = "\t", header = TRUE)

# Basic Boxplot with all pops

pdf("Plots/odd_even_boxplots.pdf", width = 8, height = 5)

boxplot(odd_out, data=odd_out, main="Observed Heterozygosity of Odd Lineage Outliers",
        xlab="Odd Lineage Populations", ylab="Observed Heterozygosity")
boxplot(odd_neut, data=odd_neut, main="Observed Heterozygosity of Odd Lineage Neutral Loci",
        xlab="Odd Lineage Populations", ylab="Observed Heterozygosity")
boxplot(odd_neut_nodup, data=odd_neut_nodup, main="Observed Heterozygosity of Odd Lineage Neutral Loci, No Duplicates",
        xlab="Odd Lineage Populations", ylab="Observed Heterozygosity")
boxplot(even_out, data=even_out, main="Observed Heterozygosity of Even Lineage Outliers",
        xlab="Odd Lineage Populations", ylab="Observed Heterozygosity")
boxplot(even_neut, data=even_neut, main="Observed Heterozygosity of Even Lineage Neutral Loci",
        xlab="Odd Lineage Populations", ylab="Observed Heterozygosity")
boxplot(even_neut_nodup, data=even_neut_nodup, main="Observed Heterozygosity of Even Lineage Neutral Loci, No Duplicates",
        xlab="Odd Lineage Populations", ylab="Observed Heterozygosity")

dev.off()

#making dataframes of hets With Beringia

odd_out_beringia <- read.table ("odd_het_beringia.txt", sep = "\t", header = TRUE)
odd_neut_beringia <- read.table ("odd_het_neut_beringia.txt", sep = "\t", header = TRUE)
odd_neut_nodup_beringia <- read.table ("odd_het_neut_nodup_beringia.txt", sep = "\t", header = TRUE)
even_out_beringia<- read.table ("even_het_beringia.txt", sep = "\t", header = TRUE)
even_neut_beringia <- read.table ("even_het_neut_beringia.txt", sep = "\t", header = TRUE)
even_neut_nodup_beringia <- read.table ("even_het_neut_nodup_beringia.txt", sep = "\t", header = TRUE)

# Basic Boxplot with all pops

pdf("Plots/Beringia_boxplots.pdf", width = 8, height = 5)

boxplot(odd_out_beringia, data=odd_out_beringia, width <- length(odd_out_beringia), main="Observed Heterozygosity of Odd Lineage Outliers",
        xlab="Odd Lineage Populations", ylab="Observed Heterozygosity")
boxplot(odd_neut_beringia, data=odd_neut_beringia, main="Observed Heterozygosity of Odd Lineage Neutral Loci",
        xlab="Odd Lineage Populations", ylab="Observed Heterozygosity")
boxplot(odd_neut_nodup_beringia, data=odd_neut_nodup_beringia, main="Observed Heterozygosity of Odd Lineage Neutral Loci, No Duplicates",
        xlab="Odd Lineage Populations", ylab="Observed Heterozygosity")
boxplot(even_out_beringia, data=even_out_beringia, main="Observed Heterozygosity of Even Lineage Outliers",
        xlab="Odd Lineage Populations", ylab="Observed Heterozygosity")
boxplot(even_neut_beringia, data=even_neut_beringia, main="Observed Heterozygosity of Even Lineage Neutral Loci",
        xlab="Odd Lineage Populations", ylab="Observed Heterozygosity")
boxplot(even_neut_nodup_beringia, data=even_neut_nodup_beringia, main="Observed Heterozygosity of Even Lineage Neutral Loci, No Duplicates",
        xlab="Odd Lineage Populations", ylab="Observed Heterozygosity")

dev.off()


#making dataframes of hets Wityh Beringia

pdf("Plots/Beringia_all_boxplots.pdf", width = 8, height = 5)

ALL_out_beringia <- read.table ("ALL_het_beringia.txt", sep = "\t", header = TRUE)
ALL_neut_beringia <- read.table ("ALL_het_neut_beringia.txt", sep = "\t", header = TRUE)
ALL_neut_nodup_beringia <- read.table ("ALL_het_neut_nodup_beringia.txt", sep = "\t", header = TRUE)

dev.off()

# Basic Boxplot with all pops

boxplot(ALL_out_beringia, data=ALL_out_beringia, main="Observed Heterozygosity of Outliers",
        xlab="Odd Lineage Populations", ylab="Observed Heterozygosity")
boxplot(ALL_neut_beringia, data=ALL_neut_beringia, main="Observed Heterozygosity of Neutral Loci",
        xlab="Odd Lineage Populations", ylab="Observed Heterozygosity")
boxplot(ALL_neut_nodup_beringia, data=ALL_neut_nodup_beringia, main="Observed Heterozygosity of Neutral Loci, No Duplicates",
        xlab="Odd Lineage Populations", ylab="Observed Heterozygosity")