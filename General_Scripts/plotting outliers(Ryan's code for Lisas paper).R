library(VennDiagram)
library(ggplot2)
library(reshape2)
library(stringr)

# read in assembled PE contigs
pink_8184_FASTA <- read.table('/home/ipseg/Data/Pink/8036/FASTA/8184_PE_contigs.tab',
                              sep = '\t',
                              header = FALSE,
                              stringsAsFactors = FALSE)
snpID <- as.numeric(str_split_fixed(as.character(pink_8184_FASTA$V1),pattern = '[:_]', 10)[,9])
queryID <- str_split_fixed(as.character(pink_8184_FASTA$V1),pattern = '>', 10)[,2]
pink_8184_FASTA$V1 = queryID
pink_8184_FASTA = cbind(snpID, pink_8184_FASTA)
names(pink_8184_FASTA) <- c('snpID', "queryId", "contig")

# 8036 subset and BLAST results
subset_8036 <- read.table ("/home/ipseg/Data/Pink/8036/8036.subset", stringsAsFactors = FALSE)
names(subset_8036) <- 'snpID'

blast_8184  <- read.table (file = '/home/ipseg/Data/Pink/8036/BLAST/8184_cap3_contigs_BLASTX.tab', header = FALSE, stringsAsFactors = FALSE)
names(blast_8184) <- c("queryId", "subjectId", "percIdentity", "alnLength", "mismatchCount", "gapOpenCount", "queryStart", "queryEnd", "subjectStart", "subjectEnd", "eVal", "bitScore") 
blast_8184$queryId <- str_split_fixed(blast_8184$queryId, pattern = 'allele', n = 3)[,1]

# subset based on 8036
contigs_8036 <- pink_8184_FASTA[pink_8184_FASTA$snpID %in% subset_8036[[1]], ]
snpID_blast <- str_split_fixed(as.character(blast_8184$queryId), pattern = '[:_]', 10)[,9]
which_in_8036 <- snpID_blast %in% as.character(subset_8036[[1]])
blast_8036  <- blast_8184[which_in_8036, ]

snpID <- str_split_fixed(as.character(blast_8036$queryId), pattern = '_', 11)[,9]
snpID <- as.numeric(str_split_fixed(as.character(snpID), pattern = ":", 2)[,1])
# len values seem wrong
len = as.numeric(str_split_fixed(as.character(blast_8036$queryId), pattern = '_', 11)[,5])
blast_8036 = cbind(snpID, len, blast_8036)

# select based on eVal <.0001
sig_blast_8036 <- blast_8036[blast_8036$eVal < 0.0001,]
seq_blast_8036 <- merge(blast_8036, contigs_8036, by = c('snpID', 'queryId'), all = TRUE)
sig_seq_blast_8036 <- merge(sig_blast_8036, contigs_8036, by = c('snpID', 'queryId'))
sorted_seq_blast_8036 = seq_blast_8036[order(seq_blast_8036$snpID, seq_blast_8036$eVal), ]
bs <- sorted_seq_blast_8036[match(x = unique(sorted_seq_blast_8036$snpID), table = sorted_seq_blast_8036$snpID), ]

# read 8036 consensus sequences
consensus_seq_8036 <- read.table(file = '/home/ipseg/Data/Pink/8036/8036_consensus_seq.txt', header = FALSE, sep = '\t')
names(consensus_seq_8036) <- c('snpID', 'consensus')

bs2 <-merge(consensus_seq_8036, bs, all = TRUE)

#write.table (blast_8036, file = "C:/Users/Lisa Seeb/Dropbox/R/BLAST Parse/blast_8036.tab", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")
#write.table (sig_blast_8036, file = "C:/Users/Lisa Seeb/Dropbox/R/BLAST Parse/sig_8036.tab", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

pink_8036_outliers <- read.table(#"Z:/WORK/PINK RAD/OUTLIERS/8036_outliers.txt",
                                "/home/ipseg/Data/Pink/8036/8036_outliers.txt",
                                 sep = "\t",
                                 header = TRUE)

bayescan_both <- pink_8036_outliers[pink_8036_outliers$Bayescan_EVEN_ODD == "both", ]$SNP_ID
bayescan_odd  <- pink_8036_outliers[pink_8036_outliers$Bayescan_EVEN_ODD %in% c("odd", "both"), ]$SNP_ID
bayescan_even <- pink_8036_outliers[pink_8036_outliers$Bayescan_EVEN_ODD %in% c("even", "both"), ]$SNP_ID
bayescan_all  <- pink_8036_outliers[pink_8036_outliers$Bayscan_ALL == "yes", ]$SNP_ID
just_bayescan_odd       <- setdiff(bayescan_odd, bayescan_even)
just_bayescan_even      <- setdiff(bayescan_even, bayescan_odd )

bayenv_snoh   <- pink_8036_outliers[pink_8036_outliers$Bayenv_SNOH == "yes", ]$SNP_ID
bayenv_range  <- pink_8036_outliers[pink_8036_outliers$Bayenv_RANGE == "yes", ]$SNP_ID
intersect_bayenv  <- intersect(bayenv_snoh, bayenv_range)
union_bayenv  <- union(bayenv_snoh, bayenv_range)

PCA_both      <- pink_8036_outliers[pink_8036_outliers$PCA == "both", ]$SNP_ID
PCA_odd       <- pink_8036_outliers[pink_8036_outliers$PCA %in% c("odd", "both"), ]$SNP_ID
PCA_even      <- pink_8036_outliers[pink_8036_outliers$PCA %in% c("even", "both"), ]$SNP_ID

bayescan_bayenv <- intersect(bayescan_both, union_bayenv)
PCA_bayenv <- intersect(union_bayenv, PCA_both)
PCA_bayescan <- intersect(bayescan_both, PCA_both)
PCA_bayescan_bayenv <- intersect(PCA_both, bayescan_bayenv)
double_outliers <- intersect(union(bayescan_odd, bayescan_even), union_bayenv)

# construct outlier data.frames
BE_snoh <- data.frame(locus = pink_8036_outliers[pink_8036_outliers$SNP_ID %in% bayenv_snoh, ]$locus, bayenv_snoh = 'yes')
BE_range <- data.frame(locus = pink_8036_outliers[pink_8036_outliers$SNP_ID %in% bayenv_range, ]$locus, bayenv_range = 'yes')

BS_all <- data.frame(locus = pink_8036_outliers[pink_8036_outliers$SNP_ID %in% bayescan_all, ]$locus, bayescan_all = 'yes')
BS_even <- data.frame(locus = pink_8036_outliers[pink_8036_outliers$SNP_ID %in% bayescan_even, ]$locus, bayescan_even = 'yes')
BS_odd <- data.frame(locus = pink_8036_outliers[pink_8036_outliers$SNP_ID %in% bayescan_odd, ]$locus, bayescan_odd = 'yes')

PC_odd <- data.frame(locus = pink_8036_outliers[pink_8036_outliers$SNP_ID %in% PCA_odd, ]$locus, PCA_odd = 'yes')
PC_even <- data.frame(locus = pink_8036_outliers[pink_8036_outliers$SNP_ID %in% PCA_even, ]$locus, PCA_even = 'yes')
#PC_all <- data.frame(locus = pink_8036_outliers[pink_8036_outliers$SNP_ID %in% PCA_both, ]$locus, PCA_both = 'yes')

outliers_merge <- merge(
                    merge(
                      merge(BE_snoh, BE_range, all = TRUE),
                      merge(PC_odd, PC_even, all = TRUE), all = TRUE),
                    merge(BS_all,
                      merge(BS_even, BS_odd, all = TRUE), all = TRUE), all = TRUE)

outliers_merge = data.frame(lapply(outliers_merge, as.character), stringsAsFactors=FALSE)
outliers_merge[is.na(outliers_merge)] <- 'no'

outliers_merge$locus <- as.character(str_split_fixed(outliers_merge$locus, pattern = 'SNP_', n=2)[,2])
names(outliers_merge) <- c('snpID', names(outliers_merge)[2:8])
outliers_full <- merge(outliers_merge, subset_8036, all  = TRUE)
outliers_full[is.na(outliers_full)] <- 'no'
bs3 <- merge(bs2, outliers_full, all = TRUE)

# read basic data on genotypes, and seqeunce variation
basic_70133 <- read.table('/home/ipseg/Data/Pink/8036/71033_basic_data.txt', header = TRUE, sep = "\t")
bs4 <- merge(bs3, basic_70133, all.x = TRUE)
# fields to remove:
# queryID, len, consensus
bs4 <- bs4[,  !(colnames(bs4) %in% c("queryId","len","consensus"))]

write.table(bs4, file = '/home/ipseg/Data/Pink/8036/supp_table_1.txt', 
            sep = "\t", col.names = TRUE, row.names= FALSE, quote = FALSE)

plot_subset_1 = unique(c(bayescan_bayenv, PCA_bayenv, PCA_bayescan, PCA_bayescan_bayenv))
plot_subset_2 = unique(c(PCA_both, union_bayenv, bayescan_both))
 
venn.diagram(x = list("ODD" = bayescan_odd, "EVEN" = bayescan_even), 
             filename = 'bayescan_EVEN_ODD.tiff',
             main = 'Bayescan Outliers',
             main.cex = 3,
             fill = c('blue', 'red'))

venn.diagram(x = list("Bayescan" = union(bayescan_odd, bayescan_even), "Bayenv" = union_bayenv), 
             filename = 'Broad bayescan_bayenv.tiff',
             main = 'Broad \nBayescan Bayenv Overlap',
             main.cex = 3,
             fill = c('purple', 'orange'),
             ext.text = FALSE,
             scaled = TRUE
             )

venn.diagram(x = list("Bayescan" = intersect(bayescan_odd, bayescan_even), "Bayenv" = union_bayenv), 
             filename = 'Narrow bayescan_bayenv.tiff',
             main = 'Narrow \nBayescan Bayenv Overlap',
             main.cex = 3,
             fill = c('purple', 'orange'),
             ext.text = FALSE
)

to_plot_5 <- melt(pink_8036_outliers[pink_8036_outliers$SNP_ID %in% plot_subset_1, ][2:8])
to_plot_14 <- melt(pink_8036_outliers[pink_8036_outliers$SNP_ID %in% plot_subset_2, ][2:8])
names(to_plot_5) <- c("locus", "pop", "freq")
names(to_plot_14) <- c("locus", "pop", "freq")

lineage_5 = rep(rep( c('odd', 'even'), each = 5), 3)
lineage_14 = rep(rep( c('odd', 'even'), each = 14), 3)

to_plot_5 <- cbind(to_plot_5, lineage_5)
to_plot_14 <- cbind(to_plot_14, lineage_14)

to_plot_5$pop <- rep(c('PS', 'PWS', 'NS'), each = 10)
to_plot_14$pop <- rep(c('PS', 'PWS', 'NS'), each = 28) 

to_plot_5$pop <- factor(to_plot_5$pop, levels = c('NS', 'PWS', 'PS'))
to_plot_14$pop <- factor(to_plot_14$pop, levels = c('NS', 'PWS', 'PS'))

# write.table(to_plot, 
#             file = "Z:/WORK/PINK RAD/OUTLIERS/6_outliers.tsv", 
#             quote = FALSE, 
#             sep = "\t",
#             row.names= FALSE)


my_plot_5 <- ggplot(to_plot_5, aes(y = freq, x = pop, color = lineage_5, fill = lineage_5))
my_plot_5 +
  facet_wrap(facets = ~ locus, nrow = 2, ncol = 3, scales="free_x") +
  geom_bar(stat="identity", width = .7, position = 'dodge') +
  theme_classic() +
  scale_color_manual(values = c('black', 'black')) +
  scale_fill_manual(values = c('gray', 'black')) +
  theme(axis.text.x = element_text(angle = 45, vjust = .9, hjust = 1, size = 14, face = 'bold'),
        axis.text.y = element_text(size = 14, face = 'bold'),
        axis.title = element_text(size = 16, face = 'bold'),
        legend.position = c(.85,.2)) + 
  labs(y = 'Allele Frequency', x = 'Site')
  
rect_left <- factor('PS')
bg_rectangles_5 <- data.frame(xmin = rect_left, xmax = rect_left   , ymin = 0, ymax = 1)

ggplot() +
  geom_rect(data = bg_rectangles_5,  aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='gray80', color='gray80') +
  geom_bar(data = to_plot_5, aes(y = freq, x = pop, color = lineage_5, fill = lineage_5), stat="identity", width = .7, position = 'dodge') +
  facet_wrap(facets = ~ locus, nrow = 2, ncol = 3, scales="free_x") +
  theme_classic() +
  scale_color_manual(values = c('black', 'black')) +
  scale_fill_manual(values = c('gray', 'black')) +
  theme(axis.text.x = element_text(angle = 45, vjust = .9, hjust = 1, size = 14, face = 'bold'),
        axis.text.y = element_text(size = 14, face = 'bold'),
        axis.title = element_text(size = 16, face = 'bold')) + 
  labs(y = 'Allele Frequency', x = 'Site') 


# alternating gray scale on background
# rename populations

my_plot_14 <- ggplot(to_plot_14, aes(y = freq, x = pop, color = lineage_14, fill = lineage_14))
my_plot_14 +
  facet_wrap(facets = ~ locus, nrow = 4, ncol = 4, scales="free_x") +
  geom_bar(stat="identity", width = .6, position = 'dodge') +
  theme_classic() +
  scale_color_manual(values = c('black', 'black')) +
  scale_fill_manual(values = c('gray', 'black')) +
  theme(axis.text.x = element_text(angle = 45, vjust = .9, hjust = 1, size = 16))


#####
rect_left = 0.5
bg_rect = data.frame(xmin= rect_left, xmax= rect_left+1, ymin=0, ymax=1)
geom_rect(data = bg_rectaes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),fill='gray80', alpha=0.8)



