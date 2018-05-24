### PLINK associtation test clump analysis visualization
###    for the pink panel of 2312 markers
###   
### Carolyn Tarpey | May 2018
# Updated may 2018 to make sure that all the populations were in their correct groups- Nome in Asia and lakel07 even etc.
### ---------------------------------------

## This requires that association tests in PLINK be run for each of the panels
## Tested the PLINK association testing with 4 categories,2 groups for each lineage: EA_case, ENA_case, OA_case, ONA_case, 
## case is the group in the lineage that was the phenotype present, so in that test we were testing which loci were associated with that region. 
## 

##These are example of the commands for plink to get the output we use here: 
# ODD_ASIA_CASE
# plink --ped Odd_panel_1197_PLINK_POS.ped --map Odd_panel_1197_PLINK_POS.map --chr-set 34 no-xy no-mt --allow-no-sex --keep-fam ODD_pops.txt 
# --assoc fisher --pheno Phenotypes.txt --mpheno 2 --out ODD_1197_PLINK_assoc_out_OA_case
# 
# plink  --ped Odd_panel_1197_PLINK_POS.ped --map Odd_panel_1197_PLINK_POS.map --chr-set 34 no-xy no-mt --allow-no-sex --keep-fam  ODD_pops.txt 
# --pheno Phenotypes.txt --mpheno 2 --clump ODD_1197_PLINK_assoc_out_OA_case.assoc.fisher --clump-best --out ODD_1197_PLINK_clump_out_OA_case

#Pink data filtering
library(ggplot2)
library(vcfR)
library(stringr)
library(dplyr)
library(gdata)

#load the PLINK association results for each of the groups: 
EA_case_clump <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/AssociationTesting/EVEN_1461_PLINK_clump_out_EA_case.clumped.best",sep="", header = TRUE)
EA_case_assoc <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/AssociationTesting/EVEN_1461_PLINK_assoc_out_EA_case.assoc.fisher",sep="", header = TRUE)
ENA_case_clump <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/AssociationTesting/EVEN_1461_PLINK_clump_out_ENA_case.clumped.best",sep="", header = TRUE)
ENA_case_assoc <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/AssociationTesting/EVEN_1461_PLINK_assoc_out_ENA_case.assoc.fisher",sep="", header = TRUE)
OA_case_clump <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/AssociationTesting/ODD_1197_PLINK_clump_out_OA_case.clumped.best",sep="", header = TRUE)
OA_case_assoc <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/AssociationTesting/ODD_1197_PLINK_assoc_out_OA_case.assoc.fisher",sep="", header = TRUE)
ONA_case_clump <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/AssociationTesting/ODD_1197_PLINK_clump_out_ONA_case.clumped.best",sep="", header = TRUE)
ONA_case_assoc <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/AssociationTesting/ODD_1197_PLINK_assoc_out_ONA_case.assoc.fisher",sep="", header = TRUE)

#Even and Odd results 
EO_clump <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/AssociationTesting/Panel_2312_PLINK_clump_EO.clumped.best",sep="", header = TRUE)
EO_assoc <-read.delim("Z:/WORK/TARPEY/Exp_Pink_Pops/Analysis/MarkerSelection/PLINK/AssociationTesting/Panel_2312_PLINK_assoc_EO.assoc.fisher",sep="", header = TRUE)









