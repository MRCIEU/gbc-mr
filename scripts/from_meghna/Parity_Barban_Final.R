## MR Analysis: Parity (Barban) ##
## European Population, n = 329,139 ##
# Meghna Sequeira #
# 2026 #

# Set working directory
setwd("/Users/meghnasequeira/MR_Data/Barban_Parity")

# Run libraries
library(data.table)
library(R.utils)
library(TwoSampleMR)
library(devtools)

######################################################################################################################################################################

#TSV to CSV

# Read the TSV
data <- read.delim("gwas-association-downloaded_2026-07-17-pubmedId_27798627.tsv", header = TRUE)

# Save as CSV
write.csv(data, "par-27798627.csv", row.names = FALSE)

# View CSV
#library(readr)
#data <- read_csv("~/Barban_Parity/par-27798627.csv")
#View(data)

# Supplementary information
#library(readxl)
#Supplementary_Tables <- read_excel("~/Barban_Parity/NIHMS841435-supplement-Supplementary_Tables.xls")
#View(Supplementary_Tables)

# Manually calculate SE in Excel

######################################################################################################################################################################

#Exposure preparation, outcome preparation, and MR

data_p <- fread("1_exp_par(in).csv")
dim(data_p)
head(data_p)

#Create exposure data=====================================================================================

#exposure_data_p <- read_exposure_data(
#  filename = "1_exp_par(in).csv",
#  clump = TRUE,
#  sep = ",",
#  phenotype_col = "PHENOTYPE", #Optional column name for the column with phenotype name corresponding to the SNP
#  snp_col = "SNPS",
#  beta_col = "BETA",
#  se_col = "SE",
#  eaf_col = "RISK.ALLELE.FREQUENCY",  
#  effect_allele_col = "STRONGEST.SNP.RISK.ALLELE",  
#  other_allele_col = "OTHER.ALLELE",  
#  pval_col = "P.VALUE",  
#  min_pval = 1e-200,  
#  log_pval = FALSE,  
#  chr_col = "CHR_ID",  
#  pos_col = "CHR_POS", 
#  clump_kb = 10000,  
# clump_r2 = 0.001,  
#  clump_p1 = 1,
#  bfile = "/Users/meghnasequeira/Kentistou_Age_at_Menarche/1000_genome_SAS/1000G_SAS",  
#  plink_bin = "/usr/local/bin/plink"
#)

# Clumping is unnecessary because all 3 SNPs are on different chromosomes; read file and proceed with harmonization
exposure_data_p <- read_exposure_data(
  filename = "1_exp_par(in).csv",
  clump = FALSE,
  sep = ",",
  phenotype_col = "PHENOTYPE",
  snp_col = "SNPS",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "RISK.ALLELE.FREQUENCY",
  effect_allele_col = "STRONGEST.SNP.RISK.ALLELE",
  other_allele_col = "OTHER.ALLELE",
  pval_col = "P.VALUE",
  chr_col = "CHR_ID",
  pos_col = "CHR_POS"
)
  
write.csv(exposure_data_p,"1_exposure_data.csv",row.names = FALSE)

#==========================================================================================================
data_2_p <- fread("3_G1_O_48_BeagleImp_Study_MichiHLA_Merged_FullyCleaned_PostImpQC_maf-geno-0.005_Logistic_Age_Gen_10PCs_OR.assoc.logistic.csv")
dim(data_2_p)
head(data_2_p)

#Outcome data==============================================================================================

outcome_data <- read_outcome_data(
  filename = "3_G1_O_48_BeagleImp_Study_MichiHLA_Merged_FullyCleaned_PostImpQC_maf-geno-0.005_Logistic_Age_Gen_10PCs_OR.assoc.logistic.csv",
  phenotype_col = "GBC",
  snps = NULL,
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1_MinorAlleleFreq",
  effect_allele_col = "A1_MinorAllele",
  other_allele_col = "A2_MajorAllele",
  pval_col = "P",
  min_pval = 1e-300,
  log_pval = FALSE,
  chr_col = "CHR",
  pos_col = "BP",
  #samplesize_col = "NMISS"
)

dim(outcome_data)

write.csv(outcome_data,"1_outcome_data.csv",row.names = FALSE)

#Harmonizaton===================================================================================================

outcome_data_p <- read_csv("1_outcome_data.csv")
dim(outcome_data_p)

dat1_p <- harmonise_data(exposure_data, outcome_data, action = 3)
write.csv(dat1_p,"harmonize_action_3.csv ", row.names = FALSE)
dim(dat1_p)

#Mr analysis=====================================================================================================

mr_p <-mr(dat1_p)
write.csv(mr_p," Final_MR_OUTPUT.csv", row.names = FALSE)

mr_bar_par <- read_csv(" Final_MR_OUTPUT.csv")
dat_bar_par <- read_csv("harmonize_action_3.csv ")

#Scatter plot of MR
p1 <- mr_scatter_plot(mr_bar_par, dat_bar_par)
p1[[1]] # Save as Barban_parity_Scatterplot.png

#Odds ratios
generate_odds_ratios(mr_bar_par)
write.csv(generate_odds_ratios(mr_bar_par), file = "generate_odds_ratios.csv", row.names=T)

#MR Heterogeneity
mr_heterogeneity(dat_bar_par)
write.csv(mr_heterogeneity(dat_bar_par), file = "mr_heterogeneity.csv", row.names=T)

#MR Forest plots (singlesnp)
res_single <- mr_singlesnp(dat_bar_par)
write.csv(mr_singlesnp(dat_bar_par), file = "res_single.csv", row.names=T)
p2 <- mr_forest_plot(res_single)
p2[[1]] # Save as Barban_parity_Forestplot1.png

#res_single2 <- mr_singlesnp(dat_bar_par, all_method=c("mr_ivw", "mr_two_sample_ml"))
#p2 <- mr_forest_plot(res_single2) 
#p2[[1]] # Save as Barban_parity_Forestplot2.png

#MR Leave one out
res_loo <- mr_leaveoneout(dat_bar_par)
write.csv(mr_leaveoneout(dat_bar_par), file = "mr_leaveoneout.csv", row.names=T)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]] # Save as Barban_parity_Loo.png

#MR Pleiotropy test
mr_pleiotropy_test(dat_bar_par)
write.csv(mr_pleiotropy_test(dat_bar_par), file = "mr_pleiotropy_test.csv", row.names=T)

#MR Report
mr_report(dat_bar_par)

