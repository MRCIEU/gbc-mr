## MR Analysis: Age at Menopause (Day) ##
## European Population, n = 69,360 ##
# Meghna Sequeira #
# 2026 #

# Set working directory
setwd("/Users/meghnasequeira/MR_Data/Day_Age_at_menopause")

# Run libraries
library(data.table)
library(R.utils)
library(TwoSampleMR)
library(devtools)

######################################################################################################################################################################

#TXT to CSV

# Read the text file
data <- read.delim("Menopause_HapMap2_DayNG2015_18112015.txt", header = TRUE)

# Save as CSV
write.csv(data, "Menopause_HapMap2_DayNG2015_18112015.csv", row.names = FALSE)

# Check CSV
reprogen <- read_csv("Menopause_HapMap2_DayNG2015_18112015.csv")
dim(reprogen)
head(reprogen)

# Read PLINK reference file
bim <- read.table("/Users/meghnasequeira/MR_Data/Kentistou_Age_at_Menarche/1000_genome_SAS/1000G_SAS.bim")

#Assign column names
colnames(bim) <- c("CHR", "MarkerName", "CM", "POS", "A1", "A2")

# Add chromosome number and position to data
hapmap <- merge(reprogen, bim[,c("MarkerName","CHR","POS")],
               by="MarkerName", all.x=TRUE)

write.csv(hapmap, "hapmap.csv", row.names = FALSE)

# Check CSV
hapmapcsv <- read_csv("hapmap.csv")
dim(hapmapcsv)
head(hapmapcsv)


#Create exposure data=====================================================================================

exposure_data <- read_exposure_data(
  filename = "hapmap.csv",
  clump = TRUE,
  sep = ",",
  phenotype_col = "PHENOTYPE", #Optional column name for the column with phenotype name corresponding to the SNP
  snp_col = "MarkerName",
  beta_col = "effect",
  se_col = "stderr",
  eaf_col = "HapMap_eaf",  
  effect_allele_col = "allele1",  
  other_allele_col = "allele2",  
  pval_col = "p",  
  min_pval = 1e-200,  
  log_pval = FALSE,  
  chr_col = "CHR",  
  pos_col = "POS", 
  clump_kb = 10000,  
  clump_r2 = 0.001,  
  clump_p1 = 1,
  bfile = "/Users/meghnasequeira/MR_Data/Kentistou_Age_at_Menarche/1000_genome_SAS/1000G_SAS",  
  plink_bin = "/usr/local/bin/plink"
)

write.csv(exposure_data,"hapmap_exposure_data.csv",row.names = FALSE)

#Create outcome data=====================================================================================
outcome_data <- read_csv("1_outcome_data.csv")

#Harmonizaton===================================================================================================
dat1_p <- harmonise_data(exposure_data, outcome_data, action = 3)
write.csv(dat1_p,"harmonize_action_3.csv ", row.names = FALSE)
dim(dat1_p)

#Mr analysis=====================================================================================================

har <- read_csv("harmonize_action_3.csv ")

mr_results <- mr( # Can't perform MR simple or MR weighted mode
  har,
  method_list = c(
    "mr_ivw",
    "mr_egger_regression",
    "mr_weighted_median"
  )
) 
write.csv(mr_results," mr_results.csv", row.names = FALSE)
mr_results <- read_csv(" mr_results.csv")

#Odds Ratio
generate_odds_ratios(mr_results)
write.csv(generate_odds_ratios(mr_results), file = "generate_odds_ratios.csv", row.names=T)

#Scatter plot of MR
p1 <- mr_scatter_plot(mr_results, har)
p1[[1]] 

#MR Heterogeneity
mr_heterogeneity(har)
write.csv(mr_heterogeneity(har), file = "mr_heterogeneity.csv", row.names=T)

#MR Forest plots (singlesnp)
res_single <- mr_singlesnp(har)
write.csv(mr_singlesnp(har), file = "res_single.csv", row.names=T)
p2 <- mr_forest_plot(res_single)
p2[[1]] 

#res_single2 <- mr_singlesnp(har, all_method=c("mr_ivw", "mr_two_sample_ml"))
#p2 <- mr_forest_plot(res_single2) 
#p2[[1]] # Save as Barban_parity_Forestplot2.png

#MR Leave one out
res_loo <- mr_leaveoneout(har)
write.csv(mr_leaveoneout(har), file = "mr_leaveoneout.csv", row.names=T)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]] 

#MR Pleiotropy test
mr_pleiotropy_test(har)
write.csv(mr_pleiotropy_test(har), file = "mr_pleiotropy_test.csv", row.names=T)

#MR Report
mr_report(har)

