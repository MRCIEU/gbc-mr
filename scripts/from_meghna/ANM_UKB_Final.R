## MR Analysis: Age At Menopause (UKB) ##
## South Asian Population, n = 2,055 ##
# Meghna Sequeira #
# 2026 #


# Set working directory
setwd("/Users/meghnasequeira/MR_Data/UKB_Age_at_menopause")

# Run libraries
library(data.table)
library(R.utils)
library(TwoSampleMR)
library(devtools)

######################################################################################################################################################################

#VCF to CSV

# Read VCF file
read_data <- fread("ukb-e-3581_CSA.vcf.gz",skip = "#CHROM",header = TRUE,sep = "\t")
head(read_data)

# Rename first column
setnames(read_data, "#CHROM", "CHROM")

# Split FORMAT names
format_names <- unlist(strsplit(read_data$FORMAT[1], ":"))

# Split sample column values
split_values <- tstrsplit(read_data[[10]], ":")

# Convert to dataframe
split_df <- as.data.table(split_values)

# Assign column names
setnames(split_df, format_names)

# Combine with original VCF columns
final_data <- cbind(read_data, split_df)

# Optional rename
setnames(final_data,
         old = c("ES", "SE", "LP", "AF","ID"),
         new = c("BETA", "SE", "LOGP", "Allele_Frequency","rsid"))

# Save CSV
data_converted <- fwrite(final_data, "ukb-e-3581_CSA.csv")

######################################################################################################################################################################

#LogP to Log conversion

data <- fread("ukb-e-3581_CSA.csv")
head(data)
dim(data)

#-----------------------------------------------------
# Read GWAS file (exposure dataset)
#-----------------------------------------------------
result_g <- fread("ukb-e-3581_CSA.csv")
dim(result_g)
cat(colnames(result_g),sep = "\n")
head(result_g)
#----------------------------------------------------------
#convert in P
#----------------------------------------------------------
result_g$P <- 10^(-result_g$LOGP)
write.csv(result_g, "ukb-e-3581_CSA_Pvalue.csv", row.names = FALSE)
head(result_g)

######################################################################################################################################################################

#Exposure preparation, outcome preparation, and MR

data_1 <- fread("ukb-e-3581_CSA_Pvalue.csv")
dim(data_1)
head(data_1)

#Create exposure data==========================================================================================
exposure_data <- read_exposure_data(
  filename = "ukb-e-3581_CSA_Pvalue.csv",
  clump = TRUE,
  sep = ",",
  phenotype_col = "Phenotype",
  snp_col = "rsid",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "Allele_Frequency",  
  effect_allele_col = "REF",  
  other_allele_col = "ALT",  
  pval_col = "P",  
  min_pval = 1e-200,  
  log_pval = FALSE,  
  chr_col = "CHROM",  
  pos_col = "POS", 
  clump_kb = 10000,  
  clump_r2 = 0.001,  
  clump_p1 = 1,
  bfile = "/Users/meghnasequeira/MR_Data/Kentistou_Age_at_Menarche/1000_genome_SAS/1000G_SAS",  
  plink_bin = "/usr/local/bin/plink")

dim(exposure_data)

write.csv(exposure_data,"1_exposure_data.csv",row.names = FALSE)


#==========================================================================================================
data_2 <- fread("3_G1_O_48_BeagleImp_Study_MichiHLA_Merged_FullyCleaned_PostImpQC_maf-geno-0.005_Logistic_Age_Gen_10PCs_OR.assoc.logistic.csv")
dim(data_2)
head(data_2)

#Create outcome data========================================================================================

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

#Harmonization==================================================================================================

dat1 <- harmonise_data(exposure_data, outcome_data, action = 3)
write.csv(dat1,"harmonize_action_3.csv ", row.names = FALSE)

#Mr analysis=====================================================================================================

mr1 <- mr(dat1)
write.csv(mr1," MR_OUTPUT.csv", row.names = FALSE)

mr_ukb_anm <- read_csv(" MR_OUTPUT.csv")
dat_ukb_anm <- read_csv("harmonize_action_3.csv ")

#Scatter plot of MR
p1 <- mr_scatter_plot(mr_ukb_anm, dat_ukb_anm)
p1[[1]] # Save as UKB_ANM_Scatterplot.png

#Odds ratios
generate_odds_ratios(mr_ukb_anm)
write.csv(generate_odds_ratios(mr_ukb_anm), file = "generate_odds_ratios.csv", row.names=T)

#MR Heterogeneity
mr_heterogeneity(dat_ukb_anm)
write.csv(mr_heterogeneity(dat_ukb_anm), file = "mr_heterogeneity.csv", row.names=T)

#MR Forest plots (singlesnp)
res_single <- mr_singlesnp(dat_ukb_anm)
write.csv(mr_singlesnp(dat_ukb_anm), file = "res_single.csv", row.names=T)
p2 <- mr_forest_plot(res_single)
p2[[1]] # Save as UKB_ANM_Forestplot1.png

#res_single2 <- mr_singlesnp(dat_ukb_anm, all_method=c("mr_ivw", "mr_two_sample_ml"))
#p2 <- mr_forest_plot(res_single2) 
#p2[[1]] # Save as UKB_ANM_Forestplot2.png

#MR Leave one out
res_loo <- mr_leaveoneout(dat_ukb_anm)
write.csv(mr_leaveoneout(dat_ukb_anm), file = "mr_leaveoneout.csv", row.names=T)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]] # Save as UKB_ANM_Loo.png

#MR Pleiotropy test
mr_pleiotropy_test(dat_ukb_anm)
write.csv(mr_pleiotropy_test(dat_ukb_anm), file = "mr_pleiotropy_test.csv", row.names=T)

#MR Report
mr_report(dat_ukb_anm)

