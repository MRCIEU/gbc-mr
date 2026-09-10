library(qqman)
library(data.table)
library(here)
library(dplyr)


gwas <- fread(here("data/G2_N_1875_withHLA_META-with_G1_O_48_withHLA_BothBeagle.meta_with-BETA.csv_With_BETA-SE-OR_n_ConfIntervals_95_with-AllelFreqs-CombinedData_WithrsID.csv"))

str(gwas)

gwas$P <- pnorm(-abs(gwas$BETA / gwas$SE))
hist(gwas$P)

gwas <- subset(gwas, !is.na(P))

png(here("figures/gwas_manhattan.png"))
manhattan(gwas, chr="CHR", bp="BP", snp="SNP", p="P")
dev.off()
