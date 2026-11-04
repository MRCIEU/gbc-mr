# make a new phewas

library(dplyr)
library(ggplot2)
library(here)
library(data.table)

phewas_results <- fread(here("data/phewas/10_phewas_results_of_8_lead_snps_GBC.csv"))

gwas <- fread(here("data/G2_N_1875_withHLA_META-with_G1_O_48_withHLA_BothBeagle.meta_with-BETA.csv_With_BETA-SE-OR_n_ConfIntervals_95_with-AllelFreqs-CombinedData_WithrsID.csv"))

gbc_tophits <- unique(phewas_results$rsid)
gwas_tophits <- gwas[gwas$SNP %in% gbc_tophits, ]

phewas_gwas_merged <- inner_join(phewas_results, gwas_tophits %>% select(SNP, A1_MinorAllele, A2_MajorAllele, BETA), by=c("rsid"="SNP"))

str(phewas_gwas_merged)
phewas_gwas_merged$beta_aligned <- phewas_gwas_merged$beta * sign(phewas_gwas_merged$BETA)
flip_index <- which(phewas_gwas_merged$A1_MinorAllele != phewas_gwas_merged$ea)
phewas_gwas_merged$beta_aligned[flip_index] <- phewas_gwas_merged$beta_aligned[flip_index] * -1


# Curated categories (Jev plus manual overrides), not the raw jev_category.
trait_categories <- fread(here("results/phewas-trait-categories-final.csv"))
trait_categories$category %>% unique()
trait_categories$category <- gsub(" measures", "", trait_categories$category)
trait_categories$category <- gsub("disease of ", "", trait_categories$category)
trait_categories$category <- gsub(" system and connective tissue", "", trait_categories$category)
trait_categories$category <- gsub(" and subcutaneous tissue", "", trait_categories$category)
trait_categories$category <- gsub(" disease", "", trait_categories$category)
trait_categories$category <- gsub(" system", "", trait_categories$category)


phewas_gwas_merged <- inner_join(phewas_gwas_merged, trait_categories %>% select(id, category), by=c("id"="id"))

# p = 0 gives an infinite -log10(p); cap at the smallest non-zero p and mark it.
phewas_gwas_merged$p_capped <- phewas_gwas_merged$p == 0
phewas_gwas_merged$p[phewas_gwas_merged$p_capped] <- min(phewas_gwas_merged$p[!phewas_gwas_merged$p_capped])

# Largest category (most traits) first.
category_size <- phewas_gwas_merged %>% distinct(category, trait) %>% count(category, sort = TRUE)
phewas_gwas_merged$category <- factor(phewas_gwas_merged$category, levels = category_size$category)

# With space="free_x" each panel's width is its number of traits plus this
# padding on each side, so single-trait panels stay wide enough for the strip.
panel_padding <- 2

p <- ggplot(phewas_gwas_merged, aes(x=trait, y=-log10(p) * sign(beta_aligned))) +
  geom_point(aes(colour=rsid)) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c(`FALSE` = "p > 0", `TRUE` = "p = 0 (capped)"), name = NULL) +
  facet_grid(. ~category, scales="free_x", space="free_x") +
  scale_x_discrete(expand = expansion(add = panel_padding)) +
  theme(
    strip.text.x = element_text(angle = 90, hjust = 0),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position="bottom"
  ) +
  labs(x="", y="Signed PheWAS -log10(p)", colour="GBC GWAS hit") +
  scale_colour_brewer(type="qual")
p

ggsave(here("figures/phewas-categories.png"), p, width = 14, height = 5, dpi = 300)
