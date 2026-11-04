# Apply manual curation to the Jev trait categories.
#
# Usage: Rscript scripts/finalise-trait-categories.R [jev.csv] [output.csv]
# Keeps Jev's answer alongside the final category and records where each final
# category came from (jev, rule or manual).

library(data.table)

args <- commandArgs(trailingOnly = TRUE)
jev_file <- if (length(args) >= 1) args[1] else "results/phewas-trait-categories.csv"
out_file <- if (length(args) >= 2) args[2] else "results/phewas-trait-categories-final.csv"

x <- fread(jev_file)
setnames(x, c("category", "probability", "confidence"),
         c("jev_category", "jev_probability", "jev_confidence"))
x[, `:=`(category = jev_category, category_source = "jev", note = "")]

set_category <- function(rows, new_category, source, new_note) {
  x[rows, `:=`(category = new_category, category_source = source, note = new_note)]
}

# eQTL, pQTL and metabolite GWAS, including the GWAS Catalog plasma protein and
# metabolite datasets that OpenGWAS imports under ebi-a ids.
set_category(grepl("^eqtl-", x$id), "molecular", "rule", "eQTL dataset")
set_category(grepl("^prot-", x$id) |
               x$id %in% c("ebi-a-GCST90019415", "ebi-a-GCST90019439", "ebi-a-GCST90019461",
                           "ebi-a-GCST90019466", "ebi-a-GCST90019471", "ebi-a-GCST90019484"),
             "molecular", "rule", "Plasma protein (pQTL) dataset")
set_category(grepl("^met-", x$id) | x$id == "ebi-a-GCST90092942",
             "molecular", "rule", "Metabolite dataset")

set_category(grepl("sex hormone-binding globulin", x$trait, ignore.case = TRUE),
             "reproductive measures", "rule", "SHBG")
set_category(grepl("allergic rhinitis", x$trait, ignore.case = TRUE),
             "immune", "rule", "Allergic rhinitis / hayfever")

# Jev confidence < 0.5, assigned by hand (rows already handled above are listed
# for completeness and keep the rule-based category).
manual <- fread(text = "
id|category|note
finn-b-RX_STATIN|metabolic disease|Statin use as a proxy for hyperlipidaemia
ukb-b-18200|other|Absence of a mixed set of doctor-diagnosed conditions
finn-b-K11_REIMB_202|disease of digestive system|FinnGen K11 (digestive) chapter reimbursement endpoint
ukb-b-8121|behavioural measures|Dietary behaviour
ukb-b-6541|other|Healthcare use (consultant speciality)
ukb-a-267|anthropometric measures|Impedance body composition
ukb-b-14540|anthropometric measures|Impedance body composition; matches ukb-a-267
ukb-b-14961|other|Unspecified medical condition
ubm-a-2016|physiological measures|Brain MRI functional connectivity measure
ukb-b-3563|behavioural measures|Supplement use
", sep = "|")
stopifnot(all(manual$id %in% x$id))
for (i in seq_len(nrow(manual))) set_category(x$id == manual$id[i], manual$category[i], "manual", manual$note[i])

setcolorder(x, c("id", "trait", "category", "category_source", "note",
                 "jev_category", "jev_probability", "jev_confidence", "model"))
fwrite(x, out_file)
message(nrow(x), " traits -> ", out_file)
print(x[, .N, by = category][order(-N)])
print(x[, .N, by = category_source])
