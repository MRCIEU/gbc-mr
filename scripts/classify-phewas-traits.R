# Assign each PheWAS trait to one broad category using TypeSafe's Jev model
# (System One Choice question; https://docs.typesafe.ai/api/system-one.md).
#
# Usage: Rscript scripts/classify-phewas-traits.R [phewas.csv] [output.csv]
# The API key is read from ~/.typesafe (or TYPESAFE_API_KEY) and never printed.
# Rows already in the output are skipped, so an interrupted run can be resumed.

library(data.table)
library(httr2)

args <- commandArgs(trailingOnly = TRUE)
phewas_file <- if (length(args) >= 1) args[1] else "data/phewas/10_phewas_results_of_8_lead_snps_GBC.csv"
out_file <- if (length(args) >= 2) args[2] else "results/phewas-trait-categories.csv"
model <- Sys.getenv("TYPESAFE_MODEL", "jev-latest")

categories <- c(
  "infectious disease",
  "neoplasm",
  "disease of blood and blood-forming organs",
  "metabolic disease",
  "neurological disease",
  "mental disorder",
  "nervous system disorders",
  "disease of eye and adnexa",
  "disease of ear and mastoid process",
  "disease of circulatory system",
  "disease of digestive system",
  "disease of skin and subcutaneous tissue",
  "disease of musculoskeletal system and connective tissue",
  "disease of genitourinary system",
  "physiological measures",
  "reproductive measures",
  "psychiatric measures",
  "anthropometric measures",
  "environmental measures",
  "socioeconomic measures",
  "behavioural measures"
)

api_key <- Sys.getenv("TYPESAFE_API_KEY")
if (!nzchar(api_key)) api_key <- trimws(readLines(path.expand("~/.typesafe"), warn = FALSE)[1])
stopifnot(nzchar(api_key))

classify_trait <- function(id, trait) {
  body <- list(
    model = model,
    # The OpenGWAS id gives context for opaque traits such as eQTL gene ids.
    state = list(trait = trait, opengwas_id = id),
    questions = list(category = list(
      type = "choice",
      instructions = "Which category best describes this GWAS trait?",
      criteria = as.list(setNames(categories, categories))
    ))
  )
  resp <- request("https://api.typesafe.ai/v1/systemone") |>
    req_auth_bearer_token(api_key) |>
    req_body_json(body, auto_unbox = TRUE) |>
    req_retry(max_tries = 5, is_transient = \(r) resp_status(r) %in% c(429, 500, 502, 503, 529)) |>
    req_perform() |>
    resp_body_json()
  answer <- resp$answers$category
  stopifnot(answer$choice %in% categories)
  data.table(id = id, trait = trait, category = answer$choice,
             probability = answer$probabilities[[answer$choice]],
             confidence = if (is.null(answer$confidence)) NA_real_ else answer$confidence,
             model = resp$model)
}

traits <- unique(fread(phewas_file)[, .(id, trait)])
done <- if (file.exists(out_file)) fread(out_file) else data.table()
todo <- if (nrow(done)) traits[!done, on = "id"] else traits
message(nrow(traits), " traits; ", nrow(done), " already classified; ", nrow(todo), " to submit")

dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
for (i in seq_len(nrow(todo))) {
  res <- tryCatch(classify_trait(todo$id[i], todo$trait[i]), error = function(e) {
    message("Failed: ", todo$id[i], " (", conditionMessage(e), ")")
    NULL
  })
  if (is.null(res)) next
  fwrite(res, out_file, append = file.exists(out_file))
  if (i %% 25 == 0) message(i, "/", nrow(todo))
}

result <- fread(out_file)
message(nrow(result), " of ", nrow(traits), " traits classified -> ", out_file)
print(result[, .N, by = category][order(-N)])
