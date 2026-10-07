# Maps visit diagnoses to the Elixhauser comorbidity categories.
#
# One row per visit (patient, episode and contact date), with 0/1 flags for the
# 31 categories and the van Walraven index. Saves elixhauser_visits.rds for ce_01.

library(data.table)
library(comorbidity)

# Paths are relative to the project root
datapath    <- "data"
resultspath <- "results"

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------
# The diagnosis table covers the whole hospital, but only surgical patients are ever used
# (filtering first skips mapping the vast majority of rows that ce_01 would discard anyway)
surgical_cohort <- unique(
    readRDS(file.path(datapath, "apotti_leikkaus_operaatio_henkilosto.rds"))[!is.na(henkilotunnus), henkilotunnus]
)

dx <- readRDS(file.path(datapath, "apotti_diagnoosi_out.rds"))[henkilotunnus %in% surgical_cohort]
rm(surgical_cohort); gc()

# ---------------------------------------------------------------------------
# Clean
# ---------------------------------------------------------------------------
# ICD-10 codes without dots or other punctuation, e.g. "I25.1" -> "I251"
dx <- (
    dx
    [, .(
        shjakso_id    = shjakso_id,
        henkilotunnus = henkilotunnus,
        visit_date    = as.IDate(kontakti_alku_hetki),
        icd10         = gsub("[^[:alnum:]]", "", diagnoosi_koodi)
    )]
)
gc()

# Some rows have an empty episode id (""), shared by tens of thousands of patients, so it isn't an episode
dx[trimws(shjakso_id) == "", shjakso_id := NA]
message("Diagnosis rows without an episode id: ", dx[is.na(shjakso_id), .N])

# A visit is one patient, episode and contact date, so a visit only gets the codes recorded on it
# (patient in the key, so no codes are ever pooled across patients!)
# Integer group number, not a pasted string id (millions of strings would exhaust memory!)
dx[, visit_id := .GRP, by = .(henkilotunnus, shjakso_id, visit_date)]
dx[, shjakso_id := NULL]

# Visit ids, patients and dates kept aside, codes alone go to the mapping
visit_meta <- unique(dx[, .(visit_id, henkilotunnus, visit_date)])
dx_codes   <- dx[, .(visit_id, icd10)]

# dx is dropped before the mapping, which is the memory peak
rm(dx); gc()

# ---------------------------------------------------------------------------
# Elixhauser categories
# ---------------------------------------------------------------------------
# One row per visit, 0/1 per category (Quan et al. ICD-10 coding algorithm)
elix_mapped <- comorbidity(
    x       = dx_codes,
    id      = "visit_id",
    code    = "icd10",
    map     = "elixhauser_icd10_quan",
    assign0 = FALSE
)
rm(dx_codes); gc()

# ---------------------------------------------------------------------------
# Comorbidity index
# ---------------------------------------------------------------------------
# van Walraven (2009) weighted score: integer weights (positive and negative) from a mortality model,
# a single numeric summary next to the category flags
# score() returns a plain vector in elix_mapped's row order, so it's attached directly (not merged)
comorbidity_index <- score(elix_mapped, weights = "vw", assign0 = FALSE)
setDT(elix_mapped)
elix_mapped[, comorbidity_index := comorbidity_index]

# Patient and date back on (ce_01 counts these rows as visits within its look-back window)
elix_visits <- elix_mapped[visit_meta, on = "visit_id"]

saveRDS(elix_visits, file.path(resultspath, "elixhauser_visits.rds"))
