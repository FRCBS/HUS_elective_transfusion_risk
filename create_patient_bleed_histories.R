# Sums recorded surgical blood loss per patient and hospital episode.
#
# Saves patient_bleed_histories.rds for ce_01 (prior bleeding predictors).

library(data.table)

# Paths are relative to the project root
datapath    <- "data"
resultspath <- "results"

# Bleeding records from the surgical follow-up form (swabs, suction, drains, dressings)
bl <- readRDS(file.path(datapath, "apotti_seurantalomake_leikkausvuodot.rds"))

# An empty episode id ("") isn't an episode (in the diagnosis table it's shared by thousands of patients)
# Such records are grouped per patient and day instead, so they don't merge into one long pseudo-episode
bl[, shjakso_id := as.character(shjakso_id)]
bl[trimws(shjakso_id) == "", shjakso_id := NA]
message("Bleeding records without an episode id: ", bl[is.na(shjakso_id), .N])
bl[is.na(shjakso_id), shjakso_id := paste0("day_", as.IDate(mittaus_hetki))]

# Total loss per episode, dated by its last measurement
epsums <- (
    bl
    [, .(blsum = sum(arvo, na.rm = TRUE), bl_date = as.IDate(max(mittaus_hetki))), by = c("henkilotunnus", "shjakso_id")]
)

# Episodes without any recorded loss carry no information
epsums <- epsums[blsum > 0]

saveRDS(epsums, file.path(resultspath, "patient_bleed_histories.rds"))
