# Extracts the lab results used as predictors from the multilab extracts.
#
# Converts censored results to numbers and derives INR from TT% where INR
# itself wasn't measured. Saves feature_labs.rds for ce_01.

library(data.table)
library(svglite)

# Paths are relative to the project root
datapath    <- "data"
resultspath <- "results"
objpath     <- file.path(resultspath, "model_objects")
figpath     <- file.path(resultspath, "figures")

# Lab codes (short name and national test number)
# Hb:          B -Hb
# Platelets:   B -Trom   (2791)
# Creatinine:  P -Krea   (4600)
# CRP:         P -CRP    (4594)
# MCV:         E -MCV    (1357)
# INR:         P -INR    (4520)
# TT%:         P -TT     (1731)  (only for deriving INR)
# Albumin:     P -Alb    (4586)
# Leukocytes:  B -Leuk   (2218)
# Ferritin:    P -Ferrit (4826)
target_ids <- c("B -Hb", "B -Trom", "P -Krea", "P -CRP", "E -MCV", "P -INR", "P -TT", "P -Alb", "B -Leuk", "P -Ferrit")

# ---------------------------------------------------------------------------
# Read and filter
# ---------------------------------------------------------------------------
# The multilab table comes in chunks of several GB each
file_paths <- list.files(path = datapath, pattern = "^multilab_laboratoriotulos_\\d+_clean\\.rds$", full.names = TRUE)

# Keeps only the target tests of one chunk, freeing the chunk before the next one is read
read_and_filter <- function(file) {
    temp_dt         <- readRDS(file)
    filtered_subset <- temp_dt[tutkimus_lyhenne %in% target_ids]
    rm(temp_dt); gc()
    filtered_subset
}

final_dt <- rbindlist(lapply(file_paths, read_and_filter))

# ---------------------------------------------------------------------------
# Censored results
# ---------------------------------------------------------------------------
# Results are text: "<X" (below detection) becomes X/2, ">X" (above range) becomes X
# Anything else that isn't a number (e.g. free-text status) becomes NA
final_dt[, threshold_val := as.numeric(gsub("[<>]", "", tulos))]
final_dt[, tulos := fcase(
    grepl("^\\s*<", tulos),                      threshold_val / 2,
    grepl("^\\s*>", tulos),                      threshold_val,
    !is.na(suppressWarnings(as.numeric(tulos))), as.numeric(tulos),
    default = NA_real_
)]
final_dt[, threshold_val := NULL]

# TRUE only for the INR values derived from TT% below (for the reported share in ce_01)
final_dt[, inr_derived := FALSE]

# ---------------------------------------------------------------------------
# INR from TT% (Supplementary Methods 3, Figure S1)
# ---------------------------------------------------------------------------
# TT% and INR side by side, one row per patient and draw time
# (mean if a test was recorded twice at the same time)
tt_inr_wide <- dcast(
    final_dt[tutkimus_lyhenne %in% c("P -TT", "P -INR")],
    henkilotunnus + naytteenotto_hetki ~ tutkimus_lyhenne,
    value.var     = "tulos",
    fun.aggregate = mean,
    na.rm         = TRUE
)

# Rename to avoid backtick pain
setnames(tt_inr_wide, c("P -TT", "P -INR"), c("tt", "inr"))

# Fit on draws with both results, excluding extreme anticoagulation (INR >= 10)
dt_paired  <- tt_inr_wide[!is.na(tt) & !is.na(inr)][inr < 10]
mod_tt_inr <- lm(log(inr) ~ log(tt), data = dt_paired)

# Fit statistics reported in Supplementary Methods 3
inr_fit_summary <- summary(mod_tt_inr)
inr_fit_stats <- list(
    coefficients    = coef(inr_fit_summary),
    n_pairs         = nrow(dt_paired),
    r_squared       = inr_fit_summary$r.squared,
    residual_se     = inr_fit_summary$sigma,
    fit_restriction = "paired same-timestamp TT% and INR, INR < 10"
)
print(inr_fit_stats)
saveRDS(inr_fit_stats, file.path(objpath, "inr_tt_regression.rds"))

# Figure S1: density of the paired values with the fitted line
# (smoothScatter, since plain points overplot into a solid block)
draw_inr_fit <- function() {
    smoothScatter(log(dt_paired$tt), log(dt_paired$inr), xlab = "log(TT%)", ylab = "log(INR)", nrpoints = 0)
    abline(mod_tt_inr, col = "red", lwd = 2)
}

# PNG for viewing, SVG for the journal (2000 px at 300 dpi; the density layer stays a raster inside the SVG)
png(file.path(figpath, "inr_tt_regression.png"), width = 2000, height = 2000, res = 300)
draw_inr_fit()
dev.off()

svglite(file.path(figpath, "inr_tt_regression.svg"), width = 2000 / 300, height = 2000 / 300)
draw_inr_fit()
dev.off()

# Draws with TT% but no INR get a predicted INR (back-transformed from the log scale)
# Rounded to 1 decimal like measured INR
inr_from_tt <- tt_inr_wide[!is.na(tt) & is.na(inr), .(
    henkilotunnus,
    naytteenotto_hetki,
    tutkimus_lyhenne = "P -INR",
    tulos            = round(exp(predict(mod_tt_inr, newdata = .SD)), 1),
    inr_derived      = TRUE
)]

# ---------------------------------------------------------------------------
# Combine and save
# ---------------------------------------------------------------------------
# Only the columns ce_01 uses
keep_cols <- c("henkilotunnus", "naytteenotto_hetki", "tutkimus_lyhenne", "tulos", "inr_derived")

# INR rows with an unreadable result (NA) are dropped first,
# or unique() below would keep the NA row and discard the derived value!
final_dt <- final_dt[!(tutkimus_lyhenne == "P -INR" & is.na(tulos)), ..keep_cols]

# One result per patient, draw time and test
final_dt <- unique(rbind(final_dt, inr_from_tt), by = c("henkilotunnus", "naytteenotto_hetki", "tutkimus_lyhenne"))

saveRDS(final_dt, file.path(datapath, "feature_labs.rds"))
