# Builds the modelling dataset (dtbin) from the raw EHR extracts.
#
# One row per elective procedure, with the transfusion outcome, patient-level
# predictors, the train/test split and the tuning folds.

library(data.table)

# Paths are relative to the project root
datapath    <- "data"
resultspath <- "results"
objpath     <- file.path(resultspath, "model_objects")

PRED_HORIZON_DAYS   <- 7L     # predictors must be recorded at least this many days before OR entry
OUTCOME_WINDOW_DAYS <- 7L     # transfusions count up to this many days after OR exit
MIN_PROC_N          <- 50L    # minimum procedures per primary NCSP code
ELIX_LOOKBACK_DAYS  <- 365L   # comorbidity look-back before the horizon
MED_LOOKBACK_DAYS   <- c(10L, 30L, 90L)  # medication look-backs, by drug class
TEST_FRAC           <- 0.20   # share of procedures targeted for the test set
N_CV_FOLDS          <- 5L     # patient-grouped tuning folds in the training set

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
# Most recent lab result drawn at least PRED_HORIZON_DAYS before each OR entry,
# plus how many days old it was (no upper age limit, see the lab notes below)
latest_lab_before_or <- function(lab_dt, or_entries, time_col, result_col, out_col = result_col) {
    lab_slim <- lab_dt[, .SD, .SDcols = c("henkilotunnus", time_col, result_col)]
    setnames(lab_slim, c(time_col, result_col), c("lab_time", "lab_value"))

    out <- (
        # Every lab of the patient against every one of their OR entries
        merge(lab_slim, or_entries, by = "henkilotunnus", allow.cartesian = TRUE)
        [, diff := difftime(potilas_saliin_hetki, lab_time, units = "hours")]
        # Only results from before the horizon
        [diff >= PRED_HORIZON_DAYS * 24L]
        # The closest one to the horizon
        [, .SD[which.min(diff)], by = c("henkilotunnus", "potilas_saliin_hetki")]
        [, .(
            henkilotunnus        = henkilotunnus,
            potilas_saliin_hetki = potilas_saliin_hetki,
            lab_value            = lab_value,
            days_since_drawn     = as.numeric(diff) / 24
        )]
    )

    # e.g. latest_hb and latest_hb_days_since_drawn
    setnames(out, c("lab_value", "days_since_drawn"), c(out_col, paste0(out_col, "_days_since_drawn")))
    out[]
}

# Fills missing values of one continuous variable with lm() predictions from complete cases
# Modifies dt in place (data.table semantics!)
impute_by_regression <- function(dt, var, predictors) {
    complete_idx <- !is.na(dt[[var]])
    missing_idx  <-  is.na(dt[[var]])

    # Nothing to impute (just in case)
    if (!any(missing_idx)) return(invisible(NULL))

    fml       <- as.formula(paste(var, "~", paste(predictors, collapse = " + ")))
    fit       <- lm(fml, data = dt[complete_idx])
    predicted <- predict(fit, newdata = dt[missing_idx])

    # Clip to the observed 1st to 99th percentile, so no extreme extrapolations
    bounds    <- quantile(dt[[var]], probs = c(0.01, 0.99), na.rm = TRUE)
    predicted <- pmin(pmax(predicted, bounds[1L]), bounds[2L])

    set(dt, which(missing_idx), var, predicted)
    message(sprintf("  imputed %s: %d missing -> regression (R2 = %.3f)", var, sum(missing_idx), summary(fit)$r.squared))
    invisible(NULL)
}

# ---------------------------------------------------------------------------
# Operations
# One row per operation and primary surgeon
# ---------------------------------------------------------------------------
# Source table has one row per operation AND assigned staff member
leikkaus_raw <- readRDS(file.path(datapath, "apotti_leikkaus_operaatio_henkilosto.rds"))

# Professional titles from most to least senior (same tiers as ce_03's title_map)
title_seniority <- c(
    "ylilääkäri",
    "osastonylilääkäri",
    "apulaisylilääkäri",
    "apulaisylihammaslääkäri",
    "erikoislääkäri",
    "erikoishammaslääkäri",
    "osastonlääkäri",
    "osastonhammaslääkäri",
    "erikoistuva lääkäri",
    "erikoistuva hammaslääkäri",
    "lääkäri",
    "päivystävä lääkäri",
    "amanuenssi",
    "lääketieteen kandidaatti",
    "lääketieteen opiskelija",
    "kätilö",
    "tutkija",
    "tutkimusassistentti"
)

# Rank 1 is most senior; unlisted titles come after every listed one, missing titles last
leikkaus_raw[, title_rank := match(toimenpide_tekija_ammatti_selite, title_seniority)]
leikkaus_raw[
    is.na(title_rank),
    title_rank := length(title_seniority) + 1L + is.na(toimenpide_tekija_ammatti_selite)
]

# Random key for breaking ties between equally senior staff
# (own seed, so it doesn't affect the train/test split later)
set.seed(20260814)
leikkaus_raw[, rand_key := runif(.N)]

# Most senior staff member first within each operation, keep that row as the primary surgeon
# Rows without an operation id can't be grouped, so they're all kept
n_before_collapse <- nrow(leikkaus_raw)
setorder(leikkaus_raw, operaatio_id, title_rank, rand_key)
keep_idx <- leikkaus_raw[, .I[if (is.na(operaatio_id[1L])) seq_len(.N) else 1L], by = operaatio_id]$V1
leikkaus_raw <- leikkaus_raw[sort(keep_idx)]
leikkaus_raw[, `:=`(title_rank = NULL, rand_key = NULL)]
message("Collapsed ", n_before_collapse - nrow(leikkaus_raw), " multi-personnel duplicate rows in leikkaus_raw")

# Every OR entry of ANY urgency class, for cutting outcome windows short at a reoperation below
# (an emergency reoperation never gets its own row, but its transfusions mustn't count for the elective one!)
all_or_entries <- leikkaus_raw[
    !is.na(henkilotunnus) & !is.na(potilas_saliin_hetki),
    .(henkilotunnus, potilas_saliin_hetki)
]

# Elective operations (urgency classes LPI to LPIII) with complete times and a procedure code
ops <- (
    leikkaus_raw
    [leikkaus_kiireellisyysluokka_selite %in% c("LPI", "LPII", "LPIII")]
    [!is.na(potilas_saliin_hetki) & !is.na(potilas_salista_hetki) & !is.na(paatoimenpide_koodi)]
    [, .(
        henkilotunnus,
        paatoimenpide_koodi,
        paatoimenpide_selite,
        potilas_saliin_hetki,
        potilas_salista_hetki,
        asa_luokka_koodi,
        toimenpide_tekija_id,
        toimenpide_tekija_ammatti_selite
    )]
)

rm(leikkaus_raw); gc()

# ---------------------------------------------------------------------------
# Transfusions
# ---------------------------------------------------------------------------
# Cleaned transfusion records from create_transfusions_clean.R, RBC units with a start time
transf <- (
    readRDS(file.path(resultspath, "transfusions_clean.rds"))
    [henkilotunnus %in% unique(ops$henkilotunnus)]
    [product_type == "RBC"]
    [!is.na(veri_alku_aika)]
    [, .(henkilotunnus, veri_alku_aika, veri_admin_tuote_numero)]
)

# ---------------------------------------------------------------------------
# Outcome window
# From OR entry to OUTCOME_WINDOW_DAYS after OR exit, cut short at discharge or the next operation
# ---------------------------------------------------------------------------
# Next OR entry of the same patient after this operation's exit (any urgency class)
# ops goes in i, so mult = "first" picks the earliest later entry for every ops row
# (the reverse join leaves most rows without a match!)
all_or_entries <- unique(all_or_entries)
setorder(all_or_entries, henkilotunnus, potilas_saliin_hetki)
ops[, next_or_entry := all_or_entries[
    ops,
    on   = .(henkilotunnus, potilas_saliin_hetki > potilas_salista_hetki),
    x.potilas_saliin_hetki,
    mult = "first"
]]
rm(all_or_entries); gc()

# Discharge time is looked up from three sources covering largely different patients:
# 1. ward stay end (most direct, but only for ward-admitted patients, ~19% of surgeries)
# 2. end of the visit containing the OR entry (most day-case surgery)
# 3. end of the diagnosis episode (near-universal last resort; agrees with 2 within 1 h in 99.9%)
disch_cohort <- unique(ops$henkilotunnus)

kaynti  <- readRDS(file.path(datapath, "apotti_kaynti_out.rds"))[henkilotunnus %in% disch_cohort]
osasto  <- readRDS(file.path(datapath, "apotti_osastohoito_out.rds"))[henkilotunnus %in% disch_cohort]
diag_dt <- readRDS(file.path(datapath, "apotti_diagnoosi_out.rds"))[henkilotunnus %in% disch_cohort]

# One row per ward stay episode
osasto_shjakso <- unique(osasto[, .(henkilotunnus, shjakso_id, shjakso_alku_hetki, shjakso_loppu_hetki)])

# Diagnosis rows overlap heavily (one row per code), so they're collapsed to one span per episode first
# Rows with an empty episode id ("") aren't an episode, and would span years per patient
diag_shjakso <- diag_dt[!is.na(shjakso_id) & trimws(shjakso_id) != "", .(
    diag_alku  = min(kontakti_alku_hetki, na.rm = TRUE),
    diag_loppu = max(kontakti_loppu_hetki, na.rm = TRUE)
), by = .(henkilotunnus, shjakso_id)]

rm(osasto, diag_dt); gc()

# End time of the interval in tbl that contains each OR entry
# (missing end == still open; shortest interval wins if several contain it)
match_discharge <- function(ops, tbl, start_col, end_col, out_col) {
    idx <- unique(ops[, .(henkilotunnus, potilas_saliin_hetki)])
    m   <- merge(idx, tbl, by = "henkilotunnus", allow.cartesian = TRUE)
    m   <- m[get(start_col) <= potilas_saliin_hetki & (is.na(get(end_col)) | get(end_col) >= potilas_saliin_hetki)]
    m[, .interval_width := get(end_col) - get(start_col)]
    out <- m[order(.interval_width), .SD[1], by = .(henkilotunnus, potilas_saliin_hetki)]
    setnames(out, end_col, out_col)
    out[, c("henkilotunnus", "potilas_saliin_hetki", out_col), with = FALSE]
}

disch_osasto <- match_discharge(ops, osasto_shjakso, "shjakso_alku_hetki", "shjakso_loppu_hetki", "disch_osasto")
disch_kaynti <- match_discharge(ops, kaynti[, .(henkilotunnus, alku_hetki, loppu_hetki)], "alku_hetki", "loppu_hetki", "disch_kaynti")
disch_diag   <- match_discharge(ops, diag_shjakso, "diag_alku", "diag_loppu", "disch_diag")

rm(kaynti, osasto_shjakso, diag_shjakso); gc()

ops <- Reduce(
    function(x, y) merge(x, y, by = c("henkilotunnus", "potilas_saliin_hetki"), all.x = TRUE),
    list(ops, disch_osasto, disch_kaynti, disch_diag)
)
rm(disch_osasto, disch_kaynti, disch_diag); gc()

# A discharge before OR exit is a data error -> blanked so the next source is used instead
ops[!is.na(disch_osasto) & disch_osasto < potilas_salista_hetki, disch_osasto := as.POSIXct(NA)]
ops[!is.na(disch_kaynti) & disch_kaynti < potilas_salista_hetki, disch_kaynti := as.POSIXct(NA)]
ops[!is.na(disch_diag)   & disch_diag   < potilas_salista_hetki, disch_diag   := as.POSIXct(NA)]

# First available source in priority order
ops[, discharge_hetki := fcoalesce(disch_osasto, disch_kaynti, disch_diag)]
ops[, `:=`(disch_osasto = NULL, disch_kaynti = NULL, disch_diag = NULL)]

# Window ends at the earliest of: fixed window, next OR entry, discharge
# Never before OR exit, the fixed window caps an implausibly late discharge
ops[, `:=`(
    start = potilas_saliin_hetki,
    end   = pmax(
        potilas_salista_hetki,
        pmin(potilas_salista_hetki + OUTCOME_WINDOW_DAYS * 24 * 3600, next_or_entry, discharge_hetki, na.rm = TRUE)
    )
)]
ops[, `:=`(next_or_entry = NULL, discharge_hetki = NULL)]

# ---------------------------------------------------------------------------
# Transfusions within the outcome window
# ---------------------------------------------------------------------------
# Duplicate ops rows would match the same transfusions twice (and inflate the outcome!)
n_dup_ops <- nrow(ops) - uniqueN(ops)
if (n_dup_ops > 0L) {
    message("Dropped ", n_dup_ops, " duplicate ops rows before overlap join")
    ops <- unique(ops)
}
setkey(ops, henkilotunnus, start, end)

# Set transfusion to a zero-length interval
transf[, `:=`(start = veri_alku_aika, end = veri_alku_aika)]
setkey(transf, henkilotunnus, start, end)

# foverlaps pairs every operation with the transfusions inside its window
optr <- (
    foverlaps(ops, transf, type = "any")
    # Count matched units, not .N (an operation without transfusions still gets one all-NA row!!!)
    [, rbc_count := sum(!is.na(veri_alku_aika)), by = .(henkilotunnus, potilas_saliin_hetki)]
    # Back to one row per operation
    [, .SD[1], by = .(henkilotunnus, potilas_saliin_hetki)]
    [, `:=`(
        start                   = NULL,
        end                     = NULL,
        i.start                 = NULL,
        i.end                   = NULL,
        veri_alku_aika          = NULL,
        veri_admin_tuote_numero = NULL
    )]
)

rm(ops, transf); gc()

# Only procedure codes with at least MIN_PROC_N operations, for stable procedure-level estimates
dt <- optr[, proc_n := .N, by = paatoimenpide_koodi][proc_n >= MIN_PROC_N]
rm(optr)

# Keys reused by every predictor block below
cohort     <- unique(dt$henkilotunnus)
surg_dates <- dt[, .(henkilotunnus, potilas_saliin_hetki, surg_date = as.IDate(potilas_saliin_hetki))]
or_entries <- dt[, .(henkilotunnus, potilas_saliin_hetki)]

# ---------------------------------------------------------------------------
# Age and sex
# ---------------------------------------------------------------------------
# Patients with missing or non-binary sex are excluded (they drop out when dtbin is built)
asiakas <- (
    readRDS(file.path(datapath, "asiakas_out.rds"))
    [henkilotunnus %in% cohort]
    [sukupuoli_selite %in% c("Mies", "Nainen")]
    [, sukupuoli_selite := factor(sukupuoli_selite, levels = c("Mies", "Nainen"))]
    [, .(henkilotunnus, sukupuoli_selite, syntymaaika)]
)

# Age in years at the day of surgery
age_sex <- (
    merge(surg_dates, asiakas, by = "henkilotunnus", all.x = TRUE)
    [, age := as.numeric(difftime(surg_date, syntymaaika, units = "weeks") / 52.143)]
    [, .(henkilotunnus, sukupuoli_selite, potilas_saliin_hetki, age)]
)

rm(asiakas)

# ---------------------------------------------------------------------------
# BMI and blood volume
# ---------------------------------------------------------------------------
bmi <- (
    readRDS(file.path(datapath, "mittaukset_bmi.rds"))
    [henkilotunnus %in% cohort]
)

# Most recent height/weight/BMI measured at least PRED_HORIZON_DAYS before surgery
bmi_preop <- (
    merge(surg_dates, bmi, by = "henkilotunnus", all.x = TRUE, allow.cartesian = TRUE)
    [, diff := difftime(surg_date, mittaus_pvm, units = "days")]
    [diff >= PRED_HORIZON_DAYS]
    [, .SD[which.min(diff)], by = .(henkilotunnus, potilas_saliin_hetki)]
)

rm(bmi)

# Blood volume (L): Nadler formula (sex-specific) for BMI <= 30, Lemmens-Bernstein-Brodsky above
# Missing BMI -> missing blood volume, both imputed later
bv <- (
    merge(
        bmi_preop,
        unique(age_sex[, .(henkilotunnus, sukupuoli_selite)], by = "henkilotunnus"),
        by    = "henkilotunnus",
        all.x = TRUE
    )
    [, bv_formula := fcase(
        bmi > 30,                                "LBB",
        bmi <= 30 & sukupuoli_selite == "Mies",   "NadlerM",
        bmi <= 30 & sukupuoli_selite == "Nainen", "NadlerF"
    )]
    # paino = weight (kg), pituus = height (cm)
    [, bv_litres := fcase(
        bv_formula == "LBB",     paino * 70 / sqrt(bmi / 22) / 1000,
        bv_formula == "NadlerM", 0.3669 * (pituus / 100)^3 + 0.03219 * paino + 0.6041,
        bv_formula == "NadlerF", 0.3561 * (pituus / 100)^3 + 0.03308 * paino + 0.1833
    )]
    [, .(henkilotunnus, bmi, bv_litres, potilas_saliin_hetki)]
)

rm(bmi_preop); gc()

# ---------------------------------------------------------------------------
# Preoperative labs
# ---------------------------------------------------------------------------
# All lab results from create_feature_labs_table.R, one analyte per tutkimus_lyhenne code
# inr_derived flags INR values converted from TT% (see create_feature_labs_table.R)
labs <- (
    readRDS(file.path(datapath, "feature_labs.rds"))
    [henkilotunnus %in% cohort]
    [!is.na(naytteenotto_hetki) & !is.na(tulos)]
    [, .(henkilotunnus, naytteenotto_hetki, tulos, tutkimus_lyhenne, inr_derived)]
)

hb         <- labs[tutkimus_lyhenne == "B -Hb",     .(henkilotunnus, naytteenotto_hetki, hb_tulos    = tulos)]
trom       <- labs[tutkimus_lyhenne == "B -Trom",   .(henkilotunnus, naytteenotto_hetki, trom_tulos  = tulos)]
creatinine <- labs[tutkimus_lyhenne == "P -Krea",   .(henkilotunnus, naytteenotto_hetki, creat_tulos = tulos)]
crp        <- labs[tutkimus_lyhenne == "P -CRP",    .(henkilotunnus, naytteenotto_hetki, crp_tulos   = tulos)]
mcv        <- labs[tutkimus_lyhenne == "E -MCV",    .(henkilotunnus, naytteenotto_hetki, mcv_tulos   = tulos)]
albumin    <- labs[tutkimus_lyhenne == "P -Alb",    .(henkilotunnus, naytteenotto_hetki, alb_tulos   = tulos)]
wbc        <- labs[tutkimus_lyhenne == "B -Leuk",   .(henkilotunnus, naytteenotto_hetki, wbc_tulos   = tulos)]
ferritin   <- labs[tutkimus_lyhenne == "P -Ferrit", .(henkilotunnus, naytteenotto_hetki, ferr_tulos  = tulos)]
inr        <- labs[tutkimus_lyhenne == "P -INR",    .(henkilotunnus, naytteenotto_hetki, inr_tulos   = tulos, inr_derived = as.numeric(inr_derived))]

rm(labs); gc()

# Latest qualifying value of each lab, plus its age in days
hb_entries    <- latest_lab_before_or(hb,         or_entries, "naytteenotto_hetki", "hb_tulos",    "latest_hb")
trom_entries  <- latest_lab_before_or(trom,       or_entries, "naytteenotto_hetki", "trom_tulos",  "latest_tromb")
creat_entries <- latest_lab_before_or(creatinine, or_entries, "naytteenotto_hetki", "creat_tulos", "latest_creatinine")
crp_entries   <- latest_lab_before_or(crp,        or_entries, "naytteenotto_hetki", "crp_tulos",   "latest_crp")
mcv_entries   <- latest_lab_before_or(mcv,        or_entries, "naytteenotto_hetki", "mcv_tulos",   "latest_mcv")
alb_entries   <- latest_lab_before_or(albumin,    or_entries, "naytteenotto_hetki", "alb_tulos",   "latest_albumin")
wbc_entries   <- latest_lab_before_or(wbc,        or_entries, "naytteenotto_hetki", "wbc_tulos",   "latest_wbc")
ferr_entries  <- latest_lab_before_or(ferritin,   or_entries, "naytteenotto_hetki", "ferr_tulos",  "latest_ferritin")
inr_entries   <- latest_lab_before_or(inr,        or_entries, "naytteenotto_hetki", "inr_tulos",   "latest_inr")

# Whether each selected INR was derived from TT% (same rows and selection as inr_entries)
# For the reported share in cohort_meta only, not a predictor
inr_derived_entries <- latest_lab_before_or(
    inr, or_entries, "naytteenotto_hetki", "inr_derived", "latest_inr_derived"
)[, .(henkilotunnus, potilas_saliin_hetki, latest_inr_derived)]

rm(hb, trom, creatinine, crp, mcv, albumin, wbc, ferritin, inr, or_entries); gc()

# ---------------------------------------------------------------------------
# ASA class
# ---------------------------------------------------------------------------
# Recorded on the operation itself (known when the surgery is scheduled)
asa <- dt[, .(
    henkilotunnus        = henkilotunnus,
    potilas_saliin_hetki = potilas_saliin_hetki,
    asa                  = asa_luokka_koodi
)]

# ---------------------------------------------------------------------------
# Bleeding history
# ---------------------------------------------------------------------------
# Recorded blood loss per earlier surgical episode, from create_patient_bleed_histories.R
bl <- readRDS(file.path(resultspath, "patient_bleed_histories.rds"))[henkilotunnus %in% cohort]

# Total loss and number of bleeding episodes before the horizon, per operation
bl_histories <- (
    merge(bl, surg_dates, by = "henkilotunnus", all.x = TRUE, allow.cartesian = TRUE)
    # Bleeds after the horizon count as 0
    [, gated_value    := fifelse(bl_date < surg_date - PRED_HORIZON_DAYS, blsum, 0)]
    [, bl_history_sum := sum(gated_value, na.rm = TRUE),     by = .(henkilotunnus, potilas_saliin_hetki)]
    [, n_prev_bl_eps  := sum(gated_value > 0, na.rm = TRUE), by = .(henkilotunnus, potilas_saliin_hetki)]
    [, .(
        total_prev_bl = first(bl_history_sum),
        n_prev_bl_eps = as.numeric(first(n_prev_bl_eps))
    ), by = .(henkilotunnus, potilas_saliin_hetki)]
)

rm(bl); gc()

# ---------------------------------------------------------------------------
# Elixhauser comorbidities
# ---------------------------------------------------------------------------
# One row per hospital visit with 0/1 flags for the 31 categories, from create_elixhauser_table.R
elxh <- readRDS(file.path(resultspath, "elixhauser_visits.rds"))[henkilotunnus %in% cohort]

# The van Walraven index is a numeric score, so it's set aside before the melt
# Missing scores are dropped here (lest we don't get -Inf from max())
score_visits <- elxh[!is.na(comorbidity_index), .(visit_id, henkilotunnus, visit_date, comorbidity_index)]

# Long format, one row per visit and category present
elxh <- melt(
    elxh[, comorbidity_index := NULL],
    id.vars       = c("visit_id", "henkilotunnus", "visit_date"),
    variable.name = "elix_category",
    value.name    = "is_present"
)[is_present == 1]

# Look-back window [horizon - ELIX_LOOKBACK_DAYS, horizon) per operation
surg_dates[, horizon_date := surg_date - PRED_HORIZON_DAYS]
surg_dates[, lookback_1yr := horizon_date - ELIX_LOOKBACK_DAYS]

# Visits within each operation's window (non-equi join)
elxh <- elxh[
    surg_dates,
    on      = .(henkilotunnus, visit_date >= lookback_1yr, visit_date < horizon_date),
    .(henkilotunnus, potilas_saliin_hetki, elix_category),
    nomatch = NULL
]

# Number of visits per category, wide with one column per category
elxh <- elxh[, .(visit_count_1yr = as.numeric(.N)), by = .(henkilotunnus, potilas_saliin_hetki, elix_category)]
elxh <- dcast(elxh, henkilotunnus + potilas_saliin_hetki ~ elix_category, value.var = "visit_count_1yr", fill = 0)

# Back to every operation, 0 where no qualifying visit
elxh <- elxh[surg_dates[, .(henkilotunnus, potilas_saliin_hetki)], on = .(henkilotunnus, potilas_saliin_hetki)]
setnafill(elxh, fill = 0, cols = setdiff(names(elxh), c("henkilotunnus", "potilas_saliin_hetki")))

# Highest van Walraven index within the same window, 0 if none
comorb_idx <- (
    score_visits[
        surg_dates,
        on      = .(henkilotunnus, visit_date >= lookback_1yr, visit_date < horizon_date),
        .(henkilotunnus, potilas_saliin_hetki, comorbidity_index),
        nomatch = NULL
    ]
    [, .(comorbidity_index = max(comorbidity_index, na.rm = TRUE)), by = .(henkilotunnus, potilas_saliin_hetki)]
)

elxh <- merge(elxh, comorb_idx, by = c("henkilotunnus", "potilas_saliin_hetki"), all.x = TRUE)
elxh[is.na(comorbidity_index), comorbidity_index := 0]

surg_dates[, `:=`(lookback_1yr = NULL, horizon_date = NULL)]

# ---------------------------------------------------------------------------
# Medications
# Drug classes, ATC codes and windows as in Supplementary Table S1
# ---------------------------------------------------------------------------
# Horizon and the three look-back starts per operation (in POSIXct seconds)
med_or <- dt[, .(henkilotunnus, op_time = potilas_saliin_hetki)]
med_or[, horizon_time := op_time - PRED_HORIZON_DAYS * 24 * 3600]
med_or[, lb_10d := horizon_time - MED_LOOKBACK_DAYS[1] * 24 * 3600]
med_or[, lb_30d := horizon_time - MED_LOOKBACK_DAYS[2] * 24 * 3600]
med_or[, lb_90d := horizon_time - MED_LOOKBACK_DAYS[3] * 24 * 3600]

# Prescriptions from the current EHR (Apotti, in chunked files) and the legacy system (Uranus)
# Medication history entries are patient-reported, not prescriptions, so they're excluded
# Apotti (many files)
ap_med <- rbindlist(lapply(
    sort(list.files(datapath, pattern = "^apotti_laakitys_\\d{3}_clean\\.rds$", full.names = TRUE)),
    function(f) {
        chunk <- readRDS(f)
        chunk[
            henkilotunnus %in% cohort & laakemaarays_luokka_selite != "Lääkehistoria",
            .(henkilotunnus, atc_koodi = as.character(atc_koodi), alku_hetki, loppu_hetki, src = "ap")
        ]
    }
))

# Uranus (single file)
ur_med <- readRDS(file.path(datapath, "uranus_laakitys.rds"))
ur_med <- ur_med[
    henkilotunnus %in% cohort,
    .(henkilotunnus, atc_koodi = as.character(atc_koodi), alku_hetki, loppu_hetki, src = "ur")
]

meds <- rbind(ap_med, ur_med)
rm(ap_med, ur_med); gc()

# Prescriptions started before the horizon and still active within the longest (90 d) window
# (missing end date == still active)
meds_win <- (
    merge(med_or, meds, by = "henkilotunnus", allow.cartesian = TRUE)
    [alku_hetki <= horizon_time & (is.na(loppu_hetki) | loppu_hetki >= lb_90d)]
)

rm(meds); gc()

# Still active within the shorter windows?
meds_win[, active_30d := is.na(loppu_hetki) | loppu_hetki >= lb_30d]
meds_win[, active_10d := is.na(loppu_hetki) | loppu_hetki >= lb_10d]

# 10-day classes
# Aspirin and P2Y12 inhibitors are flagged separately to build DAPT below
meds_win[, anticoagulant     := active_10d & grepl("^(B01AA|B01AB|B01AE|B01AF|B01AX)", atc_koodi)]
meds_win[, antiplatelet      := active_10d & grepl("^B01AC",                            atc_koodi)]
meds_win[, aspirin           := active_10d & grepl("^B01AC06",                          atc_koodi)]
meds_win[, p2y12             := active_10d & grepl("^(B01AC04|B01AC22|B01AC24)",        atc_koodi)]
meds_win[, nsaid             := active_10d & grepl("^(M01A|M01B)",                      atc_koodi)]

# 30-day classes
meds_win[, corticosteroid    := active_30d & grepl("^(H02A|D07)",                       atc_koodi)]
meds_win[, coag_factor       := active_30d & grepl("^B02BD",                            atc_koodi)]
meds_win[, hemostatic        := active_30d & grepl("^B02BX",                            atc_koodi)]

# 90-day classes (every row of meds_win is already within 90 d)
meds_win[, iron_oral         := grepl("^(B03AA|B03AB)", atc_koodi)]
meds_win[, iron_iv           := grepl("^B03AC",         atc_koodi)]
meds_win[, esa               := grepl("^B03XA",         atc_koodi)]
meds_win[, cytostatic        := grepl("^L01",           atc_koodi)]
meds_win[, immunosuppressant := grepl("^L04",           atc_koodi)]

# Any prescription of the class -> flagged, per operation
# n_med_records counts current-EHR prescriptions active within 30 d (polypharmacy proxy)
med_flags <- meds_win[, .(
    anticoagulant     = any(anticoagulant,     na.rm = TRUE),
    antiplatelet      = any(antiplatelet,      na.rm = TRUE),
    nsaid             = any(nsaid,             na.rm = TRUE),
    corticosteroid    = any(corticosteroid,    na.rm = TRUE),
    coag_factor       = any(coag_factor,       na.rm = TRUE),
    hemostatic        = any(hemostatic,        na.rm = TRUE),
    iron_oral         = any(iron_oral,         na.rm = TRUE),
    iron_iv           = any(iron_iv,           na.rm = TRUE),
    esa               = any(esa,               na.rm = TRUE),
    cytostatic        = any(cytostatic,        na.rm = TRUE),
    immunosuppressant = any(immunosuppressant, na.rm = TRUE),
    has_aspirin       = any(aspirin,           na.rm = TRUE),
    has_p2y12         = any(p2y12,             na.rm = TRUE),
    n_med_records     = sum(src == "ap" & active_30d)
), by = .(henkilotunnus, potilas_saliin_hetki = op_time)]

rm(meds_win, med_or); gc()

# Dual antiplatelet therapy == aspirin and a P2Y12 inhibitor both active
med_flags[, dapt := has_aspirin & has_p2y12]
med_flags[, `:=`(has_aspirin = NULL, has_p2y12 = NULL)]

# 0/1 integers instead of TRUE/FALSE
med_flag_cols <- c(
    "anticoagulant", "antiplatelet", "nsaid", "corticosteroid", "coag_factor", "hemostatic",
    "iron_oral", "iron_iv", "esa", "cytostatic", "immunosuppressant", "dapt"
)
med_flags[, (med_flag_cols) := lapply(.SD, as.integer), .SDcols = med_flag_cols]

# Back to every operation, 0 where no prescriptions
med_flags <- med_flags[
    surg_dates[, .(henkilotunnus, potilas_saliin_hetki)],
    on = .(henkilotunnus, potilas_saliin_hetki)
]
setnafill(med_flags, fill = 0L, cols = c(med_flag_cols, "n_med_records"))

rm(surg_dates, cohort); gc()

# ---------------------------------------------------------------------------
# Combine all patient variables
# ---------------------------------------------------------------------------
# One row per operation, every predictor block joined on patient and OR entry time
patient_vars <- Reduce(
    function(x, y) merge(x, y, by = c("henkilotunnus", "potilas_saliin_hetki"), all.x = TRUE),
    list(
        age_sex, bv, hb_entries, trom_entries,
        creat_entries, crp_entries, mcv_entries, inr_entries,
        alb_entries, wbc_entries, ferr_entries,
        bl_histories, elxh, asa, med_flags
    )
)

rm(
    age_sex, bv, hb_entries, trom_entries,
    creat_entries, crp_entries, mcv_entries, inr_entries,
    alb_entries, wbc_entries, ferr_entries,
    bl_histories, elxh, asa, med_flags
); gc()

# ---------------------------------------------------------------------------
# Build modelling dataset
# ---------------------------------------------------------------------------
# Continuous predictors, z-scored at the end 
# (everything except ids, categoricals and 0/1 flags)
cols_to_scale <- setdiff(
    names(patient_vars),
    c("henkilotunnus", "potilas_saliin_hetki", "sukupuoli_selite", "asa", med_flag_cols)
)

dtbin <- (
    merge(dt, patient_vars, by = c("henkilotunnus", "potilas_saliin_hetki"), all.x = TRUE)
    # Outcome: at least one RBC unit within the outcome window
    [, transfused      := rbc_count > 0]
    [, asa             := as.factor(asa)]
    # Procedure group is the first three characters of the NCSP code, e.g. JJB
    [, procedure_group := substr(paatoimenpide_koodi, 1, 3)]
    # Exclusions: missing or non-binary sex, ASA class V (likely a coding error for elective surgery)
    [sukupuoli_selite %in% c("Mies", "Nainen")]
    # is.na() first, since a plain asa != "5" would drop the missing-ASA rows too!
    [is.na(asa) | asa != "5"]
    [, asa := droplevels(asa)]
    [, `:=`(
        potilas_salista_hetki = NULL,
        rbc_count             = NULL,
        proc_n                = NULL,
        paatoimenpide_koodi   = NULL,
        paatoimenpide_selite  = NULL,
        asa_luokka_koodi      = NULL
    )]
    # Missing ASA is its own category
    [is.na(asa), asa := "Missing"]
    # No recorded earlier bleeding == no blood loss
    [is.na(total_prev_bl), total_prev_bl := 0]
    [is.na(n_prev_bl_eps), n_prev_bl_eps := 0]
)
# (patient id and OR entry time stay until after the split)

rm(dt, patient_vars); gc()

# ---------------------------------------------------------------------------
# Missing values
# ---------------------------------------------------------------------------
# BMI and blood volume: regression imputation (missing assumed unrelated to the outcome)
impute_predictors <- c("age", "sukupuoli_selite", "asa", "procedure_group")

message("Imputing bmi and bv_litres by regression ...")
impute_by_regression(dtbin, "bmi",       impute_predictors)
impute_by_regression(dtbin, "bv_litres", impute_predictors)

# Labs are left missing on purpose: whether a test was ordered is informative in itself
# CatBoost treats NA as smaller than any observed value (nan_mode = "Min"), so it can split on missingness

# INR converted from TT% (same draw) is a unit conversion, not imputation, so it counts as measured
# No maximum age for a lab result: an old result != no result, and the days-since-drawn columns let
# the model judge staleness per analyte
lab_cols <- c(
    "latest_hb", "latest_tromb", "latest_creatinine", "latest_crp", "latest_mcv",
    "latest_inr", "latest_albumin", "latest_wbc", "latest_ferritin"
)
for (col in lab_cols) {
    message(sprintf("  %s: %d missing (%.1f%%)", col, sum(is.na(dtbin[[col]])), 100 * mean(is.na(dtbin[[col]]))))
}

days_cols <- paste0(lab_cols, "_days_since_drawn")
for (col in days_cols) {
    message(sprintf(
        "  %s: %d missing (%.1f%%), median %.0f d",
        col,
        sum(is.na(dtbin[[col]])),
        100 * mean(is.na(dtbin[[col]])),
        median(dtbin[[col]], na.rm = TRUE)
    ))
}

# ---------------------------------------------------------------------------
# Train/test split
# ---------------------------------------------------------------------------
# By patient, so all of a patient's procedures land in the same set (no leakage!)
# Stratified by whether the patient was ever transfused
set.seed(123)

# One row per patient, in random order
pt_summary <- dtbin[, .(
    n_rows         = .N,
    any_transfused = any(transfused)
), by = henkilotunnus][sample(.N)]

# Within each stratum, take patients from the top until TEST_FRAC of its procedures is reached
test_ids <- pt_summary[, {
    n_target <- TEST_FRAC * sum(n_rows)
    sel      <- which(cumsum(n_rows) <= n_target)
    # At least one patient per stratum (just in case)
    if (length(sel) == 0L) sel <- 1L
    .(henkilotunnus = henkilotunnus[sel])
}, by = any_transfused]$henkilotunnus

dtbin[, set := fifelse(henkilotunnus %in% test_ids, "test", "train")]

# Pseudonymous patient number, kept after the ids are dropped, for ce_03's cluster bootstrap
# (not a predictor, excluded in ce_02 and ce_02a)
dtbin[, patient_idx := .GRP, by = henkilotunnus]

# Patient-grouped tuning folds for ce_02a and the class-weighting ablation, training rows only
# Same stratification as the split. Own seed.
set.seed(20260924)
train_pts <- pt_summary[!henkilotunnus %in% test_ids, .(henkilotunnus, any_transfused)][sample(.N)]

# Deal patients round-robin into folds 1..N_CV_FOLDS within each stratum
train_pts[, cv_fold := (seq_len(.N) - 1L) %% N_CV_FOLDS + 1L, by = any_transfused]

# Test rows keep cv_fold = NA
dtbin[train_pts, cv_fold := i.cv_fold, on = "henkilotunnus"]
rm(train_pts)

# ---------------------------------------------------------------------------
# Cohort metadata
# For ce_04 and ce_05 (Table 1, Methods and Results numbers)
# ---------------------------------------------------------------------------
# Needs the raw (unscaled) Hb and the ids, so computed before both are gone
cohort_meta <- list(
    n_patients  = uniqueN(dtbin$henkilotunnus),
    study_start = as.IDate(min(dtbin$potilas_saliin_hetki)),
    study_end   = as.IDate(max(dtbin$potilas_saliin_hetki)),
    hb_med_raw  = median(dtbin$latest_hb, na.rm = TRUE),
    hb_q25_raw  = as.numeric(quantile(dtbin$latest_hb, 0.25, na.rm = TRUE)),
    hb_q75_raw  = as.numeric(quantile(dtbin$latest_hb, 0.75, na.rm = TRUE))
)

# Share of the cohort's INR values that were derived from TT% (Supplementary Methods 3)
inr_src <- merge(
    dtbin[!is.na(latest_inr), .(henkilotunnus, potilas_saliin_hetki)],
    inr_derived_entries,
    by    = c("henkilotunnus", "potilas_saliin_hetki"),
    all.x = TRUE
)
cohort_meta$n_inr_measured <- nrow(inr_src)
cohort_meta$n_inr_derived  <- sum(inr_src$latest_inr_derived == 1, na.rm = TRUE)
message(sprintf(
    "INR: %d of %d qualifying values (%.1f%%) derived from TT%%",
    cohort_meta$n_inr_derived,
    cohort_meta$n_inr_measured,
    100 * cohort_meta$n_inr_derived / cohort_meta$n_inr_measured
))
rm(inr_src, inr_derived_entries)

# Identifiers aren't needed after the split
dtbin[, `:=`(henkilotunnus = NULL, potilas_saliin_hetki = NULL)]

# ---------------------------------------------------------------------------
# Scaling
# ---------------------------------------------------------------------------
# z-scores from training-set means and SDs only (no test information leaks in)
# (doesn't change tree splits, only the units)
train_means <- dtbin[set == "train", lapply(.SD, mean, na.rm = TRUE), .SDcols = cols_to_scale]
train_sds   <- dtbin[set == "train", lapply(.SD, sd,   na.rm = TRUE), .SDcols = cols_to_scale]

# Kept so ce_04 and ce_05 can turn values back into clinical units (x * sd + mean)
cohort_meta$scale_means <- as.list(train_means)
cohort_meta$scale_sds   <- as.list(train_sds)
saveRDS(cohort_meta, file.path(objpath, "cohort_meta.rds"))
message(
    "Saved cohort_meta.rds: ", cohort_meta$n_patients, " patients, study ",
    cohort_meta$study_start, " to ", cohort_meta$study_end
)

for (col in cols_to_scale) {
    # A column constant in train has SD 0 (or NA), and would turn into NaN (just in case)
    if (is.na(train_sds[[col]]) || train_sds[[col]] == 0) {
        message(sprintf("  skipping z-score for %s: zero/NA train SD", col))
        next
    }
    set(dtbin, j = col, value = (dtbin[[col]] - train_means[[col]]) / train_sds[[col]])
}

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
saveRDS(dtbin, file.path(objpath, "dtbin.rds"))
message(
    "Saved dtbin.rds: ", nrow(dtbin), " rows",
    " | train: ", sum(dtbin$set == "train"),
    " | test: ",  sum(dtbin$set == "test")
)
