# Computes the remaining numbers quoted in the manuscript, and Table 1.
#
# Prints them and saves them as manuscript_stats.rds and a Table 1 CSV.

library(data.table)

model_tag <- "final"

# Paths are relative to the project root
objpath <- file.path("results", "model_objects")

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
# dtbin is the modelling dataset from ce_01, train/test split in `set`
dtbin   <- readRDS(file.path(objpath, "dtbin.rds"))
test_dt <- dtbin[set == "test"]

# Counts and raw values saved by ce_01 before the ids were dropped and values z-scored
meta <- readRDS(file.path(objpath, "cohort_meta.rds"))

# Patient-level predictors, from ce_02
fe_predictors <- readRDS(file.path(objpath, paste0("fe_predictors_", model_tag, ".rds")))

# Per-surgeon effects (Figure 4) and Figure 3's surgeon rows, from ce_03
surgeon_shap   <- fread(file.path(objpath, paste0("surgeon_shap_", model_tag, ".csv")))
fig3_surgeons  <- fread(file.path(objpath, paste0("figure3_surgeon_rates_", model_tag, ".csv")))

# ---------------------------------------------------------------------------
# Cohort
# ---------------------------------------------------------------------------
n_procedures       <- nrow(dtbin)
n_surgeons         <- uniqueN(dtbin$toimenpide_tekija_id)
n_procedure_groups <- uniqueN(dtbin$procedure_group)
n_test             <- nrow(test_dt)

# Transfusion rate with an exact binomial 95% CI
n_transfused <- sum(dtbin$transfused)
ci           <- binom.test(n_transfused, n_procedures)$conf.int
transf_ci_lo <- round(ci[1] * 100, 1)
transf_ci_hi <- round(ci[2] * 100, 1)

# Share without a qualifying preoperative Hb
pct_hb_missing <- round(mean(is.na(dtbin$latest_hb)) * 100, 1)

n_patients  <- meta$n_patients
study_start <- format(meta$study_start, "%B %Y")
study_end   <- format(meta$study_end,   "%B %Y")

# Hb in g/L (dtbin's own Hb is z-scored)
hb_med_raw <- round(meta$hb_med_raw, 0)
hb_q25_raw <- round(meta$hb_q25_raw, 0)
hb_q75_raw <- round(meta$hb_q75_raw, 0)

# ---------------------------------------------------------------------------
# Surgeons
# ---------------------------------------------------------------------------
# Procedures per surgeon over the whole cohort, median (IQR)
surgeon_vol <- dtbin[, .N, by = toimenpide_tekija_id]$N
surg_vol_q  <- quantile(surgeon_vol, c(0.25, 0.5, 0.75))

# Surgeons seen only in the test set
# (CatBoost gives them the prior, i.e. no individual surgeon signal)
train_surgeons         <- unique(dtbin[set == "train", toimenpide_tekija_id])
test_only_surg         <- setdiff(unique(test_dt$toimenpide_tekija_id), train_surgeons)
n_test_only_surgeons   <- length(test_only_surg)
n_test_only_surg_procs <- test_dt[toimenpide_tekija_id %in% test_only_surg, .N]

# Surgeon-effect categories of Figure 4 (eligibility and 1 SD grouping done in ce_03)
n_surgeons_shap     <- nrow(surgeon_shap)
n_above_avg         <- sum(surgeon_shap$status == "Above average")
n_below_avg         <- sum(surgeon_shap$status == "Below average")
n_above_below_total <- n_above_avg + n_below_avg
n_within_1sd        <- n_surgeons_shap - n_above_below_total

# Surgeons shown in Figure 3 (legend count)
n_surgeons_plot <- uniqueN(fig3_surgeons$toimenpide_tekija_id)

# ---------------------------------------------------------------------------
# Predictors
# ---------------------------------------------------------------------------
# Share of INR values derived from TT% (Supplementary Methods 3)
n_inr_measured  <- meta$n_inr_measured
n_inr_derived   <- meta$n_inr_derived
pct_inr_derived <- round(100 * n_inr_derived / n_inr_measured, 1)

# Full model features: patient predictors plus surgeon and procedure group
n_features_full <- length(fe_predictors) + 2L

# ce_01 only creates columns for Elixhauser categories that occur in a look-back window,
# so the number of categories present comes from the predictor list (Table S2)
elix_all <- c(
    "chf", "carit", "valv", "pcd", "pvd", "hypunc", "hypc", "para",
    "ond", "cpd", "diabunc", "diabc", "hypothy", "rf", "ld", "pud",
    "aids", "lymph", "metacanc", "solidtum", "rheumd", "coag", "obes",
    "wloss", "fed", "blane", "dane", "alcohol", "drug", "psycho", "depre"
)
elix_present   <- intersect(elix_all, fe_predictors)
elix_absent    <- setdiff(elix_all, fe_predictors)
n_elix_present <- length(elix_present)

# ---------------------------------------------------------------------------
# Software versions (Methods)
# ---------------------------------------------------------------------------
# The "key packages" named in the Methods
r_version          <- paste(R.version$major, R.version$minor, sep = ".")
data_table_version <- as.character(packageVersion("data.table"))
catboost_version   <- as.character(packageVersion("catboost"))
prroc_version      <- as.character(packageVersion("PRROC"))

# ---------------------------------------------------------------------------
# Table 1
# ---------------------------------------------------------------------------
# Median (IQR) as "m (q1-q3)"
fmt_median_iqr <- function(x, digits = 0L) {
    q   <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
    fmt <- sprintf("%%.%df", digits)
    sprintf("%s (%s-%s)", sprintf(fmt, q[2]), sprintf(fmt, q[1]), sprintf(fmt, q[3]))
}

# Count with percentage, e.g. "1,234 (5.6%)"
fmt_n_pct <- function(n, N) {
    sprintf("%s (%.1f%%)", format(n, big.mark = ","), 100 * n / N)
}

# dtbin's continuous predictors are z-scored (ce_01), back to clinical units with ce_01's train mean and SD
# Columns that were never scaled come back unchanged
unscale_col <- function(x, var) {
    m <- meta$scale_means[[var]]
    s <- meta$scale_sds[[var]]
    if (is.null(m) || is.null(s)) return(x)
    x * s + m
}

# Table 1 columns: all procedures, transfused, not transfused
strata <- list(
    All              = dtbin,
    Transfused       = dtbin[transfused == 1],
    `Not transfused` = dtbin[transfused == 0]
)
n_by_stratum <- vapply(strata, nrow, integer(1))

# One Table 1 row with the median (IQR) of a continuous variable per column
row_median_iqr <- function(var, label, digits = 0L) {
    data.table(
        Variable         = label,
        All              = fmt_median_iqr(unscale_col(strata$All[[var]], var), digits),
        Transfused       = fmt_median_iqr(unscale_col(strata$Transfused[[var]], var), digits),
        `Not transfused` = fmt_median_iqr(unscale_col(strata$`Not transfused`[[var]], var), digits)
    )
}

# One row per level of a categorical variable, n (%) per column
row_n_pct_by_level <- function(var, label_prefix) {
    lvls <- sort(unique(dtbin[[var]]))
    rbindlist(lapply(
        lvls,
        function(lv) {
            data.table(
                Variable         = paste0(label_prefix, " = ", lv),
                All              = fmt_n_pct(sum(strata$All[[var]] == lv, na.rm = TRUE), n_by_stratum["All"]),
                Transfused       = fmt_n_pct(sum(strata$Transfused[[var]] == lv, na.rm = TRUE), n_by_stratum["Transfused"]),
                `Not transfused` = fmt_n_pct(sum(strata$`Not transfused`[[var]] == lv, na.rm = TRUE), n_by_stratum["Not transfused"])
            )
        }
    ))
}

# Five procedure groups with the most transfused procedures
top5_groups <- dtbin[transfused == 1, .N, by = procedure_group][order(-N)][1:5, procedure_group]

table1 <- rbindlist(list(
    data.table(
        Variable         = "N procedures",
        All              = format(n_by_stratum["All"], big.mark = ","),
        Transfused       = format(n_by_stratum["Transfused"], big.mark = ","),
        `Not transfused` = format(n_by_stratum["Not transfused"], big.mark = ",")
    ),
    row_median_iqr("age", "Age, years, median (IQR)"),
    row_n_pct_by_level("sukupuoli_selite", "Sex"),
    row_n_pct_by_level("asa", "ASA class"),
    row_median_iqr("bmi", "BMI, median (IQR)"),
    row_median_iqr("bv_litres", "Estimated blood volume, L, median (IQR)", digits = 1L),
    row_median_iqr("latest_hb", "Pre-operative haemoglobin, g/L, median (IQR)"),
    row_median_iqr("latest_tromb", "Platelet count, ×10⁹/L, median (IQR)"),
    data.table(
        Variable         = "No qualifying pre-operative Hb, n (%)",
        All              = fmt_n_pct(sum(is.na(strata$All$latest_hb)), n_by_stratum["All"]),
        Transfused       = fmt_n_pct(sum(is.na(strata$Transfused$latest_hb)), n_by_stratum["Transfused"]),
        `Not transfused` = fmt_n_pct(sum(is.na(strata$`Not transfused`$latest_hb)), n_by_stratum["Not transfused"])
    ),
    # Rate only makes sense for all procedures, so the other two cells stay empty
    data.table(
        Variable         = "Transfusion rate, n (%, 95% CI)",
        All              = sprintf(
            "%s (%.1f%%, %.1f-%.1f%%)",
            format(n_transfused, big.mark = ","),
            100 * n_transfused / n_procedures,
            transf_ci_lo,
            transf_ci_hi
        ),
        Transfused       = "",
        `Not transfused` = ""
    ),
    rbindlist(lapply(
        top5_groups,
        function(pg) {
            data.table(
                Variable         = paste0("Procedure group ", pg, ", n (%)"),
                All              = fmt_n_pct(sum(strata$All$procedure_group == pg), n_by_stratum["All"]),
                Transfused       = fmt_n_pct(sum(strata$Transfused$procedure_group == pg), n_by_stratum["Transfused"]),
                `Not transfused` = fmt_n_pct(sum(strata$`Not transfused`$procedure_group == pg), n_by_stratum["Not transfused"])
            )
        }
    ))
))

print(table1)
fwrite(table1, file.path(objpath, "table1_cohort_characteristics.csv"))
saveRDS(table1, file.path(objpath, "table1_cohort_characteristics.rds"))

# ---------------------------------------------------------------------------
# Print
# ---------------------------------------------------------------------------
cat("\n===== MANUSCRIPT STATISTICS =====\n\n")
cat(sprintf("N_PROCEDURES        = %s\n", format(n_procedures, big.mark = ",")))
cat(sprintf("N_PATIENTS          = %s\n", format(n_patients, big.mark = ",")))
cat(sprintf("N_SURGEONS          = %s\n", format(n_surgeons, big.mark = ",")))
cat(sprintf("N_PROCEDURE_GROUPS  = %s\n", format(n_procedure_groups, big.mark = ",")))
cat(sprintf("STUDY_PERIOD        = %s to %s\n", study_start, study_end))
cat(sprintf("N_TRANSFUSED        = %s\n", format(n_transfused, big.mark = ",")))
cat(sprintf("TRANSF_CI_LO        = %.1f%%\n", transf_ci_lo))
cat(sprintf("TRANSF_CI_HI        = %.1f%%\n", transf_ci_hi))
cat(sprintf("HB_MED_RAW          = %d g/L\n", hb_med_raw))
cat(sprintf("HB_Q25_RAW          = %d g/L\n", hb_q25_raw))
cat(sprintf("HB_Q75_RAW          = %d g/L\n", hb_q75_raw))
cat(sprintf("PCT_HB_MISSING      = %.1f%%\n", pct_hb_missing))
cat(sprintf("N_TEST              = %s\n", format(n_test, big.mark = ",")))
cat(sprintf("N_SURGEONS_SHAP     = %d\n", n_surgeons_shap))
cat(sprintf("N_ABOVE_AVG         = %d\n", n_above_avg))
cat(sprintf("N_BELOW_AVG         = %d\n", n_below_avg))
cat(sprintf("N_ABOVE_BELOW_TOTAL = %d\n", n_above_below_total))
cat(sprintf("N_WITHIN_1SD        = %d\n", n_within_1sd))
cat(sprintf("SURG_VOL_MEDIAN_IQR = %.0f (%.0f-%.0f) procedures per surgeon\n", surg_vol_q[2], surg_vol_q[1], surg_vol_q[3]))
cat(sprintf("N_TEST_ONLY_SURG    = %d  (%s test procedures)\n", n_test_only_surgeons, format(n_test_only_surg_procs, big.mark = ",")))
cat(sprintf(
    "INR_DERIVED         = %s of %s (%.1f%%)\n",
    format(n_inr_derived, big.mark = ","),
    format(n_inr_measured, big.mark = ","),
    pct_inr_derived
))
cat(sprintf("N_FEATURES_FULL     = %d\n", n_features_full))
cat(sprintf(
    "N_ELIX_PRESENT      = %d of 31  (absent: %s)\n",
    n_elix_present,
    if (length(elix_absent)) paste(elix_absent, collapse = ", ") else "none"
))
cat(sprintf("N_SURGEONS_PLOT     = %d\n", n_surgeons_plot))
cat(sprintf("R_VERSION           = %s\n", r_version))
cat(sprintf("DATA_TABLE_VERSION  = %s\n", data_table_version))
cat(sprintf("CATBOOST_VERSION    = %s\n", catboost_version))
cat(sprintf("PRROC_VERSION       = %s\n", prroc_version))
cat("\n=================================\n")

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
saveRDS(
    list(
        n_procedures           = n_procedures,
        n_patients             = n_patients,
        n_surgeons             = n_surgeons,
        n_procedure_groups     = n_procedure_groups,
        study_start            = meta$study_start,
        study_end              = meta$study_end,
        n_transfused           = n_transfused,
        transf_ci_lo           = transf_ci_lo,
        transf_ci_hi           = transf_ci_hi,
        hb_med_raw             = hb_med_raw,
        hb_q25_raw             = hb_q25_raw,
        hb_q75_raw             = hb_q75_raw,
        pct_hb_missing         = pct_hb_missing,
        n_test                 = n_test,
        n_surgeons_shap        = n_surgeons_shap,
        n_above_avg            = n_above_avg,
        n_below_avg            = n_below_avg,
        n_above_below_total    = n_above_below_total,
        n_within_1sd           = n_within_1sd,
        surg_vol_q             = surg_vol_q,
        n_test_only_surgeons   = n_test_only_surgeons,
        n_test_only_surg_procs = n_test_only_surg_procs,
        n_inr_measured         = n_inr_measured,
        n_inr_derived          = n_inr_derived,
        n_features_full        = n_features_full,
        elix_present           = elix_present,
        n_surgeons_plot        = n_surgeons_plot,
        r_version              = r_version,
        data_table_version     = data_table_version,
        catboost_version       = catboost_version,
        prroc_version          = prroc_version
    ),
    file.path(objpath, "manuscript_stats.rds")
)
message("Saved manuscript_stats.rds")
