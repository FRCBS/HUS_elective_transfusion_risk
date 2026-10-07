# Computes the supplement's reporting numbers (TRIPOD+AI / RECORD items).
#
# Participant flow, imputation shares, split sizes, train vs test characteristics,
# hyperparameters, calibration, subgroups and procedure-group performance.
# Read-only: refits nothing, prints everything into a log for the supplement.

library(data.table)
library(PRROC)
library(catboost)
library(parallel)

model_tag <- "final"
n_cores   <- detectCores()

PRED_HORIZON_DAYS <- 7L     # same as ce_01
MIN_PROC_N        <- 50L    # same as ce_01
N_BOOT            <- 2000L  # patient-cluster bootstrap resamples
N_BOOT_GROUP      <- 1000L  # resamples within each procedure group (Table S8)

# Paths are relative to the project root
datapath    <- "data"
resultspath <- "results"
objpath     <- file.path(resultspath, "model_objects")
outdir      <- file.path(resultspath, "reporting_supplement")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Everything printed also goes to this log (the file to bring back for the supplement)
log_file <- file.path(outdir, paste0("ce_05_output_", model_tag, ".txt"))
sink(log_file, split = TRUE)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
# Section banner in the log
section <- function(title) {
    cat("\n\n", strrep("=", 78), "\n", title, "\n", strrep("=", 78), "\n", sep = "")
}

# Tables go to outdir as CSV and into the log
save_tab <- function(dt, name) {
    fwrite(dt, file.path(outdir, paste0(name, "_", model_tag, ".csv")))
    print(dt, nrows = 200L)
}

# Count with percentage, e.g. "1,234 (5.6%)"
fmt_n_pct <- function(n, N) {
    sprintf("%s (%.1f%%)", format(n, big.mark = ","), 100 * n / N)
}

# Median (IQR) as "m (q1-q3)"
fmt_med_iqr <- function(x, digits = 0L) {
    q   <- quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
    fmt <- paste0("%.", digits, "f")
    sprintf("%s (%s-%s)", sprintf(fmt, q[2]), sprintf(fmt, q[1]), sprintf(fmt, q[3]))
}

# dtbin is the modelling dataset from ce_01, with z-scored continuous predictors
dtbin         <- readRDS(file.path(objpath, "dtbin.rds"))
meta          <- readRDS(file.path(objpath, "cohort_meta.rds"))
fe_predictors <- readRDS(file.path(objpath, paste0("fe_predictors_", model_tag, ".rds")))

# Back to clinical units with ce_01's train mean and SD (unscaled columns come back unchanged)
unscale_col <- function(x, var) {
    m <- meta$scale_means[[var]]
    s <- meta$scale_sds[[var]]
    if (is.null(m) || is.null(s)) return(x)
    x * s + m
}

# ---------------------------------------------------------------------------
# 1. Participant flow (Figure S2, Supplementary Methods 1)
# ---------------------------------------------------------------------------
section("1. PARTICIPANT FLOW")

# ce_01's exclusion chain replayed on the raw operation table, without the outcome join
# (an independent check: the final count is compared with dtbin below)
leik <- readRDS(file.path(datapath, "apotti_leikkaus_operaatio_henkilosto.rds"))

# One row per operation (ce_01 keeps the most senior staff member; any row gives the same count)
n_rows_raw <- nrow(leik)
keep_idx   <- leik[, .I[if (is.na(operaatio_id[1L])) seq_len(.N) else 1L], by = operaatio_id]$V1
ops_all    <- leik[sort(keep_idx)]
rm(leik, keep_idx); gc()

# Procedure and patient counts after each step
flow <- list()
add_step <- function(step, dt) {
    flow[[length(flow) + 1L]] <<- data.table(
        step         = step,
        n_procedures = nrow(dt),
        n_patients   = uniqueN(dt$henkilotunnus, na.rm = TRUE)
    )
}

add_step("Surgical operations in source table (after collapsing personnel rows)", ops_all)

s1 <- ops_all[leikkaus_kiireellisyysluokka_selite %in% c("LPI", "LPII", "LPIII")]
add_step("Elective (national urgency codes)", s1)

s2 <- s1[!is.na(potilas_saliin_hetki) & !is.na(potilas_salista_hetki) & !is.na(paatoimenpide_koodi)]
add_step("OR entry time, exit time and primary NCSP code recorded", s2)

# ce_01 does this inside the outcome join; here explicitly, keeping the earliest OR exit
setorder(s2, henkilotunnus, potilas_saliin_hetki, potilas_salista_hetki)
s3 <- unique(s2, by = c("henkilotunnus", "potilas_saliin_hetki"))
add_step("One record per patient and OR entry time", s3)

s4 <- s3[, proc_n := .N, by = paatoimenpide_koodi][proc_n >= MIN_PROC_N]
add_step("Primary procedure code with >= 50 cases in study period", s4)

asiakas <- readRDS(file.path(datapath, "asiakas_out.rds"))
sex_ok  <- unique(
    asiakas[sukupuoli_selite %in% c("Mies", "Nainen"), .(henkilotunnus, sukupuoli_selite, syntymaaika)],
    by = "henkilotunnus"
)
rm(asiakas); gc()
s5 <- merge(s4, sex_ok, by = "henkilotunnus")
add_step("Sex recorded as male or female", s5)

s6 <- s5[is.na(asa_luokka_koodi) | as.character(asa_luokka_koodi) != "5"]
add_step("ASA class not V (final cohort)", s6)

# Excluded counts are the drop from the previous step
flow_dt <- rbindlist(flow)
flow_dt[, excluded_procedures := shift(n_procedures) - n_procedures]
flow_dt[, excluded_patients   := shift(n_patients)   - n_patients]
cat("\nRows in raw table before collapsing multi-personnel rows: ", format(n_rows_raw, big.mark = ","), "\n", sep = "")
cat(
    "Non-elective operations with missing urgency code (within the elective exclusion): ",
    format(ops_all[is.na(leikkaus_kiireellisyysluokka_selite), .N], big.mark = ","), "\n",
    sep = ""
)
save_tab(flow_dt, "flow")

# Does the replayed chain arrive at dtbin? (a mismatch is reported, not hidden!)
cat(sprintf(
    "\nReconciliation: replicated final cohort %d procedures / %d patients; dtbin %d procedures / cohort_meta %d patients\n",
    nrow(s6),
    uniqueN(s6$henkilotunnus),
    nrow(dtbin),
    meta$n_patients
))
if (nrow(s6) != nrow(dtbin) || uniqueN(s6$henkilotunnus) != meta$n_patients) {
    cat(
        "WARNING: replicated flow does not match dtbin exactly; the difference is",
        nrow(s6) - nrow(dtbin), "procedures. Report dtbin counts as final and",
        "describe the step counts as approximate.\n"
    )
}

# Adults only? (the data lake extract was restricted to adults)
s6[, age_check := as.numeric(difftime(as.IDate(potilas_saliin_hetki), syntymaaika, units = "weeks")) / 52.143]
cat(sprintf("Age at surgery < 18 y in final cohort: %d; missing age: %d\n", s6[age_check < 18, .N], s6[is.na(age_check), .N]))

cat("\nFinal cohort outcome (dtbin):\n")
print(dtbin[, .(procedures = .N, transfused = sum(transfused), not_transfused = sum(!transfused))])

# Operating units in the final cohort (TRIPOD+AI 6a)
cat(sprintf(
    "\nHospitals (leikkaava_toimipiste_selite) in final cohort: %d distinct; %d procedures with no site recorded\n",
    uniqueN(s6$leikkaava_toimipiste_selite, na.rm = TRUE),
    s6[is.na(leikkaava_toimipiste_selite), .N]
))
site_dt <- s6[, .(procedures = .N, patients = uniqueN(henkilotunnus)), by = leikkaava_toimipiste_selite][order(-procedures)]
site_dt[, pct_procedures := round(100 * procedures / nrow(s6), 1)]
save_tab(site_dt, "hospitals")

# ---------------------------------------------------------------------------
# 2. BMI and blood volume imputation (Supplementary Methods 2)
# ---------------------------------------------------------------------------
section("2. BMI AND BLOOD VOLUME IMPUTATION")

# dtbin only has the imputed values, so the missing share is recounted with ce_01's selection rule
keys <- s6[, .(
    henkilotunnus,
    potilas_saliin_hetki,
    sukupuoli_selite,
    surg_date       = as.IDate(potilas_saliin_hetki),
    age             = age_check,
    asa             = fifelse(is.na(asa_luokka_koodi), "Missing", as.character(asa_luokka_koodi)),
    procedure_group = substr(paatoimenpide_koodi, 1, 3)
)]

bmi <- readRDS(file.path(datapath, "mittaukset_bmi.rds"))[henkilotunnus %in% keys$henkilotunnus]

# Is the source BMI weight / height^2? (backs the "calculated from height and weight" wording)
bmi_chk <- bmi[!is.na(bmi) & !is.na(paino) & !is.na(pituus) & pituus > 0]
bmi_chk[, bmi_calc := paino / (pituus / 100)^2]
cat(sprintf(
    "Source BMI rows with height and weight: %d; |BMI - weight/height^2| < 0.5 in %.1f%%\n",
    nrow(bmi_chk),
    100 * mean(abs(bmi_chk$bmi - bmi_chk$bmi_calc) < 0.5)
))
cat(sprintf(
    "Source BMI rows: %d total; BMI missing but height and weight present: %d\n",
    nrow(bmi),
    bmi[is.na(bmi) & !is.na(paino) & !is.na(pituus), .N]
))
rm(bmi_chk)

# Most recent measurement at least PRED_HORIZON_DAYS before surgery (as in ce_01)
bmi_preop <- merge(keys[, .(henkilotunnus, potilas_saliin_hetki, surg_date)], bmi, by = "henkilotunnus", allow.cartesian = TRUE)
bmi_preop[, diff := difftime(surg_date, mittaus_pvm, units = "days")]
bmi_preop <- bmi_preop[diff >= PRED_HORIZON_DAYS]
setorder(bmi_preop, henkilotunnus, potilas_saliin_hetki, diff)
bmi_preop <- unique(bmi_preop, by = c("henkilotunnus", "potilas_saliin_hetki"))
rm(bmi); gc()

# Blood volume with ce_01's formulas (Lemmens-Bernstein-Brodsky above BMI 30, Nadler below)
bv_chk <- merge(
    keys,
    bmi_preop[, .(henkilotunnus, potilas_saliin_hetki, bmi, paino, pituus)],
    by    = c("henkilotunnus", "potilas_saliin_hetki"),
    all.x = TRUE
)
bv_chk[, bv_litres := fcase(
    bmi > 30,                                 paino * 70 / sqrt(bmi / 22) / 1000,
    bmi <= 30 & sukupuoli_selite == "Mies",   0.3669 * (pituus / 100)^3 + 0.03219 * paino + 0.6041,
    bmi <= 30 & sukupuoli_selite == "Nainen", 0.3561 * (pituus / 100)^3 + 0.03308 * paino + 0.1833
)]

imp_dt <- data.table(
    variable  = c("bmi", "bv_litres"),
    n_cohort  = nrow(bv_chk),
    n_missing = c(sum(is.na(bv_chk$bmi)), sum(is.na(bv_chk$bv_litres)))
)
imp_dt[, pct_imputed := round(100 * n_missing / n_cohort, 1)]

# R2 of the imputation regression as fitted in ce_01 (full cohort, before the split)
imp_dt[, r2_imputation_model := vapply(
    variable,
    function(v) {
        fml <- as.formula(paste(v, "~ age + sukupuoli_selite + asa + procedure_group"))
        summary(lm(fml, data = bv_chk[!is.na(get(v))]))$r.squared
    },
    numeric(1)
)]
save_tab(imp_dt, "bmi_bv_imputation")
rm(bmi_preop, bv_chk, keys); gc()

# ---------------------------------------------------------------------------
# 3. Split, fold and ablation sizes (Supplementary Methods 5)
# ---------------------------------------------------------------------------
section("3. TRAIN / TEST / CV-FOLD SIZES AND EVENTS")

split_dt <- dtbin[, .(
    procedures       = .N,
    patients         = uniqueN(patient_idx),
    transfused       = sum(transfused),
    pct_transfused   = round(100 * mean(transfused), 2),
    surgeons         = uniqueN(toimenpide_tekija_id),
    procedure_groups = uniqueN(procedure_group)
), keyby = set]
save_tab(split_dt, "split_sizes")

# Tuning folds within the training set
fold_dt <- dtbin[set == "train", .(
    procedures = .N,
    patients   = uniqueN(patient_idx),
    transfused = sum(transfused)
), keyby = cv_fold]
save_tab(fold_dt, "cv_fold_sizes")
cat("Class-weighting ablation: holdout = fold 1; ablation training/tuning = folds 2-5.\n")

# Events per feature (sample size justification)
n_features_full <- length(fe_predictors) + 2L
cat(sprintf(
    "\nCandidate features: %d patient-level + surgeon + procedure group = %d in the full model\n",
    length(fe_predictors),
    n_features_full
))
cat(sprintf(
    "Training events per full-model feature: %.1f; test events: %d\n",
    split_dt["train", transfused] / n_features_full,
    split_dt["test", transfused]
))

# ---------------------------------------------------------------------------
# 4. Train vs test characteristics (Table S9)
# ---------------------------------------------------------------------------
section("4. TRAIN VS TEST CHARACTERISTICS")

sets <- list(Train = dtbin[set == "train"], Test = dtbin[set == "test"])

# Median (IQR) row, in clinical units
row_med <- function(var, label, digits = 0L) {
    data.table(
        Variable = label,
        Train    = fmt_med_iqr(unscale_col(sets$Train[[var]], var), digits),
        Test     = fmt_med_iqr(unscale_col(sets$Test[[var]],  var), digits)
    )
}

# n (%) row for any condition, given as a function of the data
row_cnt <- function(expr_fun, label) {
    data.table(
        Variable = label,
        Train    = fmt_n_pct(sum(expr_fun(sets$Train), na.rm = TRUE), nrow(sets$Train)),
        Test     = fmt_n_pct(sum(expr_fun(sets$Test),  na.rm = TRUE), nrow(sets$Test))
    )
}

# One n (%) row per level of a categorical variable
row_levels <- function(var, label) {
    lv <- sort(unique(as.character(dtbin[[var]])))
    rbindlist(lapply(lv, function(l) row_cnt(function(d) as.character(d[[var]]) == l, paste0(label, ", ", l))))
}

# Share without a qualifying result
row_na <- function(var, label) {
    row_cnt(function(d) is.na(d[[var]]), paste0("No qualifying ", label, ", n (%)"))
}

lab_cols <- c(
    latest_hb         = "haemoglobin",
    latest_tromb      = "platelet count",
    latest_creatinine = "creatinine",
    latest_crp        = "CRP",
    latest_mcv        = "MCV",
    latest_inr        = "INR",
    latest_albumin    = "albumin",
    latest_wbc        = "leucocyte count",
    latest_ferritin   = "ferritin"
)

char_dt <- rbindlist(c(
    list(
        data.table(Variable = "Procedures", Train = format(nrow(sets$Train), big.mark = ","), Test = format(nrow(sets$Test), big.mark = ",")),
        data.table(Variable = "Patients",   Train = format(uniqueN(sets$Train$patient_idx), big.mark = ","), Test = format(uniqueN(sets$Test$patient_idx), big.mark = ",")),
        row_cnt(function(d) d$transfused, "Transfused, n (%)"),
        row_med("age", "Age, years, median (IQR)"),
        row_levels("sukupuoli_selite", "Sex"),
        row_levels("asa", "ASA class"),
        row_med("bmi", "BMI, median (IQR)"),
        row_med("bv_litres", "Estimated blood volume, L, median (IQR)", 1L),
        row_med("latest_hb", "Haemoglobin, g/L, median (IQR)"),
        row_med("latest_tromb", "Platelets, x10^9/L, median (IQR)"),
        row_med("latest_creatinine", "Creatinine, umol/L, median (IQR)"),
        row_med("comorbidity_index", "van Walraven index, median (IQR)"),
        # Episode counts are z-scored too, so unscaled and rounded back to whole numbers
        row_cnt(function(d) round(unscale_col(d$n_prev_bl_eps, "n_prev_bl_eps")) > 0, "Any prior recorded bleeding, n (%)"),
        row_cnt(function(d) d$anticoagulant == 1, "Anticoagulant, n (%)"),
        row_cnt(function(d) d$antiplatelet == 1, "Antiplatelet, n (%)"),
        row_cnt(function(d) d$iron_oral == 1 | d$iron_iv == 1, "Oral or IV iron, n (%)")
    ),
    lapply(names(lab_cols), function(v) row_na(v, lab_cols[[v]]))
))
save_tab(char_dt, "train_vs_test")

# ---------------------------------------------------------------------------
# 5. Missingness of every predictor
# ---------------------------------------------------------------------------
section("5. PREDICTOR MISSINGNESS (dtbin, after ce_01 handling)")

miss_dt <- rbindlist(lapply(
    c(fe_predictors, "toimenpide_tekija_id", "procedure_group"),
    function(v) {
        data.table(
            predictor  = v,
            n_na_all   = sum(is.na(dtbin[[v]])),
            n_na_train = sum(is.na(sets$Train[[v]])),
            n_na_test  = sum(is.na(sets$Test[[v]]))
        )
    }
))
miss_dt[, pct_na_all := round(100 * n_na_all / nrow(dtbin), 1)]
save_tab(miss_dt[n_na_all > 0][order(-n_na_all)], "predictor_missingness")
cat(sprintf("Predictors with no missing values after ce_01 handling: %d of %d\n", miss_dt[n_na_all == 0, .N], nrow(miss_dt)))

# Reported in Supplementary Methods 1
cat(sprintf("ASA 'Missing' category: %s\n", fmt_n_pct(sum(dtbin$asa == "Missing"), nrow(dtbin))))
cat(sprintf(
    "No recorded prior bleeding (NA set to 0 in ce_01): %s\n",
    fmt_n_pct(sum(round(unscale_col(dtbin$n_prev_bl_eps, "n_prev_bl_eps")) == 0), nrow(dtbin))
))

# ---------------------------------------------------------------------------
# 6. Hyperparameters (Table S10, Supplementary Methods 6)
# ---------------------------------------------------------------------------
section("6. FINAL HYPERPARAMETERS")

variants <- c(patient = "Patient-only", proc = "Patient + procedure", surgeon = "Patient + surgeon", both = "Full")

# Fitted model bundles from ce_02
bundles <- lapply(
    names(variants),
    function(v) readRDS(file.path(objpath, paste0("model_cb_", v, "_", model_tag, ".rds")))
)
names(bundles) <- names(variants)

# Tuned values from ce_02a, plus the tree count actually fitted
hp_dt <- rbindlist(lapply(
    names(variants),
    function(v) {
        tuned <- readRDS(file.path(objpath, paste0("tuned_params_cb_", v, "_", model_tag, ".rds")))
        data.table(
            model            = variants[[v]],
            learning_rate    = signif(tuned$params$learning_rate, 3),
            depth            = tuned$params$depth,
            l2_leaf_reg      = tuned$params$l2_leaf_reg,
            iterations       = tuned$nrounds,
            fitted_trees     = bundles[[v]]$model$tree_count,
            cv_auprc_mean    = round(tuned$cv_score, 4),
            cv_auprc_fold_sd = round(sd(tuned$fold_scores), 4)
        )
    }
))
save_tab(hp_dt, "hyperparameters")
cat(
    "Search: 30 random draws; learning_rate log-uniform 10^-2 to 10^-0.5; depth 4-10;",
    "l2_leaf_reg in {1, 3, 5, 10}; loss Logloss; early stopping after 50 rounds",
    "(max 2000); iterations = mean best iteration across 5 patient-grouped folds.\n"
)

# Every CatBoost setting of the full model, defaults included (the defaults table in Supplementary Methods 6)
# (wrapped, since the accessor differs between package versions)
cb_params <- tryCatch(catboost.get_model_params(bundles$both$model), error = function(e) NULL)
if (is.null(cb_params)) {
    cat("catboost.get_model_params() not available in this catboost version; skipped.\n")
} else {
    flat <- unlist(cb_params)
    keep <- grep(
        "one_hot_max_size|border_count|nan_mode|max_ctr_complexity|boosting_type|bootstrap_type|subsample|bagging_temperature|random_strength|leaf_estimation|grow_policy|ctr|loss_function|model_size_reg|rsm",
        names(flat),
        value = TRUE
    )
    print(data.table(setting = keep, value = unname(flat[keep])), nrows = 200L)
}

# ---------------------------------------------------------------------------
# 7. Calibration (Table S4)
# ---------------------------------------------------------------------------
section("7. CALIBRATION (test set, patient-cluster bootstrap)")

test_dt <- copy(sets$Test)
y_test  <- as.numeric(test_dt$transfused)

# Test-set probabilities of one model (pool built as in ce_03)
predict_bundle <- function(bundle) {
    cols     <- bundle$feature_names
    cat_cols <- bundle$cat_feature_names
    test_cb  <- copy(test_dt[, ..cols])
    test_cb[, (cat_cols) := lapply(.SD, as.factor), .SDcols = cat_cols]
    pool <- catboost.load_pool(data = as.data.frame(test_cb))
    catboost.predict(bundle$model, pool, prediction_type = "Probability")
}

# Columns are models, rows test procedures
pred_mat <- sapply(bundles, predict_bundle)
colnames(pred_mat) <- names(variants)
rm(bundles); gc()

auprc_of <- function(p, y) {
    pr.curve(scores.class0 = p, weights.class0 = y)$auc.integral
}

# Calibration intercept: logistic model with the predicted logit as offset (calibration-in-the-large)
# Slope: coefficient of the predicted logit (1 == perfect, below 1 == too extreme)
# Predictions are clipped away from 0 and 1 so the logit stays finite
cal_of <- function(p, y) {
    lp    <- qlogis(pmin(pmax(p, 1e-8), 1 - 1e-8))
    int   <- suppressWarnings(glm.fit(x = matrix(1, length(y)), y = y, offset = lp, family = binomial()))$coefficients[1L]
    slope <- suppressWarnings(glm.fit(x = cbind(1, lp), y = y, family = binomial()))$coefficients[2L]
    c(
        cal_intercept = unname(int),
        cal_slope     = unname(slope),
        oe_ratio      = mean(y) / mean(p),
        brier         = mean((p - y)^2)
    )
}

# Patient-cluster bootstrap over any subset of test rows (same scheme as ce_03)
cluster_boot <- function(rows, stat_fun, n_boot = N_BOOT, seed = 2028L) {
    by_pt <- split(rows, test_dt$patient_idx[rows])
    point <- stat_fun(rows)

    # Reproducible parallel streams, old generator restored afterwards
    old_rng <- RNGkind("L'Ecuyer-CMRG")[1]
    set.seed(seed)
    bm <- do.call(
        cbind,
        mclapply(
            seq_len(n_boot),
            function(b) stat_fun(unlist(by_pt[sample.int(length(by_pt), replace = TRUE)], use.names = FALSE)),
            mc.cores = n_cores
        )
    )
    RNGkind(old_rng)

    data.table(
        statistic = names(point),
        estimate  = unname(point),
        ci_lo     = apply(bm, 1, quantile, 0.025, na.rm = TRUE),
        ci_hi     = apply(bm, 1, quantile, 0.975, na.rm = TRUE)
    )
}

# All calibration measures for all four models, e.g. both_cal_slope
cal_stats <- function(idx) {
    y <- y_test[idx]
    unlist(lapply(
        names(variants),
        function(v) {
            s <- cal_of(pred_mat[idx, v], y)
            setNames(s, paste0(v, "_", names(s)))
        }
    ))
}
cal_ci <- cluster_boot(seq_len(nrow(test_dt)), cal_stats)

# Split "both_cal_slope" into model and statistic
cal_ci[, model := unname(variants[sub("_.*", "", statistic)])]
cal_ci[, statistic := sub("^[a-z]+_", "", statistic)]
save_tab(
    cal_ci[, .(model, statistic, estimate = round(estimate, 4), ci_lo = round(ci_lo, 4), ci_hi = round(ci_hi, 4))],
    "calibration"
)

# ---------------------------------------------------------------------------
# 8. Subgroups by sex and age (Table S5)
# ---------------------------------------------------------------------------
section("8. PERFORMANCE BY SEX AND AGE GROUP (full model; gain vs patient + procedure)")

test_dt[, age_raw := unscale_col(age, "age")]
test_dt[, age_group := fcase(age_raw < 50, "<50", age_raw < 70, "50-69", default = ">=70")]
test_dt[, sex := fifelse(sukupuoli_selite == "Nainen", "Female", "Male")]

# Full-model AUPRC, lift and calibration, and the gain over the patient + procedure model
sub_stats <- function(idx) {
    y      <- y_test[idx]
    a_both <- auprc_of(pred_mat[idx, "both"], y)
    a_proc <- auprc_of(pred_mat[idx, "proc"], y)
    c(
        auprc_full               = a_both,
        lift_full                = a_both / mean(y),
        delta_auprc_full_vs_proc = a_both - a_proc,
        cal_of(pred_mat[idx, "both"], y)
    )
}

subgroups <- list(
    list(var = "sex",       levels = c("Female", "Male")),
    list(var = "age_group", levels = c("<50", "50-69", ">=70"))
)

# Bootstrap within each subgroup separately
sub_dt <- rbindlist(lapply(
    subgroups,
    function(sg) {
        rbindlist(lapply(
            sg$levels,
            function(lv) {
                rows <- which(test_dt[[sg$var]] == lv)
                ci   <- cluster_boot(rows, sub_stats)
                ci[, `:=`(
                    subgroup   = sg$var,
                    level      = lv,
                    n          = length(rows),
                    events     = sum(y_test[rows]),
                    prevalence = round(mean(y_test[rows]), 4)
                )]
            }
        ))
    }
))
setcolorder(sub_dt, c("subgroup", "level", "n", "events", "prevalence", "statistic"))
save_tab(
    sub_dt[, .(subgroup, level, n, events, prevalence, statistic, estimate = round(estimate, 4), ci_lo = round(ci_lo, 4), ci_hi = round(ci_hi, 4))],
    "subgroups"
)
cat("Note: AUPRC depends on prevalence, so compare lift (AUPRC / prevalence) across subgroups.\n")

# ---------------------------------------------------------------------------
# 9. Ten most transfused procedure groups (Table S8)
# ---------------------------------------------------------------------------
section("9. PERFORMANCE WITHIN THE TEN MOST TRANSFUSED PROCEDURE GROUPS")

# Same ten groups as Figure 3
# Procedure group is constant within a group, so full vs patient + procedure == the within-group value of surgeon identity
top10_groups <- (
    dtbin
    [, .(n_transfused = sum(transfused)), by = procedure_group]
    [order(-n_transfused)]
    [1:10, procedure_group]
)

grp_stats <- function(idx) {
    y      <- y_test[idx]
    a_both <- auprc_of(pred_mat[idx, "both"], y)
    a_proc <- auprc_of(pred_mat[idx, "proc"], y)
    c(
        auprc_proc               = a_proc,
        auprc_full               = a_both,
        delta_auprc_full_vs_proc = a_both - a_proc,
        oe_ratio_full            = mean(y) / mean(pred_mat[idx, "both"])
    )
}

grp_dt <- rbindlist(lapply(
    top10_groups,
    function(g) {
        rows <- which(test_dt$procedure_group == g)

        # Too few test transfusions for a meaningful AUPRC (just in case)
        if (sum(y_test[rows]) < 10L) {
            return(data.table(
                procedure_group = g,
                n               = length(rows),
                events          = sum(y_test[rows]),
                statistic       = "too_few_events",
                estimate        = NA_real_,
                ci_lo           = NA_real_,
                ci_hi           = NA_real_
            ))
        }

        ci <- cluster_boot(rows, grp_stats, n_boot = N_BOOT_GROUP)
        ci[, `:=`(procedure_group = g, n = length(rows), events = sum(y_test[rows]))]
    }
))
setcolorder(grp_dt, c("procedure_group", "n", "events", "statistic"))
save_tab(
    grp_dt[, .(procedure_group, n, events, statistic, estimate = round(estimate, 3), ci_lo = round(ci_lo, 3), ci_hi = round(ci_hi, 3))],
    "procedure_group_performance"
)

# ---------------------------------------------------------------------------
# Save and session
# ---------------------------------------------------------------------------
saveRDS(
    list(
        flow             = flow_dt,
        hospitals        = site_dt,
        imputation       = imp_dt,
        split            = split_dt,
        folds            = fold_dt,
        train_vs_test    = char_dt,
        missingness      = miss_dt,
        hyperparameters  = hp_dt,
        calibration      = cal_ci,
        subgroups        = sub_dt,
        procedure_groups = grp_dt
    ),
    file.path(outdir, paste0("reporting_supplement_", model_tag, ".rds"))
)

section("SESSION")
cat(R.version.string, "\n")
cat(
    "catboost", as.character(packageVersion("catboost")),
    "| data.table", as.character(packageVersion("data.table")),
    "| PRROC", as.character(packageVersion("PRROC")), "\n"
)
cat("Outputs in:", outdir, "\n")

sink()
