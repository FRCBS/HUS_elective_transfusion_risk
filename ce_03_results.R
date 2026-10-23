# Evaluates the four CatBoost models on the held-out test set.
#
# Produces the discrimination and calibration tables, PR curves,
# SHAP values, permutation importance and
# the surgeon-level figures.

library(data.table)
library(ggplot2)
library(scales)
library(parallel)
library(PRROC)
library(catboost)
library(tools)

model_tag <- "final"
n_cores   <- detectCores()

N_BOOT       <- 2000L  # patient-cluster bootstrap resamples
N_PERM       <- 20L    # shuffles per feature in permutation importance
N_CALIB_BINS <- 10L    # bins in the reliability diagram
TOP_FRAC     <- 0.05   # share of highest-risk procedures flagged

MIN_TEST_N_SURGEON  <- 15L  # test procedures a surgeon needs for Figure 4
MIN_GROUP_N_SURGEON <- 15L  # procedures within a group a surgeon needs for Figure 3
HIGH_RATE_RATIO     <- 1.5  # Figure 3 highlights surgeons above this times the group average

# Paths are relative to the project root
objpath <- file.path("results", "model_objects")
figpath <- file.path("results", "figures")

# Display names, shared by every table and figure
model_names <- c(
    patient = "Patient-only",
    proc    = "Patient + procedure",
    surgeon = "Patient + surgeon",
    both    = "Full (patient + procedure + surgeon)"
)

# Figures as PNG for drafts and SVG for the journal
save_pub_figure <- function(name, width, height, plot = last_plot()) {
    stem <- file.path(figpath, paste0(name, "_", model_tag))
    ggsave(paste0(stem, ".png"), plot = plot, bg = "white", width = width, height = height, units = "px")
    ggsave(paste0(stem, ".svg"), plot = plot, bg = "white", width = width, height = height, units = "px")
}

# Tables go to objpath as CSV, plus RDS when asked
save_table <- function(dt, name, rds = FALSE) {
    stem <- file.path(objpath, paste0(name, "_", model_tag))
    fwrite(dt, paste0(stem, ".csv"))
    if (rds) saveRDS(dt, paste0(stem, ".rds"))
}

# ---------------------------------------------------------------------------
# Load data and test set
# ---------------------------------------------------------------------------
# dtbin is the modelling dataset from ce_01, train/test split in `set`
dtbin   <- readRDS(file.path(objpath, "dtbin.rds"))
test_dt <- dtbin[set == "test"]

# Outcome as 0/1 (for PRROC and the Brier score)
actual_outcomes <- as.numeric(test_dt$transfused)

# Share transfused (no skill AUPRC)
prevalence <- mean(actual_outcomes)

# ---------------------------------------------------------------------------
# Load models and predict on the test set
# ---------------------------------------------------------------------------
# Each saved bundle holds the model plus the feature names it was trained on
load_and_predict_cb <- function(bundle_path, test_dt) {
    bundle   <- readRDS(bundle_path)
    cols     <- bundle$feature_names
    cat_cols <- bundle$cat_feature_names

    # Same feature columns, in the same order, as in training
    test_cb <- copy(test_dt[, ..cols])

    # CatBoost wants categorical features as factor columns in a data.frame
    test_cb[, (cat_cols) := lapply(.SD, as.factor), .SDcols = cat_cols]
    pool <- catboost.load_pool(data = as.data.frame(test_cb))

    list(
        bundle    = bundle,
        test_data = test_cb,
        pool      = pool,
        pred      = catboost.predict(bundle$model, pool, prediction_type = "Probability")
    )
}

# Model files are named by variant, e.g. model_cb_both_final.rds
model_file <- function(variant) {
    file.path(objpath, paste0("model_cb_", variant, "_", model_tag, ".rds"))
}

cb_patient_out <- load_and_predict_cb(model_file("patient"), test_dt)
cb_proc_out    <- load_and_predict_cb(model_file("proc"),    test_dt)
cb_surgeon_out <- load_and_predict_cb(model_file("surgeon"), test_dt)
cb_both_out    <- load_and_predict_cb(model_file("both"),    test_dt)

# Predicted transfusion probabilities, one column per model
test_dt[, p_cb_patient := cb_patient_out$pred]
test_dt[, p_cb_proc    := cb_proc_out$pred]
test_dt[, p_cb_surgeon := cb_surgeon_out$pred]
test_dt[, p_cb_both    := cb_both_out$pred]

# The full model is reused for SHAP and permutation importance
cb_both_bundle <- cb_both_out$bundle
test_cb_both   <- cb_both_out$test_data
test_pool_both <- cb_both_out$pool

# ---------------------------------------------------------------------------
# Discrimination
# ---------------------------------------------------------------------------
# AUPRC is the primary metric (about 2% outcome)
# curve = TRUE keeps the curve points for Figure 1
pr_cb_patient <- pr.curve(scores.class0 = test_dt$p_cb_patient, weights.class0 = actual_outcomes, curve = TRUE)
pr_cb_proc    <- pr.curve(scores.class0 = test_dt$p_cb_proc,    weights.class0 = actual_outcomes, curve = TRUE)
pr_cb_surgeon <- pr.curve(scores.class0 = test_dt$p_cb_surgeon, weights.class0 = actual_outcomes, curve = TRUE)
pr_cb_both    <- pr.curve(scores.class0 = test_dt$p_cb_both,    weights.class0 = actual_outcomes, curve = TRUE)

# Brier score is the mean squared error of the predicted probabilities, lower is better
brier_cb_patient <- mean((test_dt$p_cb_patient - actual_outcomes)^2)
brier_cb_proc    <- mean((test_dt$p_cb_proc    - actual_outcomes)^2)
brier_cb_surgeon <- mean((test_dt$p_cb_surgeon - actual_outcomes)^2)
brier_cb_both    <- mean((test_dt$p_cb_both    - actual_outcomes)^2)

# ---------------------------------------------------------------------------
# Patient-cluster bootstrap CIs (Table 2, Results)
# ---------------------------------------------------------------------------
# AUPRC only, without the curve (this runs thousands of times!)
auprc_of <- function(p, y) {
    pr.curve(scores.class0 = p, weights.class0 = y)$auc.integral
}

# AUROC is reported in Table 2 only, for comparison with published models
auroc_of <- function(p, y) {
    roc.curve(scores.class0 = p, weights.class0 = y)$auc
}

# Flags the top 5% of predicted risk and counts the transfusions caught
top_frac_capture <- function(p, y, frac = TOP_FRAC) {
    n_flag  <- ceiling(frac * length(p))
    flagged <- order(p, decreasing = TRUE)[seq_len(n_flag)]
    tp      <- sum(y[flagged])
    c(
        recall      = tp / sum(y),           # share of all transfusions caught
        tp_per_1000 = 1000 * tp / n_flag     # transfused cases per 1,000 flagged
    )
}

# Predictions as a matrix with one column per model, for fast row subsetting
pred_mat <- as.matrix(test_dt[, .(p_cb_patient, p_cb_proc, p_cb_surgeon, p_cb_both)])
colnames(pred_mat) <- names(model_names)

# Every statistic for one set of row indices
# All models are scored on the same rows, so model differences are paired
boot_stats <- function(idx) {
    y        <- actual_outcomes[idx]
    a        <- apply(pred_mat[idx, ], 2, auprc_of, y = y)
    r        <- apply(pred_mat[idx, ], 2, auroc_of, y = y)
    cap_proc <- top_frac_capture(pred_mat[idx, "proc"], y)
    cap_both <- top_frac_capture(pred_mat[idx, "both"], y)
    c(
        setNames(a,           paste0("auprc_", names(a))),
        # Lift is AUPRC relative to the resample's own prevalence
        setNames(a / mean(y), paste0("lift_",  names(a))),
        setNames(r,           paste0("auroc_", names(r))),
        delta_proc_vs_patient    = a[["proc"]]    - a[["patient"]],
        delta_surgeon_vs_patient = a[["surgeon"]] - a[["patient"]],
        delta_both_vs_patient    = a[["both"]]    - a[["patient"]],
        # Primary comparison: what surgeon identity adds beyond procedure group
        delta_both_vs_proc       = a[["both"]]    - a[["proc"]],
        rel_gain_both_vs_proc    = a[["both"]] / a[["proc"]] - 1,
        top5_recall_proc         = cap_proc[["recall"]],
        top5_recall_both         = cap_both[["recall"]],
        top5_tp_per_1000_proc    = cap_proc[["tp_per_1000"]],
        top5_tp_per_1000_both    = cap_both[["tp_per_1000"]],
        top5_extra_tp_per_1000   = cap_both[["tp_per_1000"]] - cap_proc[["tp_per_1000"]]
    )
}

# Patients are resampled, not procedures, since one patient's procedures aren't independent
# Each entry holds the test-set row numbers of one patient
rows_by_patient <- split(seq_len(nrow(test_dt)), test_dt$patient_idx)
n_test_patients <- length(rows_by_patient)

# Point estimates use the full test set
point_est <- boot_stats(seq_len(nrow(test_dt)))

# Parallel workers need the L'Ecuyer generator for reproducible random streams
# The previous generator is restored afterwards, so later seeds behave as usual
old_rng <- RNGkind("L'Ecuyer-CMRG")[1]
set.seed(2026)

# Rows are statistics, columns bootstrap resamples
boot_mat <- do.call(
    cbind,
    mclapply(
        seq_len(N_BOOT),
        function(b) {
            # Draw patients with replacement and bring along all their procedures
            idx <- unlist(rows_by_patient[sample.int(n_test_patients, replace = TRUE)], use.names = FALSE)
            boot_stats(idx)
        },
        mc.cores = n_cores
    )
)
RNGkind(old_rng)

# Percentile 95% CIs from the bootstrap distribution
boot_ci <- data.table(
    statistic = names(point_est),
    estimate  = unname(point_est),
    ci_lo     = apply(boot_mat, 1, quantile, 0.025),
    ci_hi     = apply(boot_mat, 1, quantile, 0.975)
)
print(boot_ci)
save_table(boot_ci, "bootstrap_ci", rds = TRUE)

# Formats a statistic's CI as "lo to hi"
fmt_ci <- function(stat, digits = 3) {
    r <- boot_ci[statistic == stat]
    sprintf(paste0("%.", digits, "f to %.", digits, "f"), r$ci_lo, r$ci_hi)
}

# ---------------------------------------------------------------------------
# Table 2
# ---------------------------------------------------------------------------
# Full test-set AUPRC per model, in model_names order
auprc_vec <- c(
    pr_cb_patient$auc.integral,
    pr_cb_proc$auc.integral,
    pr_cb_surgeon$auc.integral,
    pr_cb_both$auc.integral
)

# Full test-set AUROC per model, from the bootstrap point estimates
auroc_vec <- point_est[paste0("auroc_", names(model_names))]

# First row is the no-skill baseline, whose AUPRC equals prevalence and AUROC is 0.5
# Deltas are against the patient-only model, so its own row is left empty
discrim_dt <- data.table(
    Model           = c("No-skill baseline", unname(model_names)),
    `Grouping vars` = c("N/A", "None", "procedure", "surgeon", "surgeon + procedure"),
    AUPRC           = round(c(prevalence, auprc_vec), 3),
    AUPRC_95CI      = c(NA, vapply(paste0("auprc_", names(model_names)), fmt_ci, character(1))),
    Lift            = c(1, round(auprc_vec / prevalence, 1)),
    AUPRC_delta_vs_patient = c(NA, NA, round(auprc_vec[2:4] - auprc_vec[1], 3)),
    AUPRC_delta_95CI = c(
        NA,
        NA,
        vapply(
            c("delta_proc_vs_patient", "delta_surgeon_vs_patient", "delta_both_vs_patient"),
            fmt_ci,
            character(1)
        )
    ),
    Brier = c(NA, round(c(brier_cb_patient, brier_cb_proc, brier_cb_surgeon, brier_cb_both), 4)),
    AUROC = round(c(0.5, unname(auroc_vec)), 3),
    AUROC_95CI = c(NA, vapply(paste0("auroc_", names(model_names)), fmt_ci, character(1)))
)

print(discrim_dt)
save_table(discrim_dt, "table2_discrimination", rds = TRUE)

# ---------------------------------------------------------------------------
# Calibration of the full model (Figure S3)
# ---------------------------------------------------------------------------
# Reliability diagram: mean predicted vs observed rate within bins of predicted risk
calib_dt <- data.table(p_pred = test_dt$p_cb_both, y = actual_outcomes)

# Equal-sized bins by rank, since quantile breaks collapse on the many near-zero predictions
calib_dt[, bin := pmin(
    N_CALIB_BINS,
    ceiling(frank(p_pred, ties.method = "first") / .N * N_CALIB_BINS)
)]

# Exact binomial 95% CI for each bin's observed rate
calib_summary <- (
    calib_dt
    [, {
        ci <- binom.test(sum(y), .N)$conf.int
        .(
            mean_pred = mean(p_pred),
            observed  = mean(y),
            ci_lo     = ci[1],
            ci_hi     = ci[2],
            n         = .N
        )
    }, by = bin]
    [order(bin)]
)

save_table(calib_summary, "calibration_data_both", rds = TRUE)

# Square axes that fit the widest CI, so the diagonal sits at 45 degrees
calib_axis_max <- max(calib_summary$ci_hi, calib_summary$mean_pred) * 1.05

# Perfect calib == points exactly on the line
ggplot(calib_summary, aes(x = mean_pred, y = observed)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0, color = "#2a78d6") +
    geom_line(color = "#2a78d6") +
    geom_point(size = 2.5, color = "#2a78d6") +
    coord_equal(xlim = c(0, calib_axis_max), ylim = c(0, calib_axis_max)) +
    labs(
        x = "Mean predicted probability (decile bin)",
        y = "Observed event rate (95% CI)"
    ) +
    theme_minimal() +
    theme(text = element_text(size = 18))

save_pub_figure("calibration_both", width = 2500, height = 2500)

# ---------------------------------------------------------------------------
# PR curves (Figure 1)
# ---------------------------------------------------------------------------
# Legend labels carry each model's AUPRC, e.g. "Patient-only (0.255)"
pr_model_levels <- sprintf("%s (%.3f)", model_names, auprc_vec)
pr_model_colors <- setNames(
    c("#2a78d6", "#eb6834", "#1baf7a", "#eda100"),
    pr_model_levels
)

# Columns 1 and 2 of each PRROC curve are recall and precision
pr_curve_data <- rbind(
    data.frame(recall = pr_cb_patient$curve[, 1], precision = pr_cb_patient$curve[, 2], model = pr_model_levels[1]),
    data.frame(recall = pr_cb_proc$curve[, 1],    precision = pr_cb_proc$curve[, 2],    model = pr_model_levels[2]),
    data.frame(recall = pr_cb_surgeon$curve[, 1], precision = pr_cb_surgeon$curve[, 2], model = pr_model_levels[3]),
    data.frame(recall = pr_cb_both$curve[, 1],    precision = pr_cb_both$curve[, 2],    model = pr_model_levels[4])
)

# Factor levels keep the legend in model order instead of alphabetical
pr_curve_data$model <- factor(pr_curve_data$model, levels = pr_model_levels)

# Dashed line is the no-skill baseline, precision equal to prevalence
ggplot(pr_curve_data, aes(x = recall, y = precision, color = model)) +
    geom_line(linewidth = 0.8) +
    geom_hline(yintercept = prevalence, linetype = "dashed", color = "gray40") +
    annotate(
        "text",
        x     = 0.2,
        y     = prevalence + 0.03,
        label = sprintf("No-skill baseline (%.1f%%)", prevalence * 100),
        color = "gray40",
        size  = 6
    ) +
    scale_color_manual(values = pr_model_colors) +
    labs(
        x     = "Recall",
        y     = "Precision",
        color = "Model (AUPRC)"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom", text = element_text(size = 18)) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))

save_pub_figure("auprc", width = 3500, height = 2000)

# ---------------------------------------------------------------------------
# SHAP values of the full model
# ---------------------------------------------------------------------------
# SHAP splits each prediction into additive per-feature contributions on the log-odds scale
# Rows are test procedures, columns features
shap_mat <- catboost.get_feature_importance(
    cb_both_bundle$model,
    test_pool_both,
    type = "ShapValues"
)

# The last column is the bias term (the average prediction), not a feature
shap_mat <- shap_mat[, -ncol(shap_mat)]
colnames(shap_mat) <- cb_both_bundle$feature_names

# ---------------------------------------------------------------------------
# Permutation AUPRC importance
# ---------------------------------------------------------------------------
# A feature's importance is how much AUPRC drops when its values are shuffled
# Shuffling breaks the feature's link to the outcome while keeping its distribution
baseline_auprc <- pr_cb_both$auc.integral
perm_features  <- cb_both_bundle$feature_names

# Full-model AUPRC on a test set where one feature has been shuffled
auprc_with <- function(permuted_dt) {
    pool   <- catboost.load_pool(data = as.data.frame(permuted_dt))
    p_perm <- catboost.predict(cb_both_bundle$model, pool, prediction_type = "Probability")
    auprc_of(p_perm, actual_outcomes)
}

# Each feature is shuffled N_PERM times, since a single shuffle is noisy
# Rows are repeats, columns features, values AUPRC decreases
set.seed(42)
perm_drops <- vapply(
    perm_features,
    function(feat) {
        vapply(
            seq_len(N_PERM),
            function(r) {
                # Copy, so the shuffle never touches the original test data
                tmp <- copy(test_cb_both)
                set(tmp, j = feat, value = sample(tmp[[feat]]))
                baseline_auprc - auprc_with(tmp)
            },
            numeric(1)
        )
    },
    numeric(N_PERM)
)

# Mean drop per feature with the spread across repeats, ranked largest first
perm_dt <- data.table(
    feature    = perm_features,
    auprc_drop = colMeans(perm_drops),
    sd_drop    = apply(perm_drops, 2, sd),
    drop_lo    = apply(perm_drops, 2, quantile, 0.025),
    drop_hi    = apply(perm_drops, 2, quantile, 0.975)
)
setorder(perm_dt, -auprc_drop)
perm_dt[, rank := seq_len(.N)]

# How stable is surgeon identity's rank across repeats?
# Rank of every feature within each repeat, 1 being the largest drop
rank_by_rep         <- t(apply(perm_drops, 1, function(d) rank(-d, ties.method = "min")))
surg_rank_by_rep    <- rank_by_rep[, "toimenpide_tekija_id"]
surg_mean_rank      <- perm_dt[feature == "toimenpide_tekija_id", rank]
n_reps_at_surg_rank <- sum(surg_rank_by_rep == surg_mean_rank)
message(
    sprintf(
        "Surgeon identity: rank %d by mean drop; same rank in %d/%d repeats (range %d-%d)",
        surg_mean_rank,
        n_reps_at_surg_rank,
        N_PERM,
        min(surg_rank_by_rep),
        max(surg_rank_by_rep)
    )
)

# ---------------------------------------------------------------------------
# Surgeon permutation within procedure group (Results)
# ---------------------------------------------------------------------------
# A free shuffle creates surgeon-procedure pairs that never occur, which inflates the drop
# Shuffling only within procedure groups keeps every pair realistic
proc_vec <- test_cb_both$procedure_group

# Row order with surgeon IDs shuffled among rows of the same procedure group
within_group_perm <- function() {
    idx <- seq_along(proc_vec)
    # Groups with a single test procedure have nothing to shuffle and stay as they are
    as.integer(ave(idx, proc_vec, FUN = function(i) if (length(i) > 1L) sample(i) else i))
}

set.seed(43)
cond_surg_drops <- vapply(
    seq_len(N_PERM),
    function(r) {
        tmp <- copy(test_cb_both)
        set(tmp, j = "toimenpide_tekija_id", value = tmp$toimenpide_tekija_id[within_group_perm()])
        baseline_auprc - auprc_with(tmp)
    },
    numeric(1)
)

perm_cond_surg <- data.table(
    feature    = "toimenpide_tekija_id (within procedure group)",
    auprc_drop = mean(cond_surg_drops),
    sd_drop    = sd(cond_surg_drops),
    drop_lo    = quantile(cond_surg_drops, 0.025),
    drop_hi    = quantile(cond_surg_drops, 0.975)
)
print(perm_cond_surg)

save_table(perm_dt, "perm_auprc")
save_table(perm_cond_surg, "perm_auprc_surgeon_within_procedure")
save_table(
    data.table(repeat_id = seq_len(N_PERM), surgeon_rank = surg_rank_by_rep),
    "perm_surgeon_rank_by_repeat"
)

# ---------------------------------------------------------------------------
# Permutation importance figure (Figure 2)
# ---------------------------------------------------------------------------
# Readable axis labels for the feature column names
# Elixhauser codes are the category names of the comorbidity package
feature_labels <- c(
    chf                                = "Congestive heart failure",
    carit                              = "Cardiac arrhythmia",
    valv                               = "Valvular disease",
    pcd                                = "Pulmonary circulation disorder",
    pvd                                = "Peripheral vascular disorder",
    hypunc                             = "Hypertension, uncomplicated",
    hypc                               = "Hypertension, complicated",
    para                               = "Paralysis",
    ond                                = "Other neurological disorder",
    cpd                                = "Chronic pulmonary disease",
    diabunc                            = "Diabetes, uncomplicated",
    diabc                              = "Diabetes, complicated",
    hypothy                            = "Hypothyroidism",
    rf                                 = "Renal failure",
    ld                                 = "Liver disease",
    pud                                = "Peptic ulcer disease",
    aids                               = "HIV/AIDS",
    lymph                              = "Lymphoma",
    metacanc                           = "Metastatic cancer",
    solidtum                           = "Solid tumor (non-metastatic)",
    rheumd                             = "Rheumatoid/collagen vascular disease",
    coag                               = "Coagulopathy",
    obes                               = "Obesity",
    wloss                              = "Weight loss",
    fed                                = "Fluid and electrolyte disorder",
    blane                              = "Blood loss anaemia",
    dane                               = "Deficiency anaemia",
    alcohol                            = "Alcohol abuse",
    drug                               = "Drug abuse",
    psycho                             = "Psychoses",
    depre                              = "Depression",
    age                                = "Age",
    sukupuoli_selite                   = "Sex",
    bmi                                = "BMI",
    bv_litres                          = "Blood volume (L)",
    latest_hb                          = "Haemoglobin",
    latest_tromb                       = "Platelet count",
    latest_creatinine                  = "Creatinine",
    latest_crp                         = "C-reactive protein",
    latest_mcv                         = "Mean corpuscular volume",
    latest_inr                         = "INR",
    latest_albumin                     = "Albumin",
    latest_wbc                         = "Leucocyte count",
    latest_ferritin                    = "Ferritin",
    latest_hb_days_since_drawn         = "Days since haemoglobin result",
    latest_tromb_days_since_drawn      = "Days since platelet count result",
    latest_creatinine_days_since_drawn = "Days since creatinine result",
    latest_crp_days_since_drawn        = "Days since CRP result",
    latest_mcv_days_since_drawn        = "Days since MCV result",
    latest_inr_days_since_drawn        = "Days since INR result",
    latest_albumin_days_since_drawn    = "Days since albumin result",
    latest_wbc_days_since_drawn        = "Days since WBC result",
    latest_ferritin_days_since_drawn   = "Days since ferritin result",
    asa                                = "ASA class",
    total_prev_bl                      = "Prior bleeding burden",
    n_prev_bl_eps                      = "Prior bleeding episodes",
    comorbidity_index                  = "Comorbidity index",
    anticoagulant                      = "Anticoagulant",
    antiplatelet                       = "Antiplatelet",
    nsaid                              = "NSAID",
    corticosteroid                     = "Corticosteroid",
    coag_factor                        = "Coagulation factor concentrate",
    hemostatic                         = "Systemic haemostatic",
    iron_oral                          = "Oral iron",
    iron_iv                            = "IV iron",
    esa                                = "Erythropoiesis-stimulating agent",
    cytostatic                         = "Cytostatic",
    immunosuppressant                  = "Immunosuppressant",
    dapt                               = "Dual antiplatelet therapy",
    n_med_records                      = "Active prescription count",
    toimenpide_tekija_id               = "Surgeon identity",
    procedure_group                    = "Procedure group"
)

# Any feature not listed above gets a title-cased fallback
label_feature <- function(x) {
    known    <- feature_labels[x]
    fallback <- toTitleCase(gsub("_", " ", x))
    ifelse(is.na(known), fallback, known)
}

# perm_dt is sorted by drop, so the first 20 rows are the 20 most important features
top20_perm <- perm_dt[1:20]
top20_perm[, feature_label := label_feature(feature)]

# Reversed levels put the largest bar at the top of the plot
top20_perm[, feature_label := factor(feature_label, levels = rev(feature_label))]

# Error bars span the 2.5th to 97.5th percentile of the drop across repeats
ggplot(top20_perm, aes(x = auprc_drop, y = feature_label)) +
    geom_col(fill = "#2a78d6", alpha = 0.85) +
    geom_errorbarh(aes(xmin = drop_lo, xmax = drop_hi), height = 0.3, color = "gray30") +
    labs(
        x = sprintf("Mean AUPRC decrease when feature is shuffled (%d permutations, 95%% range)", N_PERM),
        y = NULL
    ) +
    theme_minimal()

save_pub_figure("perm_auprc", width = 2200, height = 2000)

# ---------------------------------------------------------------------------
# Surgeon effects from SHAP
# Also read by ce_04 and ce_05b
# ---------------------------------------------------------------------------
# A surgeon's effect is the mean SHAP contribution of surgeon identity over their test procedures
# Positive means the model raised predicted risk for that surgeon's cases
surgeon_shap <- (
    data.table(
        surgeon_id   = as.character(test_dt$toimenpide_tekija_id),
        shap_surgeon = shap_mat[, "toimenpide_tekija_id"]
    )
    [, .(
        mean_shap = mean(shap_surgeon),
        sd_shap   = sd(shap_surgeon),
        n         = .N
    ), by = surgeon_id]
    # Too few test procedures give an unstable mean
    [n >= MIN_TEST_N_SURGEON]
)

# ±1.96 SE across the surgeon's procedures
# Shows how consistently the model attributes risk 
# (not uncertainty in the effect itself!!)
surgeon_shap[, `:=`(
    ci_lo = mean_shap - 1.96 * sd_shap / sqrt(n),
    ci_hi = mean_shap + 1.96 * sd_shap / sqrt(n)
)]

# Above or below average means more than 1 SD of the surgeon effects from zero
# Not a statistical outlier test
shap_sd <- sd(surgeon_shap$mean_shap)
surgeon_shap[, status := fcase(
    mean_shap >  shap_sd, "Above average",
    mean_shap < -shap_sd, "Below average",
    default = "Within ±1 SD"
)]

# Surgeon title record likely includes just the newest (guess by data lake personnel)
# e.g. resident to specialist transition -> specialist across all procedures
# But just in case, each surgeon gets their most frequent title, ties broken alphabetically
surgeon_titles <- (
    dtbin
    [, .N, by = .(surgeon_id = as.character(toimenpide_tekija_id), toimenpide_tekija_ammatti_selite)]
    [order(surgeon_id, -N, toimenpide_tekija_ammatti_selite)]
    [, .SD[1L], by = surgeon_id]
    [, N := NULL]
)
surgeon_shap <- merge(surgeon_shap, surgeon_titles, by = "surgeon_id", all.x = TRUE)

# Figure 4 main panel shows only the above- and below-average surgeons
# Ranked by effect, so rank 1 is the most negative
outliers <- surgeon_shap[status != "Within ±1 SD"]
setorder(outliers, mean_shap)
outliers[, surgeon_rank := .I]

# All eligible surgeons, ranked the same way, for the Figure 4 inset
all_surgeons <- copy(surgeon_shap)
setorder(all_surgeons, mean_shap)
all_surgeons[, surgeon_rank := .I]

# ---------------------------------------------------------------------------
# Surgeon effects figure (Figure 4)
# ---------------------------------------------------------------------------
# coord_flip puts surgeons on the vertical axis, so x and y are swapped on screen
p_surgeon_effects <- ggplot(outliers, aes(x = surgeon_rank, y = mean_shap, color = status)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0, color = "black", linewidth = 0.3) +
    geom_point(size = 2) +
    coord_flip() +
    scale_color_manual(
        values = c(
            "Above average" = "#e34948",
            "Below average" = "#2a78d6"
        )
    ) +
    labs(
        y     = "Mean SHAP (log-odds)",
        x     = "Surgeon rank (by mean SHAP)",
        color = "Direction"
    ) +
    theme_classic()

# Inset: all eligible surgeons in miniature, those within 1 SD in grey
p_inset <- ggplot(all_surgeons, aes(x = surgeon_rank, y = mean_shap, color = status)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.2) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0, color = "black", linewidth = 0.15, alpha = 0.5) +
    geom_point(size = 0.4) +
    coord_flip() +
    scale_color_manual(
        values = c(
            "Above average" = "#e34948",
            "Below average" = "#2a78d6",
            "Within ±1 SD"  = "gray70"
        )
    ) +
    theme_void() +
    theme(
        legend.position  = "none",
        panel.background = element_rect(fill = "white", color = "gray50", linewidth = 0.3)
    )

# The low-rank, high-SHAP corner is empty
# The inset goes there, sized as a share of each axis range
inset_frac <- 0.32
rank_range <- range(outliers$surgeon_rank)
shap_range <- range(c(outliers$ci_lo, outliers$ci_hi))

# annotation_custom uses the axes before coord_flip: x is rank, y is SHAP
p_surgeon_effects <- p_surgeon_effects +
    annotation_custom(
        grob = ggplotGrob(p_inset),
        xmin = rank_range[1],
        xmax = rank_range[1] + inset_frac * diff(rank_range),
        ymin = shap_range[2] - inset_frac * diff(shap_range),
        ymax = shap_range[2]
    )

save_pub_figure("surgeon_effects", width = 2000, height = 3000, plot = p_surgeon_effects)

# ---------------------------------------------------------------------------
# Surgeon transfusion rates within procedure groups (Figure 3)
# ---------------------------------------------------------------------------
# Raw observed rates from the whole cohort, no model involved
# The ten groups with the most transfused cases; ranking by volume would favour low-risk groups
top_groups <- (
    dtbin
    [, .(n_transfused = sum(transfused == 1)), by = procedure_group]
    [order(-n_transfused)]
    [1:10, procedure_group]
)

# One row per surgeon and group, for surgeons with enough cases in that group
plot_dt <- (
    dtbin
    [procedure_group %in% top_groups, .(
        rate = mean(transfused == 1),
        vol  = .N
    ), by = .(toimenpide_tekija_id, procedure_group)]
    [vol >= MIN_GROUP_N_SURGEON]
)

# x-axis follows the transfused-case ranking instead of alphabetical order
plot_dt[, procedure_group := factor(procedure_group, levels = top_groups)]

# Group average is the unweighted mean of surgeon rates
group_means <- plot_dt[, .(avg_rate = mean(rate)), by = procedure_group]
plot_dt     <- merge(plot_dt, group_means, by = "procedure_group")

# The surgeons shown, for the figure legend's count in ce_04
save_table(plot_dt, "figure3_surgeon_rates")

# Boxplots show the spread across surgeons, points are individual surgeons sized by volume
ggplot(plot_dt, aes(x = procedure_group, y = rate)) +
    geom_boxplot(outlier.shape = NA, fill = "gray95", color = "gray60") +
    geom_jitter(aes(size = vol, color = rate > (avg_rate * HIGH_RATE_RATIO)), width = 0.2, alpha = 0.7) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_color_manual(
        values = c("FALSE" = "#2a78d6", "TRUE" = "#e34948"),
        labels = sprintf(c("Within %gx group average", "Above %gx group average"), HIGH_RATE_RATIO)
    ) +
    guides(size = guide_legend(nrow = 2)) +
    labs(
        y     = "Transfusion rate",
        x     = "Procedure group (NCSP, 3-character code)",
        size  = "Case volume",
        color = "Surgeon rate"
    ) +
    theme_minimal() +
    theme(
        legend.position = "bottom",
        text            = element_text(size = 18),
        legend.text     = element_text(size = 12),
        legend.title    = element_text(size = 14)
    )

save_pub_figure("surgeon_procedure_transfusion_rates", width = 3500, height = 2000)

# ---------------------------------------------------------------------------
# Surgeon seniority (Supplementary Methods 10)
# ---------------------------------------------------------------------------
# Finnish professional titles grouped as specialist or consultant (senior) versus not
# Confounded with case allocation, since ce_01 picks the primary surgeon by seniority
title_map <- c(
    "ylilääkäri"                = "Senior",
    "osastonylilääkäri"         = "Senior",
    "apulaisylilääkäri"         = "Senior",
    "apulaisylihammaslääkäri"   = "Senior",
    "erikoislääkäri"            = "Senior",
    "erikoishammaslääkäri"      = "Senior",
    "osastonlääkäri"            = "Not Senior",
    "osastonhammaslääkäri"      = "Not Senior",
    "erikoistuva lääkäri"       = "Not Senior",
    "erikoistuva hammaslääkäri" = "Not Senior",
    "lääkäri"                   = "Not Senior",
    "päivystävä lääkäri"        = "Not Senior",
    "amanuenssi"                = "Not Senior",
    "lääketieteen kandidaatti"  = "Not Senior",
    "lääketieteen opiskelija"   = "Not Senior",
    "kätilö"                    = "Not Senior",
    "tutkija"                   = "Not Senior",
    "tutkimusassistentti"       = "Not Senior"
)

# Missing or unlisted titles become "Unknown"
surgeon_shap[, experience := title_map[toimenpide_tekija_ammatti_selite]]
surgeon_shap[is.na(experience), experience := "Unknown"]

# Median (IQR) surgeon effect per seniority group
seniority_summary <- (
    surgeon_shap
    [, .(
        n_surgeons  = .N,
        median_shap = median(mean_shap),
        q25_shap    = quantile(mean_shap, 0.25),
        q75_shap    = quantile(mean_shap, 0.75)
    ), by = experience]
    [order(experience)]
)
print(seniority_summary)
save_table(seniority_summary, "surgeon_seniority_shap")

# ---------------------------------------------------------------------------
# Surgeon volume (Figure S4)
# ---------------------------------------------------------------------------
# Volume counts all of a surgeon's cohort procedures, train and test
# It measures activity in the study period, not career experience
surgeon_volume <- dtbin[, .(volume = .N), by = .(surgeon_id = as.character(toimenpide_tekija_id))]
surgeon_shap   <- merge(surgeon_shap, surgeon_volume, by = "surgeon_id", all.x = TRUE)

# Rank correlation, since volume is heavily skewed
spearman_rho <- function(d) {
    cor(d$volume, d$mean_shap, method = "spearman")
}

rho_est <- spearman_rho(surgeon_shap)

# 95% CI by resampling surgeons with replacement
set.seed(2027)
rho_boot <- vapply(
    seq_len(N_BOOT),
    function(b) spearman_rho(surgeon_shap[sample.int(.N, replace = TRUE)]),
    numeric(1)
)

volume_summary <- data.table(
    n_surgeons = nrow(surgeon_shap),
    rho        = rho_est,
    rho_ci_lo  = quantile(rho_boot, 0.025),
    rho_ci_hi  = quantile(rho_boot, 0.975)
)
print(volume_summary)
save_table(volume_summary, "surgeon_volume_shap")

# The steps in the smoother come from CatBoost binning surgeon frequency into at most 15 levels
ggplot(surgeon_shap, aes(x = volume, y = mean_shap)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    geom_point(alpha = 0.4, size = 1.2, color = "#2a78d6") +
    geom_smooth(method = "loess", formula = y ~ x, se = TRUE, color = "#e34948") +
    scale_x_log10(labels = comma) +
    labs(
        x = "Surgeon case volume in study period (log scale)",
        y = "Mean surgeon SHAP contribution (log-odds)"
    ) +
    theme_minimal() +
    theme(text = element_text(size = 14))

save_pub_figure("surgeon_volume_shap", width = 2200, height = 1800)

# Per-surgeon effects, titles and volumes, read by ce_04 and ce_05b
save_table(surgeon_shap, "surgeon_shap")

message("Done.")
