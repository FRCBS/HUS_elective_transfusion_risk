# Where does the gain from surgeon identity come from, and do surgeons differ beyond the model?
#
# Follow-up to ce_05 for Supplementary Methods 10-11 and Tables S6-S8:
# gain by subset and group size, within-procedure discrimination, random-intercept
# models for surgeon variation, and out-of-sample observed/expected by surgeon category.
# Read-only: refits no CatBoost model, prints everything into a log.

library(data.table)
library(PRROC)
library(catboost)
library(parallel)
library(lme4)

model_tag <- "final"
n_cores   <- detectCores()

N_BOOT       <- 2000L  # patient-cluster bootstrap resamples
N_BOOT_SMALL <- 1000L  # resamples for the within-group statistics (slow)
MIN_EVENTS   <- 5L     # test transfusions a group needs for within-group AUPRC

# Paths are relative to the project root
resultspath <- "results"
objpath     <- file.path(resultspath, "model_objects")
outdir      <- file.path(resultspath, "reporting_supplement")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Everything printed also goes to this log
log_file <- file.path(outdir, paste0("ce_05b_output_", model_tag, ".txt"))
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

auprc_of <- function(p, y) {
    pr.curve(scores.class0 = p, weights.class0 = y)$auc.integral
}

# Keeps probabilities away from 0 and 1, so their logit stays finite
clip <- function(p) {
    pmin(pmax(p, 1e-8), 1 - 1e-8)
}

# ---------------------------------------------------------------------------
# Data and predictions
# Patient + procedure and full models
# ---------------------------------------------------------------------------
dtbin    <- readRDS(file.path(objpath, "dtbin.rds"))
test_dt  <- dtbin[set == "test"]
train_dt <- dtbin[set == "train"]

# Probabilities of one model for any rows (pool built as in ce_03)
predict_bundle <- function(bundle, dt) {
    cols     <- bundle$feature_names
    cat_cols <- bundle$cat_feature_names
    x <- copy(dt[, ..cols])
    x[, (cat_cols) := lapply(.SD, as.factor), .SDcols = cat_cols]
    pool <- catboost.load_pool(data = as.data.frame(x))
    catboost.predict(bundle$model, pool, prediction_type = "Probability")
}
b_proc <- readRDS(file.path(objpath, paste0("model_cb_proc_", model_tag, ".rds")))
b_full <- readRDS(file.path(objpath, paste0("model_cb_both_", model_tag, ".rds")))

# Training-set predictions are in-sample (optimistic!), used only in section 6
test_dt[,  `:=`(p_proc = predict_bundle(b_proc, test_dt),  p_full = predict_bundle(b_full, test_dt))]
train_dt[, `:=`(p_proc = predict_bundle(b_proc, train_dt), p_full = predict_bundle(b_full, train_dt))]
rm(b_proc, b_full); gc()

test_dt[,  y := as.numeric(transfused)]
train_dt[, y := as.numeric(transfused)]
test_dt[,  surgeon_id := as.character(toimenpide_tekija_id)]

# Should reproduce Table 2
cat(sprintf(
    "Check: pooled test AUPRC patient + procedure %.3f, full %.3f (Table 2: 0.430, 0.470)\n",
    auprc_of(test_dt$p_proc, test_dt$y),
    auprc_of(test_dt$p_full, test_dt$y)
))

# Same ten groups as Figure 3 and Table S8
top10_groups <- (
    dtbin
    [, .(n_tr = sum(transfused)), by = procedure_group]
    [order(-n_tr)]
    [1:10, procedure_group]
)

# Patient-cluster bootstrap over any subset of test rows (same scheme as ce_03 and ce_05)
cluster_boot <- function(rows, stat_fun, n_boot = N_BOOT, seed = 2029L) {
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

# AUPRC and lift of both models, and the gain from surgeon identity
pooled_stats <- function(idx) {
    y  <- test_dt$y[idx]
    ap <- auprc_of(test_dt$p_proc[idx], y)
    af <- auprc_of(test_dt$p_full[idx], y)
    c(auprc_proc = ap, auprc_full = af, delta = af - ap, lift_proc = ap / mean(y), lift_full = af / mean(y))
}

# Size of a subset of test rows
describe_rows <- function(rows) {
    data.table(
        n          = length(rows),
        events     = sum(test_dt$y[rows]),
        groups     = uniqueN(test_dt$procedure_group[rows]),
        prevalence = round(mean(test_dt$y[rows]), 4)
    )
}

# ---------------------------------------------------------------------------
# 1. Pooled gain inside vs outside the ten most transfused groups
# ---------------------------------------------------------------------------
section("1. POOLED GAIN INSIDE VS OUTSIDE THE TEN MOST TRANSFUSED PROCEDURE GROUPS")

sets1 <- list(
    all       = seq_len(nrow(test_dt)),
    top10     = which(test_dt$procedure_group %in% top10_groups),
    other_284 = which(!test_dt$procedure_group %in% top10_groups)
)
s1 <- rbindlist(lapply(
    names(sets1),
    function(nm) cbind(data.table(subset = nm), describe_rows(sets1[[nm]]), cluster_boot(sets1[[nm]], pooled_stats))
))
save_tab(
    s1[, .(subset, n, events, groups, prevalence, statistic, estimate = round(estimate, 4), ci_lo = round(ci_lo, 4), ci_hi = round(ci_hi, 4))],
    "gain_top10_vs_other"
)

# ---------------------------------------------------------------------------
# 2. Gain by procedure-group size (Supplementary Methods 10)
# ---------------------------------------------------------------------------
section("2. POOLED GAIN BY PROCEDURE-GROUP TRAINING VOLUME (tertiles of test procedures)")

# Group size = training procedures in the group (what CatBoost's procedure statistic was estimated from)
# Tertiles are cut over test procedures, so each holds about a third of the test set
grp_size <- train_dt[, .(train_n = .N), by = procedure_group]
test_dt[grp_size, train_n := i.train_n, on = "procedure_group"]
cuts <- quantile(test_dt$train_n, c(1/3, 2/3))
test_dt[, size_tertile := fcase(train_n <= cuts[1], "1 small", train_n <= cuts[2], "2 medium", default = "3 large")]
cat(sprintf("Tertile cut points (training procedures per group): %.0f and %.0f\n", cuts[1], cuts[2]))

s2 <- rbindlist(lapply(
    sort(unique(test_dt$size_tertile)),
    function(tt) {
        rows <- which(test_dt$size_tertile == tt)
        cbind(
            data.table(size_tertile = tt, train_n_range = paste(range(test_dt$train_n[rows]), collapse = "-")),
            describe_rows(rows),
            cluster_boot(rows, pooled_stats)
        )
    }
))
save_tab(
    s2[, .(size_tertile, train_n_range, n, events, groups, prevalence, statistic, estimate = round(estimate, 4), ci_lo = round(ci_lo, 4), ci_hi = round(ci_hi, 4))],
    "gain_by_group_size"
)

# ---------------------------------------------------------------------------
# 3. Within-procedure discrimination (Supplementary Methods 10)
# ---------------------------------------------------------------------------
section("3. WITHIN-PROCEDURE DISCRIMINATION ACROSS ALL GROUPS (event-weighted)")

# Pooled AUPRC is mostly about telling procedure types apart; this looks only within each group
# Groups need at least MIN_EVENTS transfusions and at least one non-transfused procedure
elig <- test_dt[, .(ev = sum(y), n = .N), by = procedure_group][ev >= MIN_EVENTS & ev < n, procedure_group]
cat(sprintf(
    "Groups with >= %d test transfusions: %d of %d; they hold %d of %d test transfusions\n",
    MIN_EVENTS,
    length(elig),
    uniqueN(test_dt$procedure_group),
    test_dt[procedure_group %in% elig, sum(y)],
    sum(test_dt$y)
))

# Within-group AUPRC and lift per group, averaged with weights by the group's transfusions
within_stats <- function(idx) {
    d <- test_dt[idx][procedure_group %in% elig, .(procedure_group, y, p_proc, p_full)]

    # A resample can leave a group with no events (or only events), which has no AUPRC
    r <- d[, if (sum(y) > 0 && sum(y) < .N) .(
        ev   = sum(y),
        prev = mean(y),
        ap   = auprc_of(p_proc, y),
        af   = auprc_of(p_full, y)
    ), by = procedure_group]

    w <- r$ev / sum(r$ev)
    c(
        wmean_auprc_proc           = sum(w * r$ap),
        wmean_auprc_full           = sum(w * r$af),
        wmean_delta                = sum(w * (r$af - r$ap)),
        wmean_lift_proc            = sum(w * r$ap / r$prev),
        wmean_lift_full            = sum(w * r$af / r$prev),
        share_groups_gain_positive = mean(r$af > r$ap)
    )
}

s3_rows <- list(
    all_eligible = which(test_dt$procedure_group %in% elig),
    top10        = which(test_dt$procedure_group %in% top10_groups),
    other        = which(test_dt$procedure_group %in% setdiff(elig, top10_groups))
)
s3 <- rbindlist(lapply(
    names(s3_rows),
    function(nm) cbind(data.table(subset = nm), cluster_boot(s3_rows[[nm]], within_stats, n_boot = N_BOOT_SMALL))
))
save_tab(
    s3[, .(subset, statistic, estimate = round(estimate, 4), ci_lo = round(ci_lo, 4), ci_hi = round(ci_hi, 4))],
    "within_procedure_gain"
)

# Point estimates for each eligible group
per_grp <- test_dt[procedure_group %in% elig, .(
    n          = .N,
    events     = sum(y),
    train_n    = train_n[1],
    auprc_proc = auprc_of(p_proc, y),
    auprc_full = auprc_of(p_full, y)
), by = procedure_group]
per_grp[, `:=`(delta = auprc_full - auprc_proc, top10 = procedure_group %in% top10_groups)]
setorder(per_grp, -events)
save_tab(
    per_grp[, .(procedure_group, top10, n, events, train_n, auprc_proc = round(auprc_proc, 3), auprc_full = round(auprc_full, 3), delta = round(delta, 3))],
    "within_procedure_by_group"
)

# ---------------------------------------------------------------------------
# 4. Surgeon variation beyond the model (Table S7)
# ---------------------------------------------------------------------------
section("4. SURGEON HETEROGENEITY BEYOND THE MODEL (random-intercept logistic, test set)")

# Test set only, so the predictions are out of sample
# The model's logit is an offset, so the random effects measure variation the model does NOT explain:
#   (1 | grp)       leftover miscalibration of whole procedure groups
#   (1 | surg)      a surgeon's general tendency, given their patients and procedures
#   (1 | grp:surg)  the same, allowed to differ between procedure groups
td <- test_dt[, .(
    y       = as.integer(y),
    lp_proc = qlogis(clip(p_proc)),
    lp_full = qlogis(clip(p_full)),
    grp     = procedure_group,
    surg    = surgeon_id
)]
ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

fit <- function(fml, data) {
    glmer(fml, data = data, family = binomial, control = ctrl)
}

# Median odds ratio: between two random surgeons for comparable patients
mor <- function(sd) {
    exp(sqrt(2) * sd * qnorm(0.75))
}

# SD, profile CI, MOR and likelihood-ratio test of one model against the group-only model
summarise_fit <- function(label, m_base, m, data_label, offset_label) {
    vc  <- as.data.frame(VarCorr(m))
    lrt <- anova(m_base, m)

    # A variance can't go below zero (tested on the boundary), so the chi-square p-value is halved
    p_boundary <- lrt[["Pr(>Chisq)"]][2] / 2

    # Profile intervals sometimes fail to converge (the "-" cells in Table S7)
    ci <- tryCatch(
        {
            cf <- confint(m, parm = "theta_", method = "profile", oldNames = FALSE, quiet = TRUE)
            data.table(term = rownames(cf), sd_lo = cf[, 1], sd_hi = cf[, 2])
        },
        error = function(e) NULL
    )

    out <- data.table(
        model  = label,
        data   = data_label,
        offset = offset_label,
        term   = paste0("sd_(Intercept)|", vc$grp),
        sd     = vc$sdcor
    )
    if (!is.null(ci)) out <- merge(out, ci, by = "term", all.x = TRUE)
    out[, `:=`(mor = mor(sd), lrt_chisq_vs_base = lrt$Chisq[2], p_boundary = p_boundary, aic = AIC(m), n = nrow(data))]
    out[]
}

# Three surgeon models for each offset, all compared with the group-only model
run_set <- function(d, data_label) {
    res <- list()
    for (off in c("lp_proc", "lp_full")) {
        f_base <- as.formula(sprintf("y ~ 1 + offset(%s) + (1 | grp)", off))
        f_s    <- as.formula(sprintf("y ~ 1 + offset(%s) + (1 | grp) + (1 | surg)", off))
        f_gs   <- as.formula(sprintf("y ~ 1 + offset(%s) + (1 | grp) + (1 | grp:surg)", off))
        f_both <- as.formula(sprintf("y ~ 1 + offset(%s) + (1 | grp) + (1 | surg) + (1 | grp:surg)", off))
        m_base <- fit(f_base, d)
        off_lab <- if (off == "lp_proc") "patient + procedure" else "full"
        res[[length(res) + 1L]] <- summarise_fit("grp + surgeon",               m_base, fit(f_s, d),    data_label, off_lab)
        res[[length(res) + 1L]] <- summarise_fit("grp + grp:surgeon",           m_base, fit(f_gs, d),   data_label, off_lab)
        res[[length(res) + 1L]] <- summarise_fit("grp + surgeon + grp:surgeon", m_base, fit(f_both, d), data_label, off_lab)
    }
    rbindlist(res, fill = TRUE)
}

s4 <- rbind(
    run_set(td, "all test procedures"),
    run_set(td[grp %in% top10_groups], "ten most transfused groups"),
    fill = TRUE
)
num_cols <- intersect(c("sd", "sd_lo", "sd_hi", "mor", "lrt_chisq_vs_base", "aic"), names(s4))
s4[, (num_cols) := lapply(.SD, round, 4), .SDcols = num_cols]
s4[, p_boundary := signif(p_boundary, 3)]
save_tab(s4, "surgeon_heterogeneity_glmm")
cat(
    "Reading: a non-zero surgeon SD with the patient + procedure offset = surgeon variation\n",
    "within procedures beyond patient features. With the full-model offset it should shrink\n",
    "if the model has captured those surgeon effects out of sample.\n",
    sep = ""
)

# ---------------------------------------------------------------------------
# 5. Surgeon effects vs test-set outcomes (Table S6, Supplementary Methods 10)
# ---------------------------------------------------------------------------
section("5. SURGEON SHAP (ce_03) VS BREADTH, CASE MIX AND TEST-SET OBSERVED/EXPECTED")

# Per-surgeon mean SHAP and category from ce_03
ss <- fread(file.path(objpath, paste0("surgeon_shap_", model_tag, ".csv")))
ss[, surgeon_id := as.character(surgeon_id)]

# How many procedure groups a surgeon works in, and the share of their main group
breadth <- dtbin[, .(
    n_groups        = uniqueN(procedure_group),
    top_group_share = max(tabulate(factor(procedure_group))) / .N
), by = .(surgeon_id = as.character(toimenpide_tekija_id))]

# Observed vs expected transfusions per surgeon in the test set (outcomes the models never saw)
surg_test <- test_dt[, .(
    n_test      = .N,
    obs         = sum(y),
    exp_proc    = sum(p_proc),
    exp_full    = sum(p_full),
    mean_p_proc = mean(p_proc)
), by = surgeon_id]

sv <- merge(ss[, .(surgeon_id, mean_shap, status, volume)], breadth, by = "surgeon_id")
sv <- merge(sv, surg_test, by = "surgeon_id")
cat(sprintf("Surgeons in SHAP file: %d; matched: %d\n", nrow(ss), nrow(sv)))

# Exact Poisson CI for O/E, treating the expected count as fixed
oe_ci <- function(obs, ex) {
    c(
        oe = obs / ex,
        lo = if (obs > 0) qchisq(0.025, 2 * obs) / 2 / ex else 0,
        hi = qchisq(0.975, 2 * (obs + 1)) / 2 / ex
    )
}

# O/E against the patient + procedure model and the full model, per surgeon category
by_status <- sv[, {
    a <- oe_ci(sum(obs), sum(exp_proc))
    b <- oe_ci(sum(obs), sum(exp_full))
    .(
        surgeons               = .N,
        test_procedures        = sum(n_test),
        transfused             = sum(obs),
        # as.numeric, since median() of integers is integer for odd n and double for even n (data.table rejects that!)
        median_volume          = as.numeric(median(volume)),
        median_n_groups        = as.numeric(median(n_groups)),
        median_top_group_share = round(median(top_group_share), 2),
        median_case_mix_risk   = round(median(mean_p_proc), 4),
        oe_vs_patient_proc     = round(a[["oe"]], 3),
        oe_pp_lo               = round(a[["lo"]], 3),
        oe_pp_hi               = round(a[["hi"]], 3),
        oe_vs_full             = round(b[["oe"]], 3),
        oe_full_lo             = round(b[["lo"]], 3),
        oe_full_hi             = round(b[["hi"]], 3)
    )
}, by = status]
save_tab(by_status, "surgeon_status_oe")
cat(
    "Reading: if the SHAP surgeon effect reflects within-procedure practice that generalises,\n",
    "'Above average' surgeons should have O/E > 1 against the patient + procedure model in\n",
    "the TEST set (outcomes the model never saw), 'Below average' < 1, and O/E vs the full\n",
    "model should be close to 1 for all.\n",
    sep = ""
)

# Same, only procedures in the ten most transfused groups
top10_oe <- (
    merge(
        test_dt[procedure_group %in% top10_groups, .(obs = sum(y), exp_proc = sum(p_proc)), by = surgeon_id],
        sv[, .(surgeon_id, status)],
        by = "surgeon_id"
    )
    [, {
        a <- oe_ci(sum(obs), sum(exp_proc))
        .(
            surgeons           = .N,
            transfused         = sum(obs),
            oe_vs_patient_proc = round(a[["oe"]], 3),
            lo                 = round(a[["lo"]], 3),
            hi                 = round(a[["hi"]], 3)
        )
    }, by = status]
)
cat("\nWithin the ten most transfused groups only:\n")
save_tab(top10_oe, "surgeon_status_oe_top10")

# Spearman correlations of mean SHAP with surgeon characteristics, bootstrap over surgeons
# log O/E gets +0.5 on observed (many surgeons have no test transfusions)
sv[, log_oe_proc := log((obs + 0.5) / exp_proc)]

spear <- function(d, v) {
    cor(d$mean_shap, d[[v]], method = "spearman")
}

vars5 <- c("n_groups", "top_group_share", "mean_p_proc", "log_oe_proc")

# Rows are variables, columns resamples
set.seed(2030)
boot5 <- replicate(N_BOOT, {
    d <- sv[sample.int(.N, replace = TRUE)]
    vapply(vars5, function(v) spear(d, v), numeric(1))
})

cor5 <- data.table(
    variable = vars5,
    rho      = round(vapply(vars5, function(v) spear(sv, v), numeric(1)), 3),
    ci_lo    = round(apply(boot5, 1, quantile, 0.025), 3),
    ci_hi    = round(apply(boot5, 1, quantile, 0.975), 3)
)
save_tab(cor5, "surgeon_shap_correlates")
cat("log_oe_proc: test-set observed vs patient + procedure expected per surgeon (out of sample).\n")

# ---------------------------------------------------------------------------
# 6. Ten most transfused groups: train vs test (Table S8, Supplementary Methods 11)
# ---------------------------------------------------------------------------
section("6. TEN MOST TRANSFUSED GROUPS: TRAIN (in-sample) VS TEST CALIBRATION AND AUPRC")

# Training values are in-sample (optimistic!), shown only to see whether the test-set O/E > 1
# comes from selecting the groups by their transfusion counts
grp_cal <- function(d, lab) {
    d[procedure_group %in% top10_groups, .(
        set         = lab,
        n           = .N,
        events      = sum(y),
        prevalence  = round(mean(y), 4),
        mean_p_full = round(mean(p_full), 4),
        oe_proc     = round(sum(y) / sum(p_proc), 3),
        oe_full     = round(sum(y) / sum(p_full), 3),
        auprc_proc  = round(auprc_of(p_proc, y), 3),
        auprc_full  = round(auprc_of(p_full, y), 3)
    ), by = procedure_group]
}
s6 <- rbind(grp_cal(train_dt, "train"), grp_cal(test_dt, "test"))
setorder(s6, procedure_group, -set)
save_tab(s6, "top10_train_vs_test")

# Caesarean section (MCA): predictor profile vs all other test procedures
cat("\nCaesarean section (MCA) vs all other test procedures:\n")
mca_prof <- test_dt[, .(
    n             = .N,
    prevalence    = round(mean(y), 4),
    pct_asa_1_2   = round(100 * mean(asa %in% c("1", "2")), 1),
    pct_no_hb     = round(100 * mean(is.na(latest_hb)), 1),
    median_p_full = signif(median(p_full), 3),
    p90_p_full    = signif(quantile(p_full, 0.9), 3)
), by = .(group = fifelse(procedure_group == "MCA", "MCA", "other"))]
save_tab(mca_prof, "mca_profile")

# ---------------------------------------------------------------------------
# Session
# ---------------------------------------------------------------------------
section("SESSION")
cat(R.version.string, "\n")
cat("lme4", as.character(packageVersion("lme4")), "\n")
cat("Outputs in:", outdir, "\n")
sink()
