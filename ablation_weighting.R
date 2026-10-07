# Standalone ablation: does inverse-frequency class weighting actually help
# the ranking performance of the full model? It can bias the predictions.
#
# This script tunes the full model twice on the same random grid
# (weighted vs unweighted class_weights), fits their best params
# on a train/holdout split (not the real frozen split),
# then compares Brier score, mean P/O rate, and an external PRAUC.
#
# NB!!
# CatBoost's built-in "PRAUC" metric defaults to use_weights = TRUE.
# Early stopping here therefore uses "PRAUC:use_weights=false".
#
# Not part of the regular pipeline, nothing here can alter reproducibility.

library(data.table)
library(catboost)
library(parallel)
library(PRROC)

model_tag <- "final"
n_cores   <- detectCores()

N_RANDOM   <- 30L    # random hyperparameter draws per arm
MAX_ROUNDS <- 2000L  # upper limit on trees per fold
EARLY_STOP <- 50L    # stop a fold after this many rounds without improvement

# Paths are relative to the project root
objpath <- file.path("results", "model_objects")

# non_predictor_cols, model_variants() and make_cb_pool(), shared with ce_02 and ce_02a
source(file.path("src", "ele", "ce_00_shared.R"))

# ---------------------------------------------------------------------------
# Load train set, split into holdout and training folds
# ---------------------------------------------------------------------------
dtbin    <- readRDS(file.path(objpath, "dtbin.rds"))
train_dt <- dtbin[set == "train"]
rm(dtbin); gc()

# Uses 'cv_fold' column seeded in ce_01
# Fold 1 is the holdout, folds 2 to 5 train and tune (no patient in both!)
holdout_fold     <- 1L
ablation_holdout <- train_dt[cv_fold == holdout_fold]
ablation_train   <- train_dt[cv_fold != holdout_fold]
tuning_folds     <- sort(unique(ablation_train$cv_fold))

# "Fixed Effect" predictors, plus surgeon and procedure group for the full model
fe_predictors <- setdiff(names(train_dt), non_predictor_cols)
both_cols     <- model_variants(fe_predictors)$both$feature_cols
both_cats     <- model_variants(fe_predictors)$both$cat_names

# Inverse-frequency weight: a transfused case counts as much as all non-transfused cases per transfused one
y_ablation_train <- as.numeric(ablation_train$transfused)
pos_weight       <- sum(y_ablation_train == 0) / sum(y_ablation_train == 1)

# Sanity check
message(sprintf("Ablation train prevalence: %.3f%%  (pos_weight = %.2f)", 100 * mean(y_ablation_train), pos_weight))

# ---------------------------------------------------------------------------
# Pools
# ---------------------------------------------------------------------------
# All of ablation_train, for the final fit of each arm
train_pool_both <- make_cb_pool(ablation_train, both_cols, both_cats, y_ablation_train)

# Patient-grouped tuning folds within ablation_train (built once, shared by both arms)
fold_pools <- lapply(
    tuning_folds,
    function(k) {
        in_val <- ablation_train$cv_fold == k
        list(
            train = make_cb_pool(ablation_train[!in_val], both_cols, both_cats, y_ablation_train[!in_val]),
            val   = make_cb_pool(ablation_train[in_val],  both_cols, both_cats, y_ablation_train[in_val]),
            y_val = y_ablation_train[in_val]
        )
    }
)

# ---------------------------------------------------------------------------
# Shared grid ensures comparability
# ---------------------------------------------------------------------------
# Same ranges as ce_02a
set.seed(42)
grid_cb <- data.frame(
    learning_rate = 10^runif(N_RANDOM, -2.0, -0.5),
    depth         = sample(4:10,           N_RANDOM, replace = TRUE),
    l2_leaf_reg   = sample(c(1, 3, 5, 10), N_RANDOM, replace = TRUE)
)

base_params_cb <- list(
    loss_function = "Logloss",
    eval_metric   = "PRAUC:use_weights=false", # important!
    random_seed   = 123,
    thread_count  = n_cores,
    logging_level = "Silent",
    od_type       = "Iter",
    od_wait       = EARLY_STOP,
    iterations    = MAX_ROUNDS
)

# ---------------------------------------------------------------------------
# Helper: tune the full model over the shared grid
# ---------------------------------------------------------------------------
# Same procedure as ce_02a's tune_cb(), plus the arm's class weights
# A draw's score is the mean unweighted validation AUPRC (PRROC), its tree count the mean best iteration
tune_cb_arm <- function(class_weights, arm_label) {
    message("Tuning CB-both [", arm_label, "] (", N_RANDOM, " random draws, ", length(fold_pools), " patient-grouped folds)")

    results <- lapply(
        seq_len(N_RANDOM),
        function(i) {
            params <- c(
                base_params_cb,
                list(
                    learning_rate = grid_cb$learning_rate[i],
                    depth         = grid_cb$depth[i],
                    l2_leaf_reg   = grid_cb$l2_leaf_reg[i],
                    class_weights = class_weights
                )
            )

            # use_best_model cuts each fold's model back to its best iteration
            fold_params <- c(params, list(use_best_model = TRUE, allow_writing_files = FALSE))

            # Columns are folds, rows the tree count and validation AUPRC
            per_fold <- vapply(
                fold_pools,
                function(fp) {
                    model <- catboost.train(fp$train, fp$val, params = fold_params)
                    pred  <- catboost.predict(model, fp$val, prediction_type = "Probability")
                    c(
                        n_trees = model$tree_count,
                        auprc   = pr.curve(scores.class0 = pred, weights.class0 = fp$y_val)$auc.integral
                    )
                },
                numeric(2)
            )

            nrounds    <- as.integer(round(mean(per_fold["n_trees", ])))
            best_score <- mean(per_fold["auprc", ])
            message(sprintf(
                "  draw %02d/%d: nrounds=%d  prauc=%.4f  lr=%.4f  depth=%d",
                i,
                N_RANDOM,
                nrounds,
                best_score,
                grid_cb$learning_rate[i],
                grid_cb$depth[i]
            ))
            list(params = params, nrounds = nrounds, cv_score = best_score)
        }
    )

    best <- results[[which.max(sapply(results, `[[`, "cv_score"))]]
    message(sprintf("CB-both [%s] best: nrounds=%d  prauc=%.4f", arm_label, best$nrounds, best$cv_score))
    best
}

# Weight 1 for non-transfused, pos_weight for transfused (or 1 for both)
best_weighted   <- tune_cb_arm(c(1, pos_weight), "weighted")
best_unweighted <- tune_cb_arm(c(1, 1),          "unweighted")

rm(fold_pools); gc()

# ---------------------------------------------------------------------------
# Fit each arm's best params on ablation_train, predict on ablation_holdout,
# compare holdout PRAUC (ranking), Brier score, and mean predicted
# vs. observed rate (calibration).
# ---------------------------------------------------------------------------
# The label isn't used for prediction, it just comes along with the pool
y_holdout    <- as.numeric(ablation_holdout$transfused)
holdout_pool <- make_cb_pool(ablation_holdout, both_cols, both_cats, y_holdout)

fit_and_evaluate <- function(best, arm_label) {
    params            <- best$params
    params$iterations <- best$nrounds
    model <- catboost.train(train_pool_both, params = params)
    pred  <- catboost.predict(model, holdout_pool, prediction_type = "Probability")

    holdout_pr <- pr.curve(scores.class0 = pred, weights.class0 = y_holdout, curve = FALSE)

    # Weighting inflates predicted risks, which shows up as mean_pred far above observed_rate
    brier <- mean((pred - y_holdout)^2)
    data.table(
        arm           = arm_label,
        cv_prauc      = best$cv_score,
        holdout_prauc = holdout_pr$auc.integral,
        mean_pred     = mean(pred),
        observed_rate = mean(y_holdout),
        brier         = brier
    )
}

ablation_summary <- rbind(
    fit_and_evaluate(best_weighted,   "weighted"),
    fit_and_evaluate(best_unweighted, "unweighted")
)

# Print results
message("\n===== ABLATION SUMMARY (CB-both) =====")
print(ablation_summary)
message("Decision metric is holdout_prauc (external, unweighted)")
message("=======================================\n")

# Save results for reference (Supplementary Table S3)
saveRDS(ablation_summary, file.path(objpath, paste0("ablation_weighting_", model_tag, ".rds")))
fwrite(ablation_summary, file.path(objpath, paste0("ablation_weighting_", model_tag, ".csv")))
message("Saved ablation_weighting_", model_tag, ".rds/.csv to ", objpath)
