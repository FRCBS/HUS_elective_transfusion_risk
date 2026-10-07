# Tunes the hyperparameters of the four CatBoost models.
#
# Random search with patient-grouped cross-validation on the training set,
# saving the best params per model for ce_02.

library(data.table)
library(catboost)
library(parallel)
library(PRROC)

model_tag <- "final"
n_cores   <- detectCores()

N_RANDOM   <- 30L    # random hyperparameter draws per model
MAX_ROUNDS <- 2000L  # upper limit on trees per fold
EARLY_STOP <- 50L    # stop a fold after this many rounds without improvement

# Paths are relative to the project root
objpath <- file.path("results", "model_objects")

# non_predictor_cols, model_variants() and make_cb_pool(), shared with ce_02
source(file.path("src", "ele", "ce_00_shared.R"))

# ---------------------------------------------------------------------------
# Load training set
# ---------------------------------------------------------------------------
# dtbin is the modelling dataset from ce_01, train/test split in `set`
dtbin    <- readRDS(file.path(objpath, "dtbin.rds"))
train_dt <- dtbin[set == "train"]
rm(dtbin); gc()

# Patient-level predictors, the same as in ce_02
fe_predictors <- setdiff(names(train_dt), non_predictor_cols)

# Folds come from ce_01's cv_fold, so no patient is in both a training and a validation fold
# (catboost.cv() can't take custom folds, hence the manual loop below)
n_folds <- max(train_dt$cv_fold)

# Outcome as 0/1 (the label CatBoost expects)
y_train <- as.numeric(train_dt$transfused)

# ---------------------------------------------------------------------------
# Fold pools
# ---------------------------------------------------------------------------
# One training and one validation pool per fold
# Built once per model and reused by every draw (pool building is slow!)
make_fold_pools <- function(feature_cols, cat_names) {
    lapply(
        seq_len(n_folds),
        function(k) {
            in_val <- train_dt$cv_fold == k
            list(
                train = make_cb_pool(train_dt[!in_val], feature_cols, cat_names, y_train[!in_val]),
                val   = make_cb_pool(train_dt[in_val],  feature_cols, cat_names, y_train[in_val]),
                y_val = y_train[in_val]
            )
        }
    )
}

# ---------------------------------------------------------------------------
# Search space
# ---------------------------------------------------------------------------
# Settings shared by every draw
# Early stopping watches CatBoost's own PRAUC on the validation fold (no class weights, so plain AUPRC)
fixed_params_cb <- list(
    loss_function = "Logloss",
    eval_metric   = "PRAUC",
    random_seed   = 123,
    thread_count  = n_cores,
    logging_level = "Silent",
    od_type       = "Iter",
    od_wait       = EARLY_STOP,
    iterations    = MAX_ROUNDS
)

# The same N_RANDOM draws for all four models, so they're tuned over an identical grid
# Learning rate log-uniform from 0.01 to about 0.32, tree depth 4 to 10, L2 regularisation from a set
set.seed(42)
grid_cb <- data.frame(
    learning_rate = 10^runif(N_RANDOM, -2.0, -0.5),
    depth         = sample(4:10,           N_RANDOM, replace = TRUE),
    l2_leaf_reg   = sample(c(1, 3, 5, 10), N_RANDOM, replace = TRUE)
)

# ---------------------------------------------------------------------------
# Tuning
# ---------------------------------------------------------------------------
# Scores every draw over the folds and returns the best one
# A draw's score is the mean validation AUPRC (PRROC, as in ce_03), its tree count the mean best iteration
tune_cb <- function(fold_pools, variant_label) {
    message("Tuning model_cb_", variant_label, " (", N_RANDOM, " draws, ", n_folds, " folds)")

    results <- lapply(
        seq_len(N_RANDOM),
        function(i) {
            params <- c(
                fixed_params_cb,
                list(
                    learning_rate = grid_cb$learning_rate[i],
                    depth         = grid_cb$depth[i],
                    l2_leaf_reg   = grid_cb$l2_leaf_reg[i]
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
                "  draw %02d/%d: nrounds=%d  prauc=%.4f (fold sd %.4f)  lr=%.4f  depth=%d",
                i,
                N_RANDOM,
                nrounds,
                best_score,
                sd(per_fold["auprc", ]),
                grid_cb$learning_rate[i],
                grid_cb$depth[i]
            ))

            # fold_scores feed the fold SD in Supplementary Table S10
            list(params = params, nrounds = nrounds, cv_score = best_score, fold_scores = per_fold["auprc", ])
        }
    )

    best <- results[[which.max(sapply(results, `[[`, "cv_score"))]]
    message(sprintf("model_cb_%s best: nrounds=%d  prauc=%.4f", variant_label, best$nrounds, best$cv_score))
    best
}

# ---------------------------------------------------------------------------
# Tune and save
# ---------------------------------------------------------------------------
variants <- model_variants(fe_predictors)

for (v in names(variants)) {
    fold_pools <- make_fold_pools(variants[[v]]$feature_cols, variants[[v]]$cat_names)
    best       <- tune_cb(fold_pools, v)
    saveRDS(best, file.path(objpath, paste0("tuned_params_cb_", v, "_", model_tag, ".rds")))

    # Fold pools are large, free them before the next model
    rm(fold_pools); gc()
}

message("All tuning complete. Results saved to: ", objpath)
