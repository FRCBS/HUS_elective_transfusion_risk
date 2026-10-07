# Fits the four CatBoost models on the full training set.
#
# Uses the hyperparameters tuned in ce_02a and saves one bundle per model
# (the model plus the feature names it needs at prediction time).

library(data.table)
library(catboost)
library(parallel)

model_tag <- "final"
n_cores   <- detectCores()

# Paths are relative to the project root
objpath <- file.path("results", "model_objects")

# non_predictor_cols, model_variants() and make_cb_pool(), shared with ce_02a
source(file.path("src", "ele", "ce_00_shared.R"))

# ---------------------------------------------------------------------------
# Load training set
# ---------------------------------------------------------------------------
# dtbin is the modelling dataset from ce_01, train/test split in `set`
dtbin    <- readRDS(file.path(objpath, "dtbin.rds"))
train_dt <- dtbin[set == "train"]

# Only the training set is needed here, free the rest
rm(dtbin); gc()

# Outcome as 0/1 (the label CatBoost expects)
y_train <- as.numeric(train_dt$transfused)

# ---------------------------------------------------------------------------
# Predictors
# ---------------------------------------------------------------------------
# Patient-level predictors, shared by all four models
fe_predictors <- setdiff(names(train_dt), non_predictor_cols)

# Saved for ce_04, ce_05 and ce_05b
saveRDS(fe_predictors, file.path(objpath, paste0("fe_predictors_", model_tag, ".rds")))

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
# Best params from ce_02a, with the number of trees set to the mean best iteration across folds
# (stops here if ce_02a hasn't been run!)
load_tuned_cb <- function(variant) {
    tuned <- readRDS(file.path(objpath, paste0("tuned_params_cb_", variant, "_", model_tag, ".rds")))

    params               <- tuned$params
    params$iterations    <- tuned$nrounds
    params$random_seed   <- 123
    params$thread_count  <- n_cores
    params$logging_level <- "Silent"

    # The early-stopping settings from tuning do nothing here (no validation set)
    params
}

# Model files are named by variant, e.g. model_cb_both_final.rds
model_file <- function(variant) {
    file.path(objpath, paste0("model_cb_", variant, "_", model_tag, ".rds"))
}

# Feature columns and categorical columns of the four models
variants <- model_variants(fe_predictors)

# ---------------------------------------------------------------------------
# Fit and save
# ---------------------------------------------------------------------------
for (v in names(variants)) {
    message("Fitting model_cb_", v)

    feature_cols <- variants[[v]]$feature_cols
    cat_names    <- variants[[v]]$cat_names

    pool  <- make_cb_pool(train_dt, feature_cols, cat_names, y_train)
    model <- catboost.train(pool, params = load_tuned_cb(v))

    # Feature names travel with the model, so later scripts build the test pool identically
    saveRDS(
        list(
            model             = model,
            feature_names     = feature_cols,
            cat_feature_names = cat_names
        ),
        model_file(v)
    )
}

message("Done.")
