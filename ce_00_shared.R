# Definitions shared by the model-fitting scripts (ce_02a, ce_02, ablation_weighting).
#
# Sourced, not run on its own, so the tuned and the final models can't drift apart.

# Columns in dtbin that are never patient-level predictors
non_predictor_cols <- c(
    "transfused",                        # the outcome
    "set",                               # train/test label
    "toimenpide_tekija_id",              # surgeon, added only to the models that use it
    "toimenpide_tekija_ammatti_selite",  # surgeon title (descriptive, ce_03)
    "procedure_group",                   # added only to the models that use it
    "patient_idx",                       # patient id for the cluster bootstrap
    "cv_fold"                            # tuning fold id
)

# The four models of the 2 x 2 design: same patient predictors, differing only in the grouping variables
# Sex and ASA class are the categorical patient predictors
model_variants <- function(fe_predictors) {
    list(
        patient = list(
            feature_cols = fe_predictors,
            cat_names    = c("sukupuoli_selite", "asa")
        ),
        proc = list(
            feature_cols = c(fe_predictors, "procedure_group"),
            cat_names    = c("procedure_group", "sukupuoli_selite", "asa")
        ),
        surgeon = list(
            feature_cols = c(fe_predictors, "toimenpide_tekija_id"),
            cat_names    = c("toimenpide_tekija_id", "sukupuoli_selite", "asa")
        ),
        both = list(
            feature_cols = c(fe_predictors, "toimenpide_tekija_id", "procedure_group"),
            cat_names    = c("toimenpide_tekija_id", "procedure_group", "sukupuoli_selite", "asa")
        )
    )
}

# CatBoost wants categorical features as factor columns in a data.frame
make_cb_pool <- function(dt, feature_cols, cat_names, label) {
    dt_pool <- copy(dt[, ..feature_cols])
    dt_pool[, (cat_names) := lapply(.SD, as.factor), .SDcols = cat_names]
    catboost.load_pool(data = as.data.frame(dt_pool), label = label)
}
