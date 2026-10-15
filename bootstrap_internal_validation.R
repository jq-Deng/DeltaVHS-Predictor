


library(rms)
library(pROC)
library(boot)


N_BOOTSTRAP <- 1000


CANDIDATE_VARS <- c(
  "pctDeltaVHS_low",
  "Age",
  "Sex",
  "BMI",
  "Disease_location",
  "Disease_behavior",
  "CDAI",
  "SES_CD",
  "CRP",
  "ESR",
  "Platelet",
  "NLR",
  "Albumin",
  "Creatinine",
  "SII"
)


.IBD_DATA_DIR <- Sys.getenv("IBD_DATA_DIR", unset = "")
if (nzchar(.IBD_DATA_DIR)) setwd(.IBD_DATA_DIR)
input_xlsx <- "Dataset_Train_AllFeatures.xlsx"
input_csv <- "Dataset_Train_AllFeatures.csv"
if (file.exists(input_xlsx)) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is needed to load the training XLSX")
  }
  train_data <- as.data.frame(readxl::read_excel(input_xlsx), check.names = FALSE)
} else if (file.exists(input_csv)) {
  train_data <- utils::read.csv(input_csv, check.names = FALSE,
                                stringsAsFactors = FALSE, fileEncoding = "UTF-8")
} else {
  stop("Training XLSX or CSV not found in working directory")
}


required <- c("SLR","VHSpre","VHSpost")
missing <- setdiff(required, names(train_data))
if (length(missing)) stop("Missing required columns: ", paste(missing, collapse=", "))
train_data$SLR <- as.numeric(train_data$SLR)
if (!all(stats::na.omit(train_data$SLR) %in% c(0,1))) stop("SLR must be 0/1")
train_data$VHSpre <- as.numeric(train_data$VHSpre)
train_data$VHSpost <- as.numeric(train_data$VHSpost)

if ("Disease_behavior" %in% names(train_data)) {
  raw_b <- trimws(as.character(train_data$Disease_behavior))
  raw_b <- ifelse(raw_b %in% c("1","2","3"), paste0("B", raw_b), raw_b)
  unexpected <- setdiff(unique(stats::na.omit(raw_b)), c("B1","B2","B3"))
  if (length(unexpected)) stop("Unexpected disease behavior codes: ",
                               paste(unexpected, collapse=", "))
  train_data$Disease_behavior <- factor(raw_b, levels = c("B1","B2","B3"))
}
for (name in c("Sex","Disease_location")) {
  if (name %in% names(train_data)) train_data[[name]] <- factor(train_data[[name]])
}


if ("Smoking" %in% names(train_data)) {
  CANDIDATE_VARS <- c(CANDIDATE_VARS, "Smoking")
  train_data$Smoking <- factor(train_data$Smoking)
}


needed <- unique(c("SLR", "VHSpre", "VHSpost", CANDIDATE_VARS))
train_data <- train_data[, intersect(needed, names(train_data)), drop = FALSE]


compute_pct_delta_vhs <- function(vhs_pre, vhs_post) {

  vhs_pre_eff <- pmax(vhs_pre, 1e-4)
  return((vhs_pre_eff - vhs_post) / vhs_pre_eff * 100)
}


determine_youden_cutoff <- function(pct_delta, slr) {

  roc_obj <- pROC::roc(slr, pct_delta, levels = c(0,1), direction = ">", quiet = TRUE)
  youden_idx <- which.max(roc_obj$sensitivities + roc_obj$specificities - 1)
  cutoff <- roc_obj$thresholds[youden_idx]
  return(cutoff)
}


dichotomize_pctDelta <- function(pct_delta, cutoff) {

  factor(ifelse(pct_delta <= cutoff, "low", "high"),
         levels = c("high", "low"))
}


get_lrt_p_value <- function(var, data) {


  formula_null <- as.formula("SLR ~ 1")
  formula_full <- as.formula(paste("SLR ~", var))
  tryCatch({
    fit0 <- glm(formula_null, data = data, family = binomial)
    fit1 <- glm(formula_full, data = data, family = binomial)
    lrt_result <- anova(fit0, fit1, test = "Chisq")
    p_val <- lrt_result$`Pr(>Chi)`[2]
    return(p_val)
  }, error = function(e) {
    return(NA)
  })
}


fit_model <- function(data, cutoff) {

  data$pctDeltaVHS_low <- dichotomize_pctDelta(data$pctDeltaVHS, cutoff)


  if ("Disease_behavior" %in% names(data)) {
    data$Disease_behavior <- factor(data$Disease_behavior, levels = c("B1", "B2", "B3"))
  }
  if ("Disease_location" %in% names(data)) {
    data$Disease_location <- factor(data$Disease_location)
  }
  if ("Sex" %in% names(data)) {
    data$Sex <- factor(data$Sex)
  }


  retained <- c("pctDeltaVHS_low")
  for (var in CANDIDATE_VARS) {
    if (var == "pctDeltaVHS_low") next
    if (!var %in% names(data)) next

    if (length(unique(data[[var]])) < 2) next

    p_val <- get_lrt_p_value(var, data)
    if (!is.na(p_val) && p_val < 0.20) {
      retained <- c(retained, var)
    }
  }


  if (length(retained) > 1) {
    full_formula <- as.formula(paste("SLR ~", paste(retained, collapse = " + ")))
    tryCatch({
      full_fit <- glm(full_formula, data = data, family = binomial)

      lower_formula <- as.formula("SLR ~ pctDeltaVHS_low")
      backward_fit <- step(full_fit, scope = list(lower = lower_formula),
                           direction = "backward", trace = 0)


      sm <- summary(backward_fit)
      coef_table <- sm$coefficients
      significant_vars <- c("pctDeltaVHS_low")
      for (var in retained) {
        if (var == "pctDeltaVHS_low") next

        var_coef_rows <- grep(paste0("^", var), rownames(coef_table), value = TRUE)
        if (length(var_coef_rows) > 0) {
          var_p_values <- coef_table[var_coef_rows, "Pr(>|z|)"]
          if (any(!is.na(var_p_values) & var_p_values < 0.05)) {
            significant_vars <- c(significant_vars, var)
          }
        }
      }
      retained <- unique(significant_vars)
    }, error = function(e) {

      warning("AIC backward elimination failed; retaining screened variables")
    })
  }


  if (!"pctDeltaVHS_low" %in% retained) {
    retained <- c(retained, "pctDeltaVHS_low")
  }


  dd <- rms::datadist(data[, c("SLR", retained), drop = FALSE])
  options(datadist = "dd")
  final_formula <- as.formula(paste("SLR ~", paste(retained, collapse = " + ")))
  final_fit <- tryCatch(
    rms::lrm(final_formula, data = data, x = TRUE, y = TRUE),
    error = function(e) {
      warning("lrm failed; using glm with the same formula: ", conditionMessage(e))
      stats::glm(final_formula, data = data, family = binomial)
    }
  )


  return(list(fit = final_fit, vars = retained, cutoff = cutoff))
}


predict_slr_risk <- function(model, newdata = NULL) {
  if (inherits(model, "lrm")) {
    p <- if (is.null(newdata)) predict(model, type = "fitted") else
      predict(model, newdata = newdata, type = "fitted")
  } else if (inherits(model, "glm")) {
    p <- if (is.null(newdata)) stats::predict(model, type = "response") else
      stats::predict(model, newdata = newdata, type = "response")
  } else stop("Unsupported model class for predictions")
  as.numeric(p)
}


compute_calibration <- function(predictions, observed) {


  logit_pred <- qlogis(pmin(pmax(predictions, 1e-10), 1 - 1e-10))
  cal_data <- data.frame(observed = observed, logit_pred = logit_pred)
  tryCatch({
    cal_fit <- glm(observed ~ logit_pred, data = cal_data, family = binomial)
    slope <- coef(cal_fit)[2]
    intercept <- coef(cal_fit)[1]
    return(list(slope = slope, intercept = intercept))
  }, error = function(e) {
    return(list(slope = NA, intercept = NA))
  })
}


brier_score <- function(predictions, observed) {
  predictions <- pmin(pmax(predictions, 1e-10), 1 - 1e-10)
  mean((predictions - observed)^2)
}


ci95 <- function(x) {
  q <- quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
  as.numeric(q)
}


set.seed(42)


train_data$pctDeltaVHS <- compute_pct_delta_vhs(train_data$VHSpre, train_data$VHSpost)


full_cutoff <- determine_youden_cutoff(train_data$pctDeltaVHS, train_data$SLR)
full_model <- fit_model(train_data, full_cutoff)
full_pred <- predict_slr_risk(full_model$fit)
full_apparent_auc <- as.numeric(pROC::roc(train_data$SLR, full_pred, levels = c(0,1), direction = "<", quiet = TRUE)$auc)


results <- data.frame(
  replicate = 1:N_BOOTSTRAP,
  cutoff = NA,
  apparent_AUC = NA,
  test_AUC_original = NA,
  optimism = NA,
  calibration_slope = NA,
  calibration_intercept = NA,
  brier = NA,
  n_vars = NA,
  has_pctDelta = NA,
  has_disease_behavior = NA
)

selection_freq <- list()

for (b in 1:N_BOOTSTRAP) {

  boot_indices <- sample(seq_len(nrow(train_data)), nrow(train_data), replace = TRUE)
  boot_data <- train_data[boot_indices, ]


  boot_cutoff <- determine_youden_cutoff(boot_data$pctDeltaVHS, boot_data$SLR)
  results$cutoff[b] <- boot_cutoff


  model_result <- fit_model(boot_data, boot_cutoff)
  boot_fit <- model_result$fit
  boot_vars <- model_result$vars


  boot_pred <- predict_slr_risk(boot_fit)
  apparent_auc <- as.numeric(pROC::roc(boot_data$SLR, boot_pred, levels = c(0,1), direction = "<", quiet = TRUE)$auc)


  train_orig <- train_data
  train_orig$pctDeltaVHS_low <- dichotomize_pctDelta(train_orig$pctDeltaVHS, boot_cutoff)

  if ("Disease_behavior" %in% names(train_orig)) train_orig$Disease_behavior <- factor(train_orig$Disease_behavior, levels = c("B1","B2","B3"))
  if ("Disease_location" %in% names(train_orig)) train_orig$Disease_location <- factor(train_orig$Disease_location)
  if ("Sex" %in% names(train_orig)) train_orig$Sex <- factor(train_orig$Sex)
  original_pred <- predict_slr_risk(boot_fit, newdata = train_orig)
  test_auc_original <- as.numeric(pROC::roc(train_data$SLR, original_pred, levels = c(0,1), direction = "<", quiet = TRUE)$auc)


  optimism <- apparent_auc - test_auc_original


  cal <- compute_calibration(original_pred, train_data$SLR)


  brier_val <- brier_score(original_pred, train_data$SLR)


  results$apparent_AUC[b] <- apparent_auc
  results$test_AUC_original[b] <- test_auc_original
  results$optimism[b] <- optimism
  results$calibration_slope[b] <- cal$slope
  results$calibration_intercept[b] <- cal$intercept
  results$brier[b] <- brier_val
  results$n_vars[b] <- length(boot_vars)
  results$has_pctDelta[b] <- "pctDeltaVHS_low" %in% boot_vars
  results$has_disease_behavior[b] <- "Disease_behavior" %in% boot_vars


  for (v in boot_vars) {
    selection_freq[[v]] <- ifelse(is.null(selection_freq[[v]]), 1, selection_freq[[v]] + 1)
  }


  if (b %% 100 == 0) cat(sprintf("Bootstrap %d/%d complete\n", b, N_BOOTSTRAP))
}


mean_optimism <- mean(results$optimism, na.rm = TRUE)
optimism_corrected_auc <- full_apparent_auc - mean_optimism


brier_full <- brier_score(full_pred, train_data$SLR)

slope_ci  <- ci95(results$calibration_slope)
int_ci    <- ci95(results$calibration_intercept)
brier_ci  <- ci95(results$brier)

summary_df <- data.frame(
  Metric = c(
    "n_bootstrap_valid",
    "full_data_apparent_AUC",
    "apparent_AUC_mean",
    "apparent_AUC_SD",
    "test_AUC_original_mean",
    "test_AUC_original_SD",
    "optimism_mean",
    "optimism_SD",
    "optimism_corrected_AUC",
    "calibration_slope_mean",
    "calibration_slope_SD",
    "calibration_slope_95CI_low",
    "calibration_slope_95CI_high",
    "calibration_intercept_mean",
    "calibration_intercept_SD",
    "calibration_intercept_95CI_low",
    "calibration_intercept_95CI_high",
    "brier_score_full",
    "brier_score_mean",
    "brier_score_95CI_low",
    "brier_score_95CI_high",
    "mean_n_retained_vars",
    "primary_final_vars"
  ),
  Value = c(
    N_BOOTSTRAP,
    round(full_apparent_auc, 6),
    round(mean(results$apparent_AUC, na.rm = TRUE), 6),
    round(sd(results$apparent_AUC, na.rm = TRUE), 6),
    round(mean(results$test_AUC_original, na.rm = TRUE), 6),
    round(sd(results$test_AUC_original, na.rm = TRUE), 6),
    round(mean_optimism, 6),
    round(sd(results$optimism, na.rm = TRUE), 6),
    round(optimism_corrected_auc, 6),
    round(mean(results$calibration_slope, na.rm = TRUE), 6),
    round(sd(results$calibration_slope, na.rm = TRUE), 6),
    round(slope_ci[1], 6),
    round(slope_ci[2], 6),
    round(mean(results$calibration_intercept, na.rm = TRUE), 6),
    round(sd(results$calibration_intercept, na.rm = TRUE), 6),
    round(int_ci[1], 6),
    round(int_ci[2], 6),
    round(brier_full, 4),
    round(mean(results$brier, na.rm = TRUE), 4),
    round(brier_ci[1], 4),
    round(brier_ci[2], 4),
    round(mean(results$n_vars, na.rm = TRUE), 3),
    paste(full_model$vars, collapse = ", ")
  )
)


sel_freq_df <- data.frame(
  Variable = names(selection_freq),
  Selection_Count = unlist(selection_freq),
  Selection_Frequency_pct = round(unlist(selection_freq) / N_BOOTSTRAP * 100, 1)
)
rownames(sel_freq_df) <- NULL
sel_freq_df <- sel_freq_df[order(-sel_freq_df$Selection_Count), ]


write.csv(results, "Comment9_Bootstrap_Results.csv", row.names = FALSE)
write.csv(summary_df, "Comment9_Bootstrap_Summary.csv", row.names = FALSE)
write.csv(sel_freq_df, "Comment9_Bootstrap_SelectionFrequency.csv", row.names = FALSE)

cat("\n=== Bootstrap Validation Complete ===\n")
cat(sprintf("Full-data apparent AUC (reconstructed): %.4f\n", full_apparent_auc))
cat(sprintf("Mean optimism: %.4f\n", mean_optimism))
cat(sprintf("Optimism-corrected AUC: %.4f\n", optimism_corrected_auc))
cat(sprintf("Mean calibration slope: %.4f (95%% CI %.3f-%.3f)\n", mean(results$calibration_slope, na.rm = TRUE), slope_ci[1], slope_ci[2]))
cat(sprintf("Mean calibration intercept: %.4f (95%% CI %.3f-%.3f)\n", mean(results$calibration_intercept, na.rm = TRUE), int_ci[1], int_ci[2]))
cat(sprintf("Brier score (full-data): %.4f (95%% CI %.3f-%.3f across resamples)\n", brier_full, brier_ci[1], brier_ci[2]))
