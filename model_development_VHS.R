


pkgs <- c("pROC","dplyr","readxl","gtsummary","rms","dcurves","ggplot2","gridExtra","tibble","rlang")
for(p in pkgs){
  if(!require(p, character.only = TRUE)) stop(paste("Missing required R package:", p))
  library(p, character.only = TRUE)
}


.IBD_DATA_DIR <- Sys.getenv("IBD_DATA_DIR", unset = "")
if (nzchar(.IBD_DATA_DIR)) setwd(.IBD_DATA_DIR)

OUT_DIR   <- file.path(getwd(), "Outputs_Nomogram_IBD")
TABLE_DIR <- file.path(OUT_DIR, "tables")
FIG_DIR   <- file.path(OUT_DIR, "figures")
dir.create(OUT_DIR,   showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR,   showWarnings = FALSE, recursive = TRUE)

cat("\n[INFO] Output folders:\n")
cat(" - ", basename(TABLE_DIR), "\n")
cat(" - ", basename(FIG_DIR), "\n\n")


export_lrm_to_csv <- function(lrm_fit, out_csv, digits = 4) {
  beta <- stats::coef(lrm_fit)
  V <- stats::vcov(lrm_fit)
  se <- sqrt(diag(V))

  common <- intersect(names(beta), names(se))
  beta <- beta[common]
  se   <- se[common]

  zval <- beta / se
  pval <- 2 * (1 - stats::pnorm(abs(zval)))

  z <- stats::qnorm(0.975)
  lower <- beta - z * se
  upper <- beta + z * se

  out <- data.frame(
    term    = names(beta),
    beta    = round(as.numeric(beta), digits),
    SE      = round(as.numeric(se), digits),
    OR      = round(exp(as.numeric(beta)), digits),
    CI_low  = round(exp(as.numeric(lower)), digits),
    CI_high = round(exp(as.numeric(upper)), digits),
    z       = round(as.numeric(zval), digits),
    p_value = signif(as.numeric(pval), 3),
    stringsAsFactors = FALSE
  )

  out <- out[!out$term %in% c("Intercept","(Intercept)"), , drop = FALSE]
  utils::write.csv(out, out_csv, row.names = FALSE, fileEncoding = "UTF-8")
  out
}


.extract_rms_calibrate_curve <- function(cal_obj, verbose = TRUE) {
  m <- as.matrix(cal_obj)
  cn <- colnames(m)
  if(is.null(cn)) stop("calibrate() returned no column names.")

  pred_col <- if("predy" %in% cn) "predy" else grep("pred", cn, ignore.case = TRUE, value = TRUE)[1]
  corr_col <- if("calibrated.corrected" %in% cn) "calibrated.corrected" else grep("corrected", cn, ignore.case = TRUE, value = TRUE)[1]
  if(is.na(pred_col) || is.na(corr_col)) stop("Unable to identify predy / calibrated.corrected columns.")

  if(verbose) message(sprintf("[calibrate] pred=%s | corrected=%s", pred_col, corr_col))

  curve <- data.frame(
    pred = as.numeric(m[, pred_col]),
    obs  = as.numeric(m[, corr_col])
  ) %>%
    dplyr::filter(is.finite(pred) & is.finite(obs)) %>%
    dplyr::arrange(pred)

  curve$pred <- pmin(pmax(curve$pred, 0), 1)
  curve$obs  <- pmin(pmax(curve$obs, 0), 1)
  curve
}


.make_binned_calibration <- function(df, pred_col = "nom_prob", n_bins = 10) {
  tmp <- df %>%
    dplyr::select(group, dplyr::all_of(pred_col)) %>%
    dplyr::filter(!is.na(group) & is.finite(.data[[pred_col]])) %>%
    as.data.frame()

  if(nrow(tmp) == 0) {
    return(data.frame(pred=numeric(0), obs=numeric(0), lower=numeric(0), upper=numeric(0), n=integer(0)))
  }

  probs <- tmp[[pred_col]]
  qs <- quantile(probs, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE, type = 8)
  qs[1] <- 0
  qs[length(qs)] <- 1

  brks <- unique(as.numeric(qs))
  if(length(brks) < 3){
    brks <- seq(0, 1, length.out = min(n_bins, 5) + 1)
  }

  tmp$binq <- cut(probs, breaks = brks, include.lowest = TRUE)

  out <- tmp %>%
    dplyr::group_by(binq) %>%
    dplyr::summarise(
      pred = mean(.data[[pred_col]], na.rm = TRUE),
      obs  = mean(group == 1, na.rm = TRUE),
      n    = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(pred) & is.finite(obs)) %>%
    dplyr::arrange(pred)

  z <- 1.96
  p <- out$obs
  n <- out$n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2*n)) / denom
  half <- (z * sqrt((p*(1-p)/n) + (z^2/(4*n^2)))) / denom

  out$lower <- pmax(0, center - half)
  out$upper <- pmin(1, center + half)
  out
}


.make_loess_calibration <- function(df, pred_col = "nom_prob", span = 0.9, grid_n = 200) {
  tmp <- df %>%
    dplyr::select(group, dplyr::all_of(pred_col)) %>%
    dplyr::filter(!is.na(group) & is.finite(.data[[pred_col]])) %>%
    as.data.frame()

  if(nrow(tmp) < 10){
    return(data.frame(pred = numeric(0), obs = numeric(0)))
  }

  fit <- tryCatch(
    stats::loess(
      as.formula(paste("group ~", pred_col)),
      data = tmp,
      span = span,
      degree = 1
    ),
    error = function(e) NULL
  )

  if(is.null(fit)){
    return(data.frame(pred = numeric(0), obs = numeric(0)))
  }

  x_grid <- seq(
    max(0, min(tmp[[pred_col]], na.rm = TRUE)),
    min(1, max(tmp[[pred_col]], na.rm = TRUE)),
    length.out = grid_n
  )

  newd <- data.frame(tmp_pred = x_grid)
  names(newd) <- pred_col

  y_hat <- tryCatch(as.numeric(predict(fit, newdata = newd)),
                    error = function(e) rep(NA_real_, length(x_grid)))

  out <- data.frame(pred = x_grid, obs = y_hat) %>%
    dplyr::filter(is.finite(pred) & is.finite(obs))

  out$obs <- pmin(pmax(out$obs, 0), 1)
  out
}


align_factor_levels_for_rms <- function(train_df, test_df, ref_df, vars) {

  fac_vars <- vars[sapply(vars, function(v) v %in% colnames(ref_df) && is.factor(ref_df[[v]]))]
  level_map <- lapply(fac_vars, function(v) levels(ref_df[[v]]))
  names(level_map) <- fac_vars

  align_one <- function(df){
    out <- as.data.frame(df)

    for(v in names(level_map)){
      if(v %in% colnames(out)){
        out[[v]] <- trimws(as.character(out[[v]]))
        out[[v]] <- factor(out[[v]], levels = level_map[[v]])
      }
    }
    out
  }

  list(train = align_one(train_df), test = align_one(test_df))
}


.calc_auc_ci <- function(y, p){
  df <- data.frame(y = y, p = p)
  df <- df[!is.na(df$y) & is.finite(df$p), , drop = FALSE]
  roc_obj <- pROC::roc(df$y, df$p, levels = c(0,1), direction = "<", quiet = TRUE)
  ci_obj  <- pROC::ci.auc(roc_obj)
  list(roc = roc_obj, auc = as.numeric(pROC::auc(roc_obj)), ci = ci_obj, n = nrow(df))
}


train_data <- readxl::read_excel("Dataset_Train_AllFeatures.xlsx")
test_data  <- readxl::read_excel("Dataset_Test_AllFeatures.xlsx")

process_data <- function(data) {
  required_columns <- c(
    "VHSpre", "VHSpost", "SLR", "Age", "Sex", "BMI",
    "Disease_location", "Disease_behavior", "CDAI", "SES_CD", "CRP",
    "ESR", "Platelet", "NLR", "Albumin", "Creatinine", "SII"
  )
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns)) {
    stop("Missing required English column(s): ", paste(missing_columns, collapse = ", "))
  }

  data <- data %>%
    dplyr::rename(
      Baseline_ITH = VHSpre,
      Post_ITH = VHSpost,
      group = SLR
    ) %>%
    dplyr::mutate(
      pct_Delta_ITH = (pmax(Baseline_ITH, 1e-4) - Post_ITH) /
        pmax(Baseline_ITH, 1e-4) * 100
    )

  data$group <- as.numeric(data$group)
  data$pct_Delta_ITH[!is.finite(data$pct_Delta_ITH)] <- NA
  as.data.frame(data)
}

train_data <- process_data(train_data)
test_data  <- process_data(test_data)


roc_obj_cutoff <- pROC::roc(train_data$group, train_data$pct_Delta_ITH, levels = c(0,1), direction = ">", quiet = TRUE)
best_cutoff <- as.numeric(pROC::coords(roc_obj_cutoff, "best", ret="threshold", best.method="youden")[1])
cat(sprintf("\n=== Cutoff pct_Delta_ITH: %.2f%% ===\n", best_cutoff))


train_data$pct_Delta_ITH_cat <- ifelse(train_data$pct_Delta_ITH >= best_cutoff, 1, 0)
test_data$pct_Delta_ITH_cat  <- ifelse(test_data$pct_Delta_ITH  >= best_cutoff, 1, 0)
train_data$pct_Delta_ITH_cat <- factor(train_data$pct_Delta_ITH_cat, levels = c(1,0))
test_data$pct_Delta_ITH_cat  <- factor(test_data$pct_Delta_ITH_cat,  levels = c(1,0))


base_factor_vars <- c("Sex","Disease_location","Disease_behavior","pct_Delta_ITH_cat")
if("Smoking" %in% colnames(train_data)) base_factor_vars <- c(base_factor_vars, "Smoking")

for(v in base_factor_vars){
  if(v %in% colnames(train_data)){
    train_data[[v]] <- factor(trimws(as.character(train_data[[v]])))
    test_data[[v]]  <- factor(trimws(as.character(test_data[[v]])), levels = levels(train_data[[v]]))
  }
}


pkgs <- c("dplyr","gtsummary","tibble","rlang","flextable","officer")
for(p in pkgs){
  if(!require(p, character.only = TRUE)) stop(paste("Missing required R package:", p))
  library(p, character.only = TRUE)
}


table1_vars <- c("Age","Sex","BMI","Disease_location","Disease_behavior","CDAI","SES_CD",
                 "CRP","ESR","Platelet","NLR","Albumin","Creatinine","SII",
                 "pct_Delta_ITH","pct_Delta_ITH_cat")
if("Smoking" %in% colnames(train_data)) table1_vars <- c(table1_vars, "Smoking")


train_t1 <- train_data %>%
  dplyr::mutate(Cohort = "Training") %>%
  dplyr::select(Cohort, dplyr::all_of(table1_vars))

test_t1 <- test_data %>%
  dplyr::mutate(Cohort = "Testing") %>%
  dplyr::select(Cohort, dplyr::all_of(table1_vars))

t1_all <- dplyr::bind_rows(train_t1, test_t1) %>% as.data.frame()
t1_all$Cohort <- factor(t1_all$Cohort, levels = c("Training","Testing"))


cat_vars <- c("Sex","Disease_location","Disease_behavior","pct_Delta_ITH_cat","Smoking")
for(v in intersect(cat_vars, colnames(t1_all))){
  t1_all[[v]] <- as.factor(t1_all[[v]])
}


.is_normal_two_groups <- function(df, var, group_col = "Cohort", p_cut = 0.05){
  if(!is.numeric(df[[var]])) return(FALSE)
  glev <- levels(df[[group_col]])
  if(length(glev) != 2) return(FALSE)

  pvals <- c()
  for(g in glev){
    x <- df[df[[group_col]] == g, var]
    x <- x[is.finite(x) & !is.na(x)]


    if(length(x) < 3) return(FALSE)
    if(length(x) > 5000) x <- sample(x, 5000)

    pv <- tryCatch(stats::shapiro.test(x)$p.value, error = function(e) NA_real_)
    if(!is.finite(pv)) return(FALSE)
    pvals <- c(pvals, pv)
  }
  all(pvals > p_cut)
}


cont_vars <- setdiff(table1_vars, intersect(table1_vars, cat_vars))
cont_vars <- cont_vars[cont_vars %in% colnames(t1_all)]
cont_vars <- cont_vars[sapply(cont_vars, function(v) is.numeric(t1_all[[v]]))]


normal_vars <- cont_vars[sapply(cont_vars, function(v) .is_normal_two_groups(t1_all, v, "Cohort", 0.05))]
nonnormal_vars <- setdiff(cont_vars, normal_vars)


stat_list <- list()
test_list <- list()

for(v in cont_vars){
  if(v %in% normal_vars){
    stat_list[[v]] <- "{mean} ± {sd}"
    test_list[[v]] <- "t.test"
  } else {
    stat_list[[v]] <- "{median} ({p25}, {p75})"
    test_list[[v]] <- "wilcox.test"
  }
}


stat_list[["all_categorical()"]] <- "{n} ({p}%)"

test_list[["all_categorical()"]] <- "chisq.test"

tbl1_tt <- t1_all %>%
  gtsummary::tbl_summary(
    by = Cohort,
    missing = "ifany",
    statistic = stat_list
  ) %>%
  gtsummary::add_p(test = test_list) %>%
  gtsummary::modify_caption("**Table 1. Baseline Characteristics (Training vs Testing Cohorts)**")


ft_tbl1 <- gtsummary::as_flex_table(tbl1_tt)
docx_path <- file.path(TABLE_DIR, "Table1_Baseline_Train_vs_Test.docx")
flextable::save_as_docx(`Table 1` = ft_tbl1, path = docx_path)

cat("[OK] Word exported:\n - ", basename(docx_path), "\n", sep = "")


cat("\n=== Logistic regression (Univ -> Multiv stepwise -> Significant) ===\n")

candidate_vars <- c("pct_Delta_ITH_cat","Age","Sex","BMI","Disease_location","Disease_behavior",
                    "CDAI","SES_CD","CRP","ESR","Platelet","NLR","Albumin","Creatinine","SII")
if("Smoking" %in% colnames(train_data)) candidate_vars <- c(candidate_vars, "Smoking")


univ_table <- gtsummary::tbl_uvregression(
  data = train_data %>% dplyr::select(all_of(c("group", candidate_vars))),
  method = glm,
  y = group,
  method.args = list(family = binomial),
  exponentiate = TRUE
) %>% gtsummary::bold_p(t = 0.05)

write.csv(tibble::as_tibble(univ_table, col_labels = TRUE),
          file.path(TABLE_DIR, "Univariate_Logistic.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("[OK] Univariate logistic saved.\n")


get_p_value <- function(var_name) {
  fit1 <- glm(as.formula(paste("group ~", var_name)), data = train_data, family = binomial)
  fit0 <- glm(group ~ 1, data = train_data, family = binomial)
  anova(fit0, fit1, test = "Chisq")$`Pr(>Chi)`[2]
}

p_values <- sapply(candidate_vars, function(v) if(v %in% colnames(train_data)) get_p_value(v) else NA)
p_values <- p_values[!is.na(p_values)]
selected_vars <- names(p_values[p_values < 0.20])
if(!("pct_Delta_ITH_cat" %in% selected_vars)) selected_vars <- c("pct_Delta_ITH_cat", selected_vars)

complete_train_data <- train_data %>%
  dplyr::select(group, all_of(selected_vars)) %>%
  na.omit() %>%
  as.data.frame()

full_model  <- glm(group ~ ., data = complete_train_data, family = binomial)
lower_model <- glm(group ~ pct_Delta_ITH_cat, data = complete_train_data, family = binomial)

step_model <- step(full_model,
                   scope = list(lower = formula(lower_model), upper = formula(full_model)),
                   direction = "backward",
                   trace = 0)

multi_table_step <- gtsummary::tbl_regression(step_model, exponentiate = TRUE) %>%
  gtsummary::bold_p(t = 0.05)

write.csv(tibble::as_tibble(multi_table_step, col_labels = TRUE),
          file.path(TABLE_DIR, "Multivariable_Logistic_AIC.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("[OK] Multivariable(AIC) logistic saved.\n")


sm <- summary(step_model)
coef_tab <- as.data.frame(sm$coefficients)
coef_tab$term <- rownames(coef_tab)
coef_tab <- coef_tab[coef_tab$term != "(Intercept)", , drop = FALSE]
sig_terms <- coef_tab$term[coef_tab$`Pr(>|z|)` < 0.05]

collapse_to_vars <- function(terms, vars){
  out <- c()
  for(v in vars){
    hit <- terms[grepl(paste0("^", v), terms)]
    if(length(hit) > 0) out <- c(out, v)
  }
  unique(out)
}
sig_vars <- collapse_to_vars(sig_terms, selected_vars)
if(!("pct_Delta_ITH_cat" %in% sig_vars)) sig_vars <- c("pct_Delta_ITH_cat", sig_vars)

cat("\n[INFO] FINAL vars used in Nomogram (P<0.05) + forced ΔITH:\n")
print(sig_vars)


cat("\n=== Step 3: Final Model & Professional Nomogram ===\n")

final_train_data <- train_data %>%
  dplyr::select(group, all_of(sig_vars)) %>%
  na.omit() %>%
  as.data.frame()


if ("pct_Delta_ITH_cat" %in% colnames(final_train_data)) {
  final_train_data$pct_Delta_ITH_cat <- factor(final_train_data$pct_Delta_ITH_cat, levels = c("1", "0"))
}

if ("Sex" %in% colnames(final_train_data)) {
  final_train_data$Sex <- factor(final_train_data$Sex, levels = c("0", "1"))
}

if ("Disease_location" %in% colnames(final_train_data)) {
  final_train_data$Disease_location <- factor(as.character(final_train_data$Disease_location), levels = c("1", "2", "3"))
}

if ("Disease_behavior" %in% colnames(final_train_data)) {
  final_train_data$Disease_behavior <- factor(as.character(final_train_data$Disease_behavior), levels = c("1", "2", "3"))
}


prof_labels <- c(
  "pct_Delta_ITH_cat" = "%ΔVHS",
  "Age"               = "Age (years)",
  "Sex"               = "Sex",
  "BMI"               = "BMI",
  "Disease_location"  = "Disease Location",
  "Disease_behavior"  = "Disease Behavior",
  "CDAI"              = "CDAI",
  "SES_CD"            = "SES-CD",
  "CRP"               = "CRP (mg/L)",
  "ESR"               = "ESR (mm/h)",
  "Platelet"          = "Platelet",
  "NLR"               = "NLR",
  "Albumin"           = "Albumin (g/L)",
  "Creatinine"        = "Creatinine",
  "SII"               = "SII"
)


for (v in colnames(final_train_data)) {
  if (v %in% names(prof_labels)) {
    Hmisc::label(final_train_data[[v]]) <- prof_labels[[v]]
  }
}


dd <- rms::datadist(final_train_data)
options(datadist = "dd")

final_lrm <- rms::lrm(
  as.formula(paste("group ~", paste(sig_vars, collapse = " + "))),
  data = final_train_data, x = TRUE, y = TRUE
)


final_lrm_csv <- file.path(TABLE_DIR, "Multivariable_Logistic_Final_Significant.csv")
export_lrm_to_csv(final_lrm, final_lrm_csv)
cat("[OK] Final lrm table saved.\n")


clean_train <- train_data[!is.na(train_data$group), ] %>% as.data.frame()
clean_test  <- test_data[!is.na(test_data$group), ] %>% as.data.frame()


aligned <- align_factor_levels_for_rms(clean_train, clean_test, final_train_data, sig_vars)
clean_train <- aligned$train
clean_test  <- aligned$test


clean_train$nom_prob <- as.numeric(predict(final_lrm, newdata = clean_train, type = "fitted"))
clean_test$nom_prob  <- as.numeric(predict(final_lrm, newdata = clean_test, type = "fitted"))


nom <- rms::nomogram(
  final_lrm,
  fun = plogis,
  fun.at = c(0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.95),
  funlabel = "Risk of Loss of Response",
  lp = FALSE
)


pdf(file.path(FIG_DIR, "Nomogram_Professional.pdf"), width = 11, height = 8)
par(family = "serif", mar = c(3, 2, 3, 2), lwd = 1.3)
plot(nom, xfrac = 0.35, cex.var = 1.05, cex.axis = 0.95, tck = -0.015, lmgp = 0.25, col.grid = c("#D3D8DE", "#E9ECEF"))
dev.off()

png(file.path(FIG_DIR, "Nomogram_Professional.png"), width = 3300, height = 2400, res = 300)
par(family = "serif", mar = c(3, 2, 3, 2), lwd = 1.3)
plot(nom, xfrac = 0.35, cex.var = 1.05, cex.axis = 0.95, tck = -0.015, lmgp = 0.25, col.grid = c("#D3D8DE", "#E9ECEF"))
dev.off()

cat("[OK] Professional Nomogram exported.\n")


cat("\n=== Step 4-5: Calibration Curves ===\n")


clean_train$p_bnd <- pmin(pmax(clean_train$nom_prob, 1e-6), 1 - 1e-6)
clean_test$p_bnd  <- pmin(pmax(clean_test$nom_prob,  1e-6), 1 - 1e-6)


cal_train_obj <- tryCatch(
  rms::calibrate(final_lrm, method = "boot", B = 1000),
  error = function(e) rms::calibrate(final_lrm, method = "boot", B = 200)
)
curve_train <- .extract_rms_calibrate_curve(cal_train_obj, verbose = FALSE)
bin_train   <- .make_binned_calibration(clean_train, pred_col = "nom_prob", n_bins = 10)


curve_test <- .make_loess_calibration(clean_test, pred_col = "nom_prob", span = 0.9, grid_n = 200)
bin_test   <- .make_binned_calibration(clean_test, pred_col = "nom_prob", n_bins = 10)


brier_tr <- tryCatch(rms::val.prob(p = clean_train$p_bnd, y = clean_train$group, pl = FALSE, statloc = FALSE)["Brier"], error = function(e) NA)
brier_te <- tryCatch(rms::val.prob(p = clean_test$p_bnd,  y = clean_test$group,  pl = FALSE, statloc = FALSE)["Brier"], error = function(e) NA)


plot_calibration_npj <- function(curve_df, bin_df, rug_df, brier_val, title, tag) {
  ggplot2::ggplot() +

    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.9, color = "#8895A5") +

    ggplot2::geom_line(data = curve_df, ggplot2::aes(x = pred, y = obs), color = "#E64B35", linewidth = 1.3) +

    ggplot2::geom_errorbar(data = bin_df, ggplot2::aes(x = pred, ymin = lower, ymax = upper), width = 0.02, color = "#E64B35", linewidth = 0.8) +
    ggplot2::geom_point(data = bin_df, ggplot2::aes(x = pred, y = obs), color = "#E64B35", size = 2.8) +

    ggplot2::geom_rug(data = rug_df, ggplot2::aes(x = nom_prob), sides = "b", alpha = 0.25, color = "#222222") +

    ggplot2::annotate("text", x = 0.96, y = 0.04, hjust = 1, vjust = 0, size = 4.2, family = "serif", color = "#222222", fontface="bold",
                      label = sprintf("Brier = %.3f", brier_val)) +

    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip="off") +
    ggplot2::labs(x = "Predicted Probability", y = "Observed Probability", title = title, tag = tag) +
    ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "#F4F7F9", color = NA),
      panel.background = ggplot2::element_rect(fill = "#F4F7F9", color = NA),
      panel.grid.major = ggplot2::element_line(color = "#FFFFFF", linewidth = 0.9),
      panel.grid.minor = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 18, color = "black"),
      plot.tag.position = c(0.015, 0.98),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13, color = "#111111", margin=ggplot2::margin(b=14)),
      axis.text = ggplot2::element_text(size = 11, color = "#222222", face = "bold"),
      axis.title = ggplot2::element_text(size = 12, color = "#111111", face = "bold"),
      plot.margin = ggplot2::margin(20, 15, 20, 15)
    )
}


p_cal_train <- plot_calibration_npj(curve_train, bin_train, clean_train, brier_tr, "Calibration (Training cohort)", "A")
p_cal_test  <- plot_calibration_npj(curve_test,  bin_test,  clean_test,  brier_te, "Calibration (Testing cohort)",  "B")

cat("[OK] Calibration plots generated.\n")


cat("\n=== Step 6: Decision Curve Analysis (Advanced Custom Plot) ===\n")


VHS_RELATED <- c("pct_Delta_VHS","pct_Delta_VHS_cat","Delta_VHS","VHS_pre","VHS_post",
                 "pct_Delta_ITH","pct_Delta_ITH_cat","Baseline_ITH","Post_ITH")
vars_clinical_only <- setdiff(sig_vars, intersect(sig_vars, VHS_RELATED))


fit_clin <- stats::glm(as.formula(paste("group ~", paste(vars_clinical_only, collapse = " + "))),
                       data = train_data, family = binomial)

clean_train$score_clin <- as.numeric(predict(fit_clin, newdata = clean_train, type = "response"))
clean_test$score_clin  <- as.numeric(predict(fit_clin, newdata = clean_test, type = "response"))


clean_train$score_dith <- (1 - suppressWarnings(as.numeric(as.character(clean_train$pct_Delta_ITH_cat))))
clean_test$score_dith  <- (1 - suppressWarnings(as.numeric(as.character(clean_test$pct_Delta_ITH_cat))))


dca_train <- dcurves::dca(
  group ~ nom_prob + score_clin + score_dith,
  data = clean_train,
  thresholds = seq(0, 0.90, by = 0.01),
  label = list(nom_prob="Nom-%ΔVHS", score_clin="Clin_pre", score_dith="%ΔVHS")
)

dca_test <- dcurves::dca(
  group ~ nom_prob + score_clin + score_dith,
  data = clean_test,
  thresholds = seq(0, 0.90, by = 0.01),
  label = list(nom_prob="Nom-%ΔVHS", score_clin="Clin_pre", score_dith="%ΔVHS")
)


plot_advanced_dca <- function(dca_obj, title, tag) {

  df <- as.data.frame(dca_obj$dca)


  max_nb <- max(df$net_benefit, na.rm = TRUE)
  y_lower <- -0.02 * max_nb


  df$label <- factor(df$label, levels = c("Nom-%ΔVHS", "Clin_pre", "%ΔVHS", "Treat All", "Treat None"))

  ggplot2::ggplot(df, ggplot2::aes(x = threshold, y = net_benefit, color = label, linetype = label)) +

    ggplot2::geom_hline(yintercept = 0, color = "#7B8B9A", linewidth = 0.6) +

    ggplot2::geom_line(linewidth = 1.3) +

    ggplot2::scale_color_manual(values = c(
      "Nom-%ΔVHS"  = "#E64B35",
      "Clin_pre"   = "#00A087",
      "%ΔVHS"      = "#3C5488",
      "Treat All"  = "#A0A8B4",
      "Treat None" = "#222222"
    )) +
    ggplot2::scale_linetype_manual(values = c(
      "Nom-%ΔVHS"  = "solid",
      "Clin_pre"   = "longdash",
      "%ΔVHS"      = "dotdash",
      "Treat All"  = "dashed",
      "Treat None" = "solid"
    )) +

    ggplot2::coord_cartesian(xlim = c(0, 0.9), ylim = c(y_lower, max_nb * 1.05), expand = FALSE) +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::labs(title = title, x = "Threshold Probability", y = "Net Benefit", tag = tag) +

    ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "#F4F7F9", color = NA),
      panel.background = ggplot2::element_rect(fill = "#F4F7F9", color = NA),
      panel.grid.major = ggplot2::element_line(color = "#FFFFFF", linewidth = 0.8),
      panel.grid.minor = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 18, color = "black"),
      plot.tag.position = c(0.015, 0.97),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13, color = "#111111", margin=ggplot2::margin(b=14)),
      axis.text = ggplot2::element_text(size = 11, color = "#222222", face = "bold"),
      axis.title = ggplot2::element_text(size = 12, color = "#111111", face = "bold"),

      legend.position = c(0.80, 0.72),
      legend.background = ggplot2::element_rect(fill = ggplot2::alpha("#F4F7F9", 0.85), color = NA),
      legend.key = ggplot2::element_blank(),
      legend.key.width = ggplot2::unit(1.5, "cm"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10, family = "serif", face = "bold"),
      plot.margin = ggplot2::margin(20, 20, 20, 20)
    )
}

p_dca_train <- plot_advanced_dca(dca_train, "Decision Curve (Training cohort)", "C")
p_dca_test  <- plot_advanced_dca(dca_test,  "Decision Curve (Testing cohort)",  "D")


fig_2x2 <- gridExtra::arrangeGrob(
  p_cal_train, p_cal_test, p_dca_train, p_dca_test,
  ncol = 2, padding = grid::unit(1.2, "line")
)


bg_rect <- grid::rectGrob(gp = grid::gpar(fill = "white", col = NA))
fig_2x2_final <- gridExtra::arrangeGrob(bg_rect, fig_2x2, heights = grid::unit(c(0, 1), c("null", "null")))

pdf(file.path(FIG_DIR, "Fig_Calibration_DCA_2x2.pdf"), width = 13.5, height = 11.5)
grid::grid.draw(fig_2x2_final)
dev.off()

png(file.path(FIG_DIR, "Fig_Calibration_DCA_2x2.png"), width = 3000, height = 2600, res = 240)
grid::grid.draw(fig_2x2_final)
dev.off()

cat("[OK] Advanced Calibration & DCA 2x2 Plots exported:\n",
    " - ", basename(file.path(FIG_DIR, "Fig_Calibration_DCA_2x2.pdf")), "\n", sep = "")


cat("\n=== Step 7-9.x: Unified diagnostic performance table ===\n")


.calc_cls_metrics <- function(y, score, thr){
  df <- data.frame(y = y, score = score)
  df <- df[!is.na(df$y) & is.finite(df$score), , drop = FALSE]
  if(nrow(df) == 0){
    return(list(n = 0, acc = NA, sen = NA, spe = NA, ppv = NA, npv = NA, f1 = NA))
  }
  pred <- ifelse(df$score >= thr, 1, 0)
  y0 <- ifelse(df$y == 1, 1, 0)
  TP <- sum(pred == 1 & y0 == 1)
  TN <- sum(pred == 0 & y0 == 0)
  FP <- sum(pred == 1 & y0 == 0)
  FN <- sum(pred == 0 & y0 == 1)

  acc <- (TP + TN) / (TP + TN + FP + FN)
  sen <- if((TP + FN) == 0) NA else TP / (TP + FN)
  spe <- if((TN + FP) == 0) NA else TN / (TN + FP)
  ppv <- if((TP + FP) == 0) NA else TP / (TP + FP)
  npv <- if((TN + FN) == 0) NA else TN / (TN + FN)
  f1  <- if(is.na(ppv) || is.na(sen) || (ppv + sen) == 0) NA else 2 * ppv * sen / (ppv + sen)

  list(n = nrow(df), acc = acc, sen = sen, spe = spe, ppv = ppv, npv = npv, f1 = f1)
}


.calc_auc_ci <- function(y, score){
  df <- data.frame(y = y, score = score)
  df <- df[!is.na(df$y) & is.finite(df$score), , drop = FALSE]
  if(nrow(df) < 3){
    return(list(auc = NA, ci = c(NA, NA, NA), roc = NULL, n = nrow(df)))
  }
  roc_obj <- pROC::roc(df$y, df$score, levels = c(0,1), direction = "<", quiet = TRUE)
  ci_obj  <- pROC::ci.auc(roc_obj)
  list(auc = as.numeric(pROC::auc(roc_obj)),
       ci  = as.numeric(ci_obj),
       roc = roc_obj,
       n   = nrow(df))
}


.get_youden_thr <- function(y, score){
  df <- data.frame(y = y, score = score)
  df <- df[!is.na(df$y) & is.finite(df$score), , drop = FALSE]
  if(nrow(df) < 5) return(NA_real_)
  roc_obj <- pROC::roc(df$y, df$score, levels = c(0,1), direction = "<", quiet = TRUE)
  as.numeric(pROC::coords(roc_obj, "best", ret = "threshold", best.method = "youden")[1])
}


.eval_model_train_test <- function(model_name, y_tr, s_tr, y_te, s_te, thr = NULL, thr_rule = "Youden(Train)"){
  if(is.null(thr)){
    thr <- .get_youden_thr(y_tr, s_tr)
    thr_rule <- "Youden(Train)"
  }

  a_tr <- .calc_auc_ci(y_tr, s_tr)
  a_te <- .calc_auc_ci(y_te, s_te)

  m_tr <- .calc_cls_metrics(y_tr, s_tr, thr)
  m_te <- .calc_cls_metrics(y_te, s_te, thr)

  sum_df <- rbind(
    data.frame(
      Model = model_name, Cohort = "Training",
      Threshold = thr, Threshold_Rule = thr_rule,
      AUC = a_tr$auc, Lower_95_CI = a_tr$ci[1], Upper_95_CI = a_tr$ci[3],
      Accuracy = m_tr$acc, Sensitivity = m_tr$sen, Specificity = m_tr$spe,
      PPV = m_tr$ppv, NPV = m_tr$npv, F1 = m_tr$f1,
      Sample_Size = m_tr$n,
      stringsAsFactors = FALSE
    ),
    data.frame(
      Model = model_name, Cohort = "Testing",
      Threshold = thr, Threshold_Rule = thr_rule,
      AUC = a_te$auc, Lower_95_CI = a_te$ci[1], Upper_95_CI = a_te$ci[3],
      Accuracy = m_te$acc, Sensitivity = m_te$sen, Specificity = m_te$spe,
      PPV = m_te$ppv, NPV = m_te$npv, F1 = m_te$f1,
      Sample_Size = m_te$n,
      stringsAsFactors = FALSE
    )
  )
  list(summary = sum_df, thr = thr)
}


fit_glm_predict <- function(train_df, test_df, vars){
  df_fit <- train_df[, c("group", vars), drop = FALSE]
  df_fit <- df_fit[stats::complete.cases(df_fit), , drop = FALSE]
  df_fit <- as.data.frame(df_fit)

  known_cat <- c("Sex","Disease_location","Disease_behavior","Smoking","pct_Delta_ITH_cat")
  for(v in intersect(vars, colnames(df_fit))){
    if(v %in% known_cat){
      df_fit[[v]] <- factor(trimws(as.character(df_fit[[v]])))
    }
  }

  align_to_fit <- function(df){
    out <- as.data.frame(df)
    for(v in intersect(vars, colnames(df_fit))){
      if(is.factor(df_fit[[v]])){
        out[[v]] <- factor(trimws(as.character(out[[v]])), levels = levels(df_fit[[v]]))
      }
    }
    out
  }

  tr <- align_to_fit(train_df)
  te <- align_to_fit(test_df)

  form <- stats::as.formula(paste("group ~", paste(vars, collapse = " + ")))
  fit <- stats::glm(form, data = df_fit, family = stats::binomial())

  tr$score <- as.numeric(stats::predict(fit, newdata = tr, type = "response"))
  te$score <- as.numeric(stats::predict(fit, newdata = te, type = "response"))
  list(tr = tr, te = te)
}


stopifnot(all(c("group","nom_prob") %in% colnames(clean_train)))
stopifnot(all(c("group","nom_prob") %in% colnames(clean_test)))

res_nom <- .eval_model_train_test(
  model_name = "Nomogram (Final lrm)",
  y_tr = clean_train$group, s_tr = clean_train$nom_prob,
  y_te = clean_test$group,  s_te = clean_test$nom_prob,
  thr = NULL
)


nom_auc_only <- res_nom$summary[, c("Model","Cohort","AUC","Lower_95_CI","Upper_95_CI","Sample_Size")]
write.csv(nom_auc_only,
          file.path(TABLE_DIR, "Nomogram_Performance_AUC.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("[OK] Nomogram AUC-only CSV saved.\n")


delta_train_df <- train_data[, c("group","pct_Delta_ITH_cat"), drop = FALSE]
delta_test_df  <- test_data[,  c("group","pct_Delta_ITH_cat"), drop = FALSE]
delta_train_df <- delta_train_df[stats::complete.cases(delta_train_df), , drop = FALSE]
delta_test_df  <- delta_test_df[stats::complete.cases(delta_test_df),  , drop = FALSE]

delta_train_score <- (1 - suppressWarnings(as.numeric(as.character(delta_train_df$pct_Delta_ITH_cat))))
delta_test_score  <- (1 - suppressWarnings(as.numeric(as.character(delta_test_df$pct_Delta_ITH_cat))))

res_delta <- .eval_model_train_test(
  model_name = "ΔITH-only (0/1)",
  y_tr = delta_train_df$group, s_tr = delta_train_score,
  y_te = delta_test_df$group,  s_te = delta_test_score,
  thr = 0.5, thr_rule = "Fixed(0.5)"
)

delta_auc_only <- res_delta$summary[, c("Model","Cohort","AUC","Lower_95_CI","Upper_95_CI","Sample_Size")]
write.csv(delta_auc_only,
          file.path(TABLE_DIR, "Delta_ITH_Binary_AUC.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
cat("[OK] Delta ITH-only AUC-only CSV saved.\n")


cmp_train <- train_data[!is.na(train_data$group), , drop = FALSE]
cmp_test  <- test_data[!is.na(test_data$group),  , drop = FALSE]
cmp_train <- as.data.frame(cmp_train)
cmp_test  <- as.data.frame(cmp_test)

perf_list <- list()
perf_list[[1]] <- res_nom$summary
perf_list[[2]] <- res_delta$summary


vars_clinical_only <- setdiff(sig_vars, c("pct_Delta_ITH_cat","pct_Delta_ITH"))
if(length(vars_clinical_only) >= 1){
  pr <- fit_glm_predict(cmp_train, cmp_test, vars_clinical_only)
  res_clin <- .eval_model_train_test(
    model_name = "Clinical-only (No ΔITH)",
    y_tr = pr$tr$group, s_tr = pr$tr$score,
    y_te = pr$te$group, s_te = pr$te$score,
    thr = NULL
  )
  perf_list[[length(perf_list)+1]] <- res_clin$summary
} else {
  cat("[WARN] Clinical-only vars empty, skip.\n")
}


pr_pre <- fit_glm_predict(cmp_train, cmp_test, c("Baseline_ITH"))
res_pre <- .eval_model_train_test(
  model_name = "Pre-ITH only",
  y_tr = pr_pre$tr$group, s_tr = pr_pre$tr$score,
  y_te = pr_pre$te$group, s_te = pr_pre$te$score,
  thr = NULL
)
perf_list[[length(perf_list)+1]] <- res_pre$summary


pr_post <- fit_glm_predict(cmp_train, cmp_test, c("Post_ITH"))
res_post <- .eval_model_train_test(
  model_name = "Post-ITH only",
  y_tr = pr_post$tr$group, s_tr = pr_post$tr$score,
  y_te = pr_post$te$group, s_te = pr_post$te$score,
  thr = NULL
)
perf_list[[length(perf_list)+1]] <- res_post$summary


RAD_XLSX <- file.path(getwd(), "Radiomics_Predictions.xlsx")
if(!file.exists(RAD_XLSX)){
  stop(paste0("[ERROR] Cannot find file: ", RAD_XLSX,
              "\nPlace Radiomics_Predictions.xlsx in the input directory."))
}

rad_df <- readxl::read_excel(RAD_XLSX)
rad_df <- as.data.frame(rad_df)

need_cols <- c("True_Group","Cohort","Radscore_Pre","Radscore_Post","Radscore_Delta")
miss_cols <- setdiff(need_cols, colnames(rad_df))
if(length(miss_cols) > 0){
  stop(paste0("[ERROR] Missing columns in Radiomics_Predictions.xlsx: ",
              paste(miss_cols, collapse=", ")))
}

rad_df$True_Group <- as.numeric(rad_df$True_Group)
rad_df$Cohort <- trimws(as.character(rad_df$Cohort))
rad_df$Cohort <- ifelse(tolower(rad_df$Cohort) %in% c("train","training"), "Training",
                        ifelse(tolower(rad_df$Cohort) %in% c("test","testing","validation","val"), "Testing", rad_df$Cohort))

rad_tr <- rad_df[rad_df$Cohort == "Training", , drop = FALSE]
rad_te <- rad_df[rad_df$Cohort == "Testing",  , drop = FALSE]

res_rad_pre <- .eval_model_train_test(
  model_name = "Radiomics (Pre) Radscore",
  y_tr = rad_tr$True_Group, s_tr = rad_tr$Radscore_Pre,
  y_te = rad_te$True_Group, s_te = rad_te$Radscore_Pre,
  thr = NULL
)
perf_list[[length(perf_list)+1]] <- res_rad_pre$summary

res_rad_post <- .eval_model_train_test(
  model_name = "Radiomics (Post) Radscore",
  y_tr = rad_tr$True_Group, s_tr = rad_tr$Radscore_Post,
  y_te = rad_te$True_Group, s_te = rad_te$Radscore_Post,
  thr = NULL
)
perf_list[[length(perf_list)+1]] <- res_rad_post$summary

res_rad_delta <- .eval_model_train_test(
  model_name = "Radiomics (Delta) Radscore",
  y_tr = rad_tr$True_Group, s_tr = rad_tr$Radscore_Delta,
  y_te = rad_te$True_Group, s_te = rad_te$Radscore_Delta,
  thr = NULL
)
perf_list[[length(perf_list)+1]] <- res_rad_delta$summary


all_perf <- do.call(rbind, perf_list)
all_perf <- as.data.frame(all_perf)


num_cols <- c("AUC","Lower_95_CI","Upper_95_CI",
              "Accuracy","Sensitivity","Specificity","PPV","NPV","F1",
              "Sample_Size","Threshold")
for(cc in intersect(num_cols, colnames(all_perf))){
  all_perf[[cc]] <- suppressWarnings(as.numeric(all_perf[[cc]]))
}

all_perf$AUC <- round(all_perf$AUC, 3)
all_perf$Lower_95_CI <- round(all_perf$Lower_95_CI, 3)
all_perf$Upper_95_CI <- round(all_perf$Upper_95_CI, 3)

all_perf$Accuracy <- round(all_perf$Accuracy, 3)
all_perf$Sensitivity <- round(all_perf$Sensitivity, 3)
all_perf$Specificity <- round(all_perf$Specificity, 3)
all_perf$PPV <- round(all_perf$PPV, 3)
all_perf$NPV <- round(all_perf$NPV, 3)
all_perf$F1 <- round(all_perf$F1, 3)

all_perf$Threshold <- round(all_perf$Threshold, 6)


keep_cols <- c("Model","Cohort",
               "AUC","Lower_95_CI","Upper_95_CI",
               "Accuracy","Sensitivity","Specificity","PPV","NPV","F1",
               "Sample_Size","Threshold","Threshold_Rule")
all_perf <- all_perf[, keep_cols, drop = FALSE]


out_csv <- file.path(TABLE_DIR, "Compare_Diagnostic_Performance_AllModels.csv")
write.csv(all_perf, out_csv, row.names = FALSE, fileEncoding = "UTF-8")
cat("[OK] Full diagnostic performance table saved:\n", basename(out_csv), "\n")


thr_df <- all_perf[all_perf$Cohort == "Training", c("Model","Threshold","Threshold_Rule"), drop = FALSE]
thr_df <- unique(thr_df)

thr_csv <- file.path(TABLE_DIR, "Model_Thresholds_TrainYouden.csv")
write.csv(thr_df, thr_csv, row.names = FALSE, fileEncoding = "UTF-8")
cat("[OK] Threshold summary saved:\n", basename(thr_csv), "\n")

print(all_perf)


cat("\n=== Step 14: Color ROC (Nom-%ΔVHS vs %ΔVHS vs Clin_pre, NPJ Style) ===\n")


stopifnot(all(c("group","nom_prob") %in% colnames(clean_train)))
stopifnot(all(c("group","nom_prob") %in% colnames(clean_test)))
stopifnot("pct_Delta_ITH_cat" %in% colnames(clean_train))
stopifnot("pct_Delta_ITH_cat" %in% colnames(clean_test))


roc_nom_tr <- pROC::roc(clean_train$group, clean_train$nom_prob, levels = c(0,1), direction = "<", quiet = TRUE)
roc_nom_te <- pROC::roc(clean_test$group, clean_test$nom_prob, levels = c(0,1), direction = "<", quiet = TRUE)


score_dith_tr <- (1 - suppressWarnings(as.numeric(as.character(clean_train$pct_Delta_ITH_cat))))
score_dith_te <- (1 - suppressWarnings(as.numeric(as.character(clean_test$pct_Delta_ITH_cat))))

roc_dith_tr <- pROC::roc(clean_train$group, score_dith_tr, levels = c(0,1), direction = "<", quiet = TRUE)
roc_dith_te <- pROC::roc(clean_test$group,  score_dith_te, levels = c(0,1), direction = "<", quiet = TRUE)


VHS_RELATED <- c("pct_Delta_VHS","pct_Delta_VHS_cat","Delta_VHS","VHS_pre","VHS_post",
                 "pct_Delta_ITH","pct_Delta_ITH_cat","Baseline_ITH","Post_ITH")
vars_clinical_only <- setdiff(sig_vars, intersect(sig_vars, VHS_RELATED))

fit_glm_predict_one <- function(train_df, test_df, vars){
  df_fit <- train_df[, c("group", vars), drop = FALSE]
  df_fit <- df_fit[stats::complete.cases(df_fit), , drop = FALSE]
  df_fit <- as.data.frame(df_fit)

  known_cat <- c("Sex","Disease_location","Disease_behavior","Smoking","pct_Delta_ITH_cat","pct_Delta_VHS_cat")
  for(v in intersect(vars, colnames(df_fit))){
    if(v %in% known_cat) df_fit[[v]] <- factor(trimws(as.character(df_fit[[v]])))
  }

  align_to_fit <- function(df){
    out <- as.data.frame(df)
    for(v in intersect(vars, colnames(df_fit))){
      if(is.factor(df_fit[[v]])) out[[v]] <- factor(trimws(as.character(out[[v]])), levels = levels(df_fit[[v]]))
    }
    out
  }

  tr <- align_to_fit(train_df); te <- align_to_fit(test_df)
  form <- stats::as.formula(paste("group ~", paste(vars, collapse = " + ")))
  fit <- stats::glm(form, data = df_fit, family = stats::binomial())

  tr$score <- as.numeric(stats::predict(fit, newdata = tr, type = "response"))
  te$score <- as.numeric(stats::predict(fit, newdata = te, type = "response"))
  list(tr = tr, te = te, vars = vars)
}

cmp_train <- as.data.frame(train_data[!is.na(train_data$group), , drop = FALSE])
cmp_test  <- as.data.frame(test_data[!is.na(test_data$group),  , drop = FALSE])

if(length(vars_clinical_only) < 1) stop("[ERROR] Clinical-only variable set is empty.")

clin_pred <- fit_glm_predict_one(cmp_train, cmp_test, vars_clinical_only)
clin_tr_clean <- clin_pred$tr[!is.na(clin_pred$tr$group) & is.finite(clin_pred$tr$score), , drop = FALSE]
clin_te_clean <- clin_pred$te[!is.na(clin_pred$te$group) & is.finite(clin_pred$te$score), , drop = FALSE]

roc_clin_tr <- pROC::roc(clin_tr_clean$group, clin_tr_clean$score, levels = c(0,1), direction = "<", quiet = TRUE)
roc_clin_te <- pROC::roc(clin_te_clean$group, clin_te_clean$score, levels = c(0,1), direction = "<", quiet = TRUE)

auc_nom_tr <- as.numeric(pROC::auc(roc_nom_tr)); auc_nom_te <- as.numeric(pROC::auc(roc_nom_te))
auc_dit_tr <- as.numeric(pROC::auc(roc_dith_tr)); auc_dit_te <- as.numeric(pROC::auc(roc_dith_te))
auc_cln_tr <- as.numeric(pROC::auc(roc_clin_tr)); auc_cln_te <- as.numeric(pROC::auc(roc_clin_te))


extract_roc_df <- function(r, model_name) {
  data.frame(FPR = 1 - r$specificities, TPR = r$sensitivities, Model = model_name)
}

labels_tr <- c(sprintf("Nom-%%ΔVHS (AUC=%.3f)", auc_nom_tr),
               sprintf("%%ΔVHS (AUC=%.3f)", auc_dit_tr),
               sprintf("Clin_pre (AUC=%.3f)", auc_cln_tr))

labels_te <- c(sprintf("Nom-%%ΔVHS (AUC=%.3f)", auc_nom_te),
               sprintf("%%ΔVHS (AUC=%.3f)", auc_dit_te),
               sprintf("Clin_pre (AUC=%.3f)", auc_cln_te))

df_roc_tr <- rbind(extract_roc_df(roc_nom_tr, labels_tr[1]),
                   extract_roc_df(roc_dith_tr, labels_tr[2]),
                   extract_roc_df(roc_clin_tr, labels_tr[3]))
df_roc_tr$Model <- factor(df_roc_tr$Model, levels = labels_tr)

df_roc_te <- rbind(extract_roc_df(roc_nom_te, labels_te[1]),
                   extract_roc_df(roc_dith_te, labels_te[2]),
                   extract_roc_df(roc_clin_te, labels_te[3]))
df_roc_te$Model <- factor(df_roc_te$Model, levels = labels_te)


plot_npj_roc <- function(df, title, tag) {
  ggplot2::ggplot(df, ggplot2::aes(x = FPR, y = TPR, color = Model, linetype = Model)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "#7B8B9A", linewidth = 0.8) +
    ggplot2::labs(title = title, x = "1 - Specificity", y = "Sensitivity", tag = tag) +
    ggplot2::scale_color_manual(values = c("#3C5488", "#E64B35", "#00A087")) +
    ggplot2::scale_linetype_manual(values = c("solid", "longdash", "dotdash")) +
    ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "#F4F7F9", color = NA),
      panel.background = ggplot2::element_rect(fill = "#F4F7F9", color = NA),
      panel.grid.major = ggplot2::element_line(color = "#FFFFFF", linewidth = 0.8),
      panel.grid.minor = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 18, color = "black"),
      plot.tag.position = c(0.015, 0.97),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13, color = "#111111", margin=ggplot2::margin(b=14)),
      axis.text = ggplot2::element_text(size = 11, color = "#222222", face = "bold"),
      axis.title = ggplot2::element_text(size = 12, color = "#111111", face = "bold"),
      legend.position = c(0.70, 0.22),
      legend.background = ggplot2::element_rect(fill = ggplot2::alpha("#F4F7F9", 0.7), color = NA),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10, color = "#222222", family="serif")
    ) +
    ggplot2::coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
    ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 20, 20))
}

p_roc_tr <- plot_npj_roc(df_roc_tr, "ROC: Training cohort", "A")
p_roc_te <- plot_npj_roc(df_roc_te, "ROC: Testing cohort", "B")

fig_roc_1x2 <- gridExtra::arrangeGrob(p_roc_tr, p_roc_te, ncol=2, padding=grid::unit(1.2, "line"))
bg_rect <- grid::rectGrob(gp = grid::gpar(fill = "white", col = NA))
fig_roc_final <- gridExtra::arrangeGrob(bg_rect, fig_roc_1x2, heights = grid::unit(c(0, 1), c("null", "null")))

pdf(file.path(FIG_DIR, "ROC_Color_NPJ_Style.pdf"), width = 11.5, height = 5.5)
grid::grid.draw(fig_roc_final)
dev.off()

png(file.path(FIG_DIR, "ROC_Color_NPJ_Style.png"), width = 2600, height = 1200, res = 220)
grid::grid.draw(fig_roc_final)
dev.off()

cat("[OK] NPJ Style ROC exported:\n")


cat("\n=== Step 15: Export professional Excel ===\n")

if(!exists("all_perf")){
  stop("[ERROR] all_perf not found. Please run Step 7-9.x unified performance block first.")
}

perf_prof <- as.data.frame(all_perf)


perf_prof$Model <- trimws(as.character(perf_prof$Model))
perf_prof$Model <- dplyr::recode(
  perf_prof$Model,
  "Nomogram (Final lrm)"       = "Nom-%ΔVHS",
  "Clinical-only (No ΔITH)"    = "Clin_pre",
  "ΔITH-only (0/1)"            = "%ΔVHS",
  "Pre-ITH only"               = "VHS_pre",
  "Post-ITH only"              = "VHS_post",
  "Radiomics (Pre) Radscore"   = "VAT-RadScore_pre",
  "Radiomics (Post) Radscore"  = "VAT-RadScore_post",
  "Radiomics (Delta) Radscore" = "ΔVAT-RadScore"
)

perf_prof$AUC_95CI <- sprintf("%.3f (%.3f–%.3f)",
                              as.numeric(perf_prof$AUC),
                              as.numeric(perf_prof$Lower_95_CI),
                              as.numeric(perf_prof$Upper_95_CI))

keep_excel_cols <- c("Model","Cohort","AUC_95CI",
                     "Accuracy","Sensitivity","Specificity","PPV","NPV","F1",
                     "Sample_Size","Threshold","Threshold_Rule")
keep_excel_cols <- keep_excel_cols[keep_excel_cols %in% colnames(perf_prof)]
perf_prof <- perf_prof[, keep_excel_cols, drop = FALSE]

thr_prof <- perf_prof[perf_prof$Cohort == "Training", c("Model","Threshold","Threshold_Rule"), drop = FALSE]
thr_prof <- unique(thr_prof)

model_defs <- data.frame(
  Term = c("Nom-%ΔVHS", "%ΔVHS", "Clin_pre", "VHS_pre", "VHS_post",
           "VAT-RadScore_pre", "VAT-RadScore_post", "ΔVAT-RadScore"),
  Full_Name_EN = c(
    "Nomogram integrating %ΔVHS and pre-treatment clinical/imaging factors",
    "Percentage change in Visceral Adipose Tissue Heterogeneity Score",
    "Pre-treatment clinical/imaging factors",
    "Pre-treatment Visceral Adipose Tissue Heterogeneity Score",
    "Post-treatment Visceral Adipose Tissue Heterogeneity Score",
    "Pre-treatment Visceral Adipose Tissue Radiomics Score",
    "Post-treatment Visceral Adipose Tissue Radiomics Score",
    "Delta Visceral Adipose Tissue Radiomics Score"
  ),
  Definition_EN = c(
    "Primary nomogram combining %DeltaVHS with pretreatment clinical and imaging predictors",
    "Percentage change in VAT heterogeneity score: (VHS_pre - VHS_post) / VHS_pre x 100%",
    "Clinical-only model (pretreatment predictors without %DeltaVHS)",
    "Pretreatment VAT heterogeneity score",
    "Post-induction VAT heterogeneity score",
    "Pretreatment VAT radiomics score",
    "Post-induction VAT radiomics score",
    "Change in VAT radiomics score (VAT-RadScore_post - VAT-RadScore_pre)"
  ),
  stringsAsFactors = FALSE
)

if(!requireNamespace("writexl", quietly = TRUE)) stop("Missing required R package: writexl")
writexl::write_xlsx(
  list("Diagnostic_Performance" = perf_prof,
       "Thresholds_Training"    = thr_prof,
       "Terminology"            = model_defs),
  path = file.path(TABLE_DIR, "Diagnostic_Performance_Professional.xlsx")
)

cat("[OK] Excel exported: Diagnostic_Performance_Professional.xlsx\n")


cat("\n=== Step 16: Visualize diagnostic metrics (Nom-%ΔVHS baseline heatmap) ===\n")

plot_num <- as.data.frame(all_perf)
plot_num$Cohort <- factor(trimws(as.character(plot_num$Cohort)), levels = c("Training","Testing"))

plot_num$Model <- trimws(as.character(plot_num$Model))
plot_num$Model <- dplyr::recode(
  plot_num$Model,
  "Nomogram (Final lrm)"       = "Nom-%ΔVHS",
  "Clinical-only (No ΔITH)"    = "Clin_pre",
  "ΔITH-only (0/1)"            = "%ΔVHS",
  "Pre-ITH only"               = "VHS_pre",
  "Post-ITH only"              = "VHS_post",
  "Radiomics (Pre) Radscore"   = "VAT-RadScore_pre",
  "Radiomics (Post) Radscore"  = "VAT-RadScore_post",
  "Radiomics (Delta) Radscore" = "ΔVAT-RadScore"
)

need_num_metrics <- c("AUC","Accuracy","Sensitivity","Specificity","PPV","NPV","F1")
for(m in need_num_metrics) plot_num[[m]] <- suppressWarnings(as.numeric(plot_num[[m]]))

long_list <- vector("list", length = nrow(plot_num) * length(need_num_metrics))
k <- 1
for(i in seq_len(nrow(plot_num))){
  for(m in need_num_metrics){
    long_list[[k]] <- data.frame(Model=plot_num$Model[i], Cohort=plot_num$Cohort[i], Metric=m, Raw=plot_num[[m]][i], stringsAsFactors=FALSE)
    k <- k + 1
  }
}
plot_long <- do.call(rbind, long_list)
plot_long <- plot_long[is.finite(plot_long$Raw) & !is.na(plot_long$Raw) & !is.na(plot_long$Cohort), , drop = FALSE]

BASE_NAME <- "Nom-%ΔVHS"
base_df <- plot_long[plot_long$Model == BASE_NAME, c("Cohort","Metric","Raw"), drop = FALSE]
colnames(base_df)[colnames(base_df)=="Raw"] <- "BaseRaw"

plot_long <- merge(plot_long, base_df, by = c("Cohort","Metric"), all.x = TRUE)
plot_long$Delta <- plot_long$Raw - plot_long$BaseRaw

plot_long$Label <- ifelse(
  plot_long$Model == BASE_NAME,
  sprintf("0.000\n(%.3f)", plot_long$Raw),
  sprintf("%+.3f\n(%.3f)", plot_long$Delta, plot_long$Raw)
)

plot_long$Metric <- factor(plot_long$Metric, levels = c("AUC","Accuracy","Sensitivity","Specificity","PPV","NPV","F1"))
model_levels <- c(BASE_NAME, setdiff(unique(as.character(plot_long$Model)), BASE_NAME))
plot_long$Model <- factor(plot_long$Model, levels = rev(model_levels))

max_abs <- max(abs(plot_long$Delta), na.rm = TRUE)
if(!is.finite(max_abs) || max_abs == 0) max_abs <- 0.001

p_heat <- ggplot2::ggplot(plot_long, ggplot2::aes(x = Metric, y = Model, fill = Delta)) +
  ggplot2::geom_tile(color = "white", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = Label), size = 3.0, family = "serif", lineheight = 0.9) +
  ggplot2::facet_wrap(~ Cohort, ncol = 2) +
  ggplot2::scale_fill_gradient2(limits = c(-max_abs, max_abs), midpoint = 0, low="#3C5488", high="#E64B35") +
  ggplot2::labs(title = "Diagnostic Performance vs Nom-%ΔVHS (Δ = Model − Nom-%ΔVHS)", x = NULL, y = NULL, fill = "ΔMetric") +
  ggplot2::theme_bw(base_family = "serif") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5, face="bold"),
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5),
    legend.position = "right",
    strip.background = ggplot2::element_rect(fill = "#F4F7F9", color = NA),
    strip.text = ggplot2::element_text(face="bold", size=11)
  )

pdf(file.path(FIG_DIR, "Heatmap_Diagnostic_Performance.pdf"), width = 13, height = 7.8)
print(p_heat)
dev.off()

png(file.path(FIG_DIR, "Heatmap_Diagnostic_Performance.png"), width = 2800, height = 1700, res = 220)
print(p_heat)
dev.off()

cat("[OK] Heatmap exported.\n")


cat("\n=== Step 17: Calculate and Visualize NRI & IDI (Baseline = Nom-%ΔVHS) ===\n")


calc_nri_idi <- function(p.std, p.new, y) {
  idx <- is.finite(p.std) & is.finite(p.new) & is.finite(y)
  p.std <- p.std[idx]; p.new <- p.new[idx]; y <- as.numeric(y[idx])

  if(length(y) < 2) return(data.frame(IDI=NA, IDI_lower=NA, IDI_upper=NA, IDI_p=NA, NRI=NA, NRI_lower=NA, NRI_upper=NA, NRI_p=NA))

  p.std.1 <- p.std[y == 1]; p.std.0 <- p.std[y == 0]
  p.new.1 <- p.new[y == 1]; p.new.0 <- p.new[y == 0]

  if(length(p.std.1) < 2 || length(p.std.0) < 2)
    return(data.frame(IDI=NA, IDI_lower=NA, IDI_upper=NA, IDI_p=NA, NRI=NA, NRI_lower=NA, NRI_upper=NA, NRI_p=NA))

  idi.event <- mean(p.new.1) - mean(p.std.1)
  idi.nonevent <- mean(p.new.0) - mean(p.std.0)
  idi <- idi.event - idi.nonevent

  diff.1 <- p.new.1 - p.std.1; diff.0 <- p.new.0 - p.std.0
  se.idi <- sqrt(var(diff.1)/length(diff.1) + var(diff.0)/length(diff.0))
  z.idi <- idi / se.idi
  p.idi <- 2 * (1 - pnorm(abs(z.idi)))

  up.event <- sum(p.new.1 > p.std.1); down.event <- sum(p.new.1 < p.std.1)
  up.nonevent <- sum(p.new.0 > p.std.0); down.nonevent <- sum(p.new.0 < p.std.0)

  nri.event <- (up.event - down.event) / length(p.new.1)
  nri.nonevent <- (down.nonevent - up.nonevent) / length(p.new.0)
  nri <- nri.event + nri.nonevent

  se.nri <- sqrt((1 - nri.event^2)/length(p.new.1) + (1 - nri.nonevent^2)/length(p.new.0))
  z.nri <- nri / se.nri
  p.nri <- 2 * (1 - pnorm(abs(z.nri)))

  data.frame(IDI = idi, IDI_lower = idi - 1.96*se.idi, IDI_upper = idi + 1.96*se.idi, IDI_p = p.idi,
             NRI = nri, NRI_lower = nri - 1.96*se.nri, NRI_upper = nri + 1.96*se.nri, NRI_p = p.nri)
}


get_record_key <- function(df) {
  keys <- c("ID", "id", "CaseID", "case_id", "PatientID", "patient_id", "SubjectID", "subject_id")
  found <- intersect(keys, names(df))
  if (length(found) == 0L) stop("A shared record identifier column is required to match input tables")
  values <- as.character(df[[found[1]]])
  if (anyNA(values) || anyDuplicated(values)) stop("Record identifiers must be unique and nonmissing")
  values
}
train_data$ID <- get_record_key(train_data)
test_data$ID <- get_record_key(test_data)
rad_df$ID <- get_record_key(rad_df)

train_merged <- merge(train_data, rad_df[, c("ID", "Radscore_Pre", "Radscore_Post", "Radscore_Delta")], by="ID", all.x=TRUE)
test_merged  <- merge(test_data,  rad_df[, c("ID", "Radscore_Pre", "Radscore_Post", "Radscore_Delta")], by="ID", all.x=TRUE)

train_merged$y_num <- as.numeric(as.character(train_merged$group))
test_merged$y_num  <- as.numeric(as.character(test_merged$group))


train_comp_base <- train_merged[complete.cases(train_merged[, c("y_num", sig_vars)]), , drop=FALSE]
form_base <- as.formula(paste("y_num ~", paste(sig_vars, collapse=" + ")))
fit_base <- glm(form_base, data = train_comp_base, family = binomial)

prob_base_tr <- suppressWarnings(as.numeric(predict(fit_base, newdata = train_merged, type = "response")))
prob_base_te <- suppressWarnings(as.numeric(predict(fit_base, newdata = test_merged,  type = "response")))


model_configs <- list(
  "ΔITH-only (0/1)"            = c("pct_Delta_ITH_cat"),
  "Clinical-only (No ΔITH)"    = vars_clinical_only,
  "Pre-ITH only"               = c("Baseline_ITH"),
  "Post-ITH only"              = c("Post_ITH"),
  "Radiomics (Pre) Radscore"   = c("Radscore_Pre"),
  "Radiomics (Post) Radscore"  = c("Radscore_Post"),
  "Radiomics (Delta) Radscore" = c("Radscore_Delta")
)

nri_idi_results <- list()

for(m_name in names(model_configs)) {
  v <- model_configs[[m_name]]
  train_comp <- train_merged[complete.cases(train_merged[, c("y_num", v)]), , drop=FALSE]
  if(nrow(train_comp) == 0) next

  form <- as.formula(paste("y_num ~", paste(v, collapse=" + ")))
  fit_m <- glm(form, data = train_comp, family = binomial)

  p_m_tr <- suppressWarnings(as.numeric(predict(fit_m, newdata = train_merged, type = "response")))
  p_m_te <- suppressWarnings(as.numeric(predict(fit_m, newdata = test_merged,  type = "response")))

  res_tr <- calc_nri_idi(prob_base_tr, p_m_tr, train_merged$y_num)
  res_te <- calc_nri_idi(prob_base_te, p_m_te, test_merged$y_num)
  res_tr$Model <- m_name; res_tr$Cohort <- "Training"
  res_te$Model <- m_name; res_te$Cohort <- "Testing"

  nri_idi_results[[length(nri_idi_results) + 1]] <- res_tr
  nri_idi_results[[length(nri_idi_results) + 1]] <- res_te
}
df_nri_idi <- do.call(rbind, nri_idi_results)


df_nri_idi$Model_Prof <- dplyr::recode(df_nri_idi$Model,
                                       "ΔITH-only (0/1)"            = "%ΔVHS",
                                       "Clinical-only (No ΔITH)"    = "Clin_pre",
                                       "Pre-ITH only"               = "VHS_pre",
                                       "Post-ITH only"              = "VHS_post",
                                       "Radiomics (Pre) Radscore"   = "VAT-RadScore_pre",
                                       "Radiomics (Post) Radscore"  = "VAT-RadScore_post",
                                       "Radiomics (Delta) Radscore" = "ΔVAT-RadScore"
)

level_order <- rev(c("%ΔVHS", "Clin_pre", "VHS_pre", "VHS_post", "VAT-RadScore_pre", "VAT-RadScore_post", "ΔVAT-RadScore"))
df_nri_idi$Model_Prof <- factor(df_nri_idi$Model_Prof, levels = level_order)

write.csv(df_nri_idi, file.path(TABLE_DIR, "NRI_IDI_Results_vs_Nomogram.csv"), row.names = FALSE)
cat("[OK] NRI/IDI Table saved.\n")


plot_nri_idi_forest <- function(df, metric_col, lower_col, upper_col, p_col, title, x_label, plot_tag,
                                bg_color = "#F4F7F9", line_color = "#3C5488", point_color = "#E64B35") {
  df$p_lbl <- ifelse(is.na(df[[p_col]]), "NA", ifelse(df[[p_col]] < 0.001, "<0.001", sprintf("=%.3f", df[[p_col]])))
  df$lbl <- sprintf("%.3f (p%s)", df[[metric_col]], df$p_lbl)

  ggplot2::ggplot(df, ggplot2::aes(x = .data[[metric_col]], y = Model_Prof)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "#7B8B9A", linewidth = 0.8) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data[[lower_col]], xmax = .data[[upper_col]]), height = 0.15, color = line_color, linewidth = 0.9) +
    ggplot2::geom_point(color = point_color, size = 3.5) +
    ggplot2::geom_text(ggplot2::aes(label = lbl), vjust = -1.1, size = 3.6, family = "serif", color = "#222222") +
    ggplot2::labs(title = title, x = x_label, y = NULL, tag = plot_tag) +
    ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.grid.major.y = ggplot2::element_line(color = "#FFFFFF", linewidth = 0.8),
      panel.grid.major.x = ggplot2::element_line(color = "#FFFFFF", linewidth = 0.8),
      panel.grid.minor = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 18, color = "black"),
      plot.tag.position = c(0.015, 0.97),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13, color = "#111111", margin=ggplot2::margin(b=14)),
      axis.text.y = ggplot2::element_text(size = 11, color = "#222222", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10, color = "#222222"),
      axis.title.x = ggplot2::element_text(size = 11, color = "#111111", face = "bold", margin=ggplot2::margin(t=10)),
      plot.margin = ggplot2::margin(20, 20, 15, 10)
    ) + ggplot2::coord_cartesian(clip = "off")
}

df_tr <- df_nri_idi[df_nri_idi$Cohort == "Training", ]; df_te <- df_nri_idi[df_nri_idi$Cohort == "Testing", ]

p_nri_tr <- plot_nri_idi_forest(df_tr, "NRI", "NRI_lower", "NRI_upper", "NRI_p", "Continuous NRI (Training cohort)", "Continuous NRI (vs. Nom-%ΔVHS)", "A", "#F4F7F9")
p_nri_te <- plot_nri_idi_forest(df_te, "NRI", "NRI_lower", "NRI_upper", "NRI_p", "Continuous NRI (Testing cohort)",  "Continuous NRI (vs. Nom-%ΔVHS)", "B", "#F4F7F9")
p_idi_tr <- plot_nri_idi_forest(df_tr, "IDI", "IDI_lower", "IDI_upper", "IDI_p", "IDI (Training cohort)", "IDI (vs. Nom-%ΔVHS)", "C", "#FCF8F5")
p_idi_te <- plot_nri_idi_forest(df_te, "IDI", "IDI_lower", "IDI_upper", "IDI_p", "IDI (Testing cohort)",  "IDI (vs. Nom-%ΔVHS)", "D", "#FCF8F5")

fig_nri_idi_2x2 <- gridExtra::arrangeGrob(p_nri_tr, p_nri_te, p_idi_tr, p_idi_te, ncol = 2, padding = grid::unit(1.2, "line"))
bg_rect <- grid::rectGrob(gp = grid::gpar(fill = "white", col = NA))
fig_nri_idi_2x2_final <- gridExtra::arrangeGrob(bg_rect, fig_nri_idi_2x2, heights = grid::unit(c(0, 1), c("null", "null")))

pdf(file.path(FIG_DIR, "Fig_NRI_IDI_vs_Nomogram_2x2.pdf"), width = 14.5, height = 11)
grid::grid.draw(fig_nri_idi_2x2_final)
dev.off()

png(file.path(FIG_DIR, "Fig_NRI_IDI_vs_Nomogram_2x2.png"), width = 3400, height = 2600, res = 240)
grid::grid.draw(fig_nri_idi_2x2_final)
dev.off()

cat("[OK] NRI/IDI 2x2 Plots (Lancet + npj style) exported:\n",
    " - ", basename(file.path(FIG_DIR, "Fig_NRI_IDI_vs_Nomogram_2x2.pdf")), "\n",
    " - ", basename(file.path(FIG_DIR, "Fig_NRI_IDI_vs_Nomogram_2x2.png")), "\n", sep = "")
cat("\n=== ALL PIPELINE STEPS COMPLETED SUCCESSFULLY ===\n")


cat("\n=== Step 16 & 17: Build Mega-Figure (Heatmap, NRI, IDI in 2x3 Grid) ===\n")


stopifnot(exists("all_perf"))
plot_num <- as.data.frame(all_perf)
plot_num$Cohort <- factor(trimws(as.character(plot_num$Cohort)), levels = c("Training","Testing"))

plot_num$Model <- trimws(as.character(plot_num$Model))
plot_num$Model <- dplyr::recode(
  plot_num$Model,
  "Nomogram (Final lrm)"       = "Nom-%ΔVHS",
  "Clinical-only (No ΔITH)"    = "Clin_pre",
  "ΔITH-only (0/1)"            = "%ΔVHS",
  "Pre-ITH only"               = "VHS_pre",
  "Post-ITH only"              = "VHS_post",
  "Radiomics (Pre) Radscore"   = "VAT-RadScore_pre",
  "Radiomics (Post) Radscore"  = "VAT-RadScore_post",
  "Radiomics (Delta) Radscore" = "ΔVAT-RadScore"
)

need_num_metrics <- c("AUC","Accuracy","Sensitivity","Specificity","PPV","NPV","F1")
for(m in need_num_metrics) plot_num[[m]] <- suppressWarnings(as.numeric(plot_num[[m]]))

long_list <- vector("list", length = nrow(plot_num) * length(need_num_metrics))
k <- 1
for(i in seq_len(nrow(plot_num))){
  for(m in need_num_metrics){
    long_list[[k]] <- data.frame(Model=plot_num$Model[i], Cohort=plot_num$Cohort[i], Metric=m, Raw=plot_num[[m]][i], stringsAsFactors=FALSE)
    k <- k + 1
  }
}
plot_long <- do.call(rbind, long_list)
plot_long <- plot_long[is.finite(plot_long$Raw) & !is.na(plot_long$Raw) & !is.na(plot_long$Cohort), , drop = FALSE]

BASE_NAME <- "Nom-%ΔVHS"
base_df <- plot_long[plot_long$Model == BASE_NAME, c("Cohort","Metric","Raw"), drop = FALSE]
colnames(base_df)[colnames(base_df)=="Raw"] <- "BaseRaw"

plot_long <- merge(plot_long, base_df, by = c("Cohort","Metric"), all.x = TRUE)
plot_long$Delta <- plot_long$Raw - plot_long$BaseRaw

plot_long$Label <- ifelse(
  plot_long$Model == BASE_NAME,
  sprintf("0.000\n(%.3f)", plot_long$Raw),
  sprintf("%+.3f\n(%.3f)", plot_long$Delta, plot_long$Raw)
)

plot_long$Metric <- factor(plot_long$Metric, levels = c("AUC","Accuracy","Sensitivity","Specificity","PPV","NPV","F1"))
model_levels <- c(BASE_NAME, setdiff(unique(as.character(plot_long$Model)), BASE_NAME))
plot_long$Model <- factor(plot_long$Model, levels = rev(model_levels))

max_abs_heat <- max(abs(plot_long$Delta), na.rm = TRUE)
if(!is.finite(max_abs_heat) || max_abs_heat == 0) max_abs_heat <- 0.001

df_heat_tr <- plot_long[plot_long$Cohort == "Training", ]
df_heat_te <- plot_long[plot_long$Cohort == "Testing", ]


calc_nri_idi <- function(p.std, p.new, y) {
  idx <- is.finite(p.std) & is.finite(p.new) & is.finite(y)
  p.std <- p.std[idx]; p.new <- p.new[idx]; y <- as.numeric(y[idx])

  if(length(y) < 2) return(data.frame(IDI=NA, IDI_lower=NA, IDI_upper=NA, IDI_p=NA, NRI=NA, NRI_lower=NA, NRI_upper=NA, NRI_p=NA))

  p.std.1 <- p.std[y == 1]; p.std.0 <- p.std[y == 0]
  p.new.1 <- p.new[y == 1]; p.new.0 <- p.new[y == 0]

  if(length(p.std.1) < 2 || length(p.std.0) < 2)
    return(data.frame(IDI=NA, IDI_lower=NA, IDI_upper=NA, IDI_p=NA, NRI=NA, NRI_lower=NA, NRI_upper=NA, NRI_p=NA))

  idi.event <- mean(p.new.1) - mean(p.std.1)
  idi.nonevent <- mean(p.new.0) - mean(p.std.0)
  idi <- idi.event - idi.nonevent

  diff.1 <- p.new.1 - p.std.1; diff.0 <- p.new.0 - p.std.0
  se.idi <- sqrt(var(diff.1)/length(diff.1) + var(diff.0)/length(diff.0))
  z.idi <- idi / se.idi; p.idi <- 2 * (1 - pnorm(abs(z.idi)))

  up.event <- sum(p.new.1 > p.std.1); down.event <- sum(p.new.1 < p.std.1)
  up.nonevent <- sum(p.new.0 > p.std.0); down.nonevent <- sum(p.new.0 < p.std.0)

  nri.event <- (up.event - down.event) / length(p.new.1)
  nri.nonevent <- (down.nonevent - up.nonevent) / length(p.new.0)
  nri <- nri.event + nri.nonevent

  se.nri <- sqrt((1 - nri.event^2)/length(p.new.1) + (1 - nri.nonevent^2)/length(p.new.0))
  z.nri <- nri / se.nri; p.nri <- 2 * (1 - pnorm(abs(z.nri)))

  data.frame(IDI = idi, IDI_lower = idi - 1.96*se.idi, IDI_upper = idi + 1.96*se.idi, IDI_p = p.idi,
             NRI = nri, NRI_lower = nri - 1.96*se.nri, NRI_upper = nri + 1.96*se.nri, NRI_p = p.nri)
}

train_data$ID <- get_record_key(train_data)
test_data$ID <- get_record_key(test_data)
rad_df$ID <- get_record_key(rad_df)

train_merged <- merge(train_data, rad_df[, c("ID", "Radscore_Pre", "Radscore_Post", "Radscore_Delta")], by="ID", all.x=TRUE)
test_merged  <- merge(test_data,  rad_df[, c("ID", "Radscore_Pre", "Radscore_Post", "Radscore_Delta")], by="ID", all.x=TRUE)
train_merged$y_num <- as.numeric(as.character(train_merged$group))
test_merged$y_num  <- as.numeric(as.character(test_merged$group))

train_comp_base <- train_merged[complete.cases(train_merged[, c("y_num", sig_vars)]), , drop=FALSE]
fit_base <- glm(as.formula(paste("y_num ~", paste(sig_vars, collapse=" + "))), data = train_comp_base, family = binomial)
prob_base_tr <- suppressWarnings(as.numeric(predict(fit_base, newdata = train_merged, type = "response")))
prob_base_te <- suppressWarnings(as.numeric(predict(fit_base, newdata = test_merged,  type = "response")))

model_configs <- list(
  "ΔITH-only (0/1)"            = c("pct_Delta_ITH_cat"),
  "Clinical-only (No ΔITH)"    = vars_clinical_only,
  "Pre-ITH only"               = c("Baseline_ITH"),
  "Post-ITH only"              = c("Post_ITH"),
  "Radiomics (Pre) Radscore"   = c("Radscore_Pre"),
  "Radiomics (Post) Radscore"  = c("Radscore_Post"),
  "Radiomics (Delta) Radscore" = c("Radscore_Delta")
)

nri_idi_results <- list()
for(m_name in names(model_configs)) {
  v <- model_configs[[m_name]]
  train_comp <- train_merged[complete.cases(train_merged[, c("y_num", v)]), , drop=FALSE]
  if(nrow(train_comp) == 0) next

  fit_m <- glm(as.formula(paste("y_num ~", paste(v, collapse=" + "))), data = train_comp, family = binomial)
  p_m_tr <- suppressWarnings(as.numeric(predict(fit_m, newdata = train_merged, type = "response")))
  p_m_te <- suppressWarnings(as.numeric(predict(fit_m, newdata = test_merged,  type = "response")))

  res_tr <- calc_nri_idi(prob_base_tr, p_m_tr, train_merged$y_num)
  res_te <- calc_nri_idi(prob_base_te, p_m_te, test_merged$y_num)
  res_tr$Model <- m_name; res_tr$Cohort <- "Training"
  res_te$Model <- m_name; res_te$Cohort <- "Testing"

  nri_idi_results[[length(nri_idi_results) + 1]] <- res_tr
  nri_idi_results[[length(nri_idi_results) + 1]] <- res_te
}
df_nri_idi <- do.call(rbind, nri_idi_results)

df_nri_idi$Model_Prof <- dplyr::recode(df_nri_idi$Model,
                                       "ΔITH-only (0/1)"            = "%ΔVHS",
                                       "Clinical-only (No ΔITH)"    = "Clin_pre",
                                       "Pre-ITH only"               = "VHS_pre",
                                       "Post-ITH only"              = "VHS_post",
                                       "Radiomics (Pre) Radscore"   = "VAT-RadScore_pre",
                                       "Radiomics (Post) Radscore"  = "VAT-RadScore_post",
                                       "Radiomics (Delta) Radscore" = "ΔVAT-RadScore"
)

level_order_nri <- rev(c("%ΔVHS", "Clin_pre", "VHS_pre", "VHS_post", "VAT-RadScore_pre", "VAT-RadScore_post", "ΔVAT-RadScore"))
df_nri_idi$Model_Prof <- factor(df_nri_idi$Model_Prof, levels = level_order_nri)
df_nri_tr <- df_nri_idi[df_nri_idi$Cohort == "Training", ]
df_nri_te <- df_nri_idi[df_nri_idi$Cohort == "Testing", ]


plot_heatmap_single <- function(df, title, tag) {
  ggplot2::ggplot(df, ggplot2::aes(x = Metric, y = Model, fill = Delta)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = Label), size = 3.2, family = "serif", lineheight = 0.9) +
    ggplot2::scale_fill_gradient2(limits = c(-max_abs_heat, max_abs_heat), midpoint = 0, low="#3C5488", high="#E64B35") +
    ggplot2::labs(title = title, x = NULL, y = NULL, fill = "ΔMetric", tag = tag) +
    ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 18, color = "black"),
      plot.tag.position = c(0.015, 0.97),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13, color = "#111111", margin=ggplot2::margin(b=14)),
      axis.text.y = ggplot2::element_text(size = 11, color = "#222222", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10, color = "#222222", angle = 45, hjust = 1, vjust=1),
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(1.2, "cm"),
      legend.title = ggplot2::element_text(size=10, face="bold", vjust=0.8),
      plot.margin = ggplot2::margin(20, 10, 5, 10)
    ) + ggplot2::coord_cartesian(clip = "off")
}

plot_nri_idi_forest <- function(df, metric_col, lower_col, upper_col, p_col, title, x_label, plot_tag, bg_color, line_color="#3C5488", point_color="#E64B35") {
  df$p_lbl <- ifelse(is.na(df[[p_col]]), "NA", ifelse(df[[p_col]] < 0.001, "<0.001", sprintf("=%.3f", df[[p_col]])))
  df$lbl <- sprintf("%.3f (p%s)", df[[metric_col]], df$p_lbl)

  ggplot2::ggplot(df, ggplot2::aes(x = .data[[metric_col]], y = Model_Prof)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "#7B8B9A", linewidth = 0.8) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data[[lower_col]], xmax = .data[[upper_col]]), height = 0.15, color = line_color, linewidth = 0.9) +
    ggplot2::geom_point(color = point_color, size = 3.5) +
    ggplot2::geom_text(ggplot2::aes(label = lbl), vjust = -1.1, size = 3.6, family = "serif", color = "#222222") +
    ggplot2::labs(title = title, x = x_label, y = NULL, tag = plot_tag) +
    ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.grid.major.y = ggplot2::element_line(color = "#FFFFFF", linewidth = 0.8),
      panel.grid.major.x = ggplot2::element_line(color = "#FFFFFF", linewidth = 0.8),
      panel.grid.minor = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 18, color = "black"),
      plot.tag.position = c(0.015, 0.97),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13, color = "#111111", margin=ggplot2::margin(b=14)),
      axis.text.y = ggplot2::element_text(size = 11, color = "#222222", face = "bold"),
      axis.text.x = ggplot2::element_text(size = 10, color = "#222222"),
      axis.title.x = ggplot2::element_text(size = 11, color = "#111111", face = "bold", margin=ggplot2::margin(t=10)),
      plot.margin = ggplot2::margin(20, 10, 15, 10)
    ) + ggplot2::coord_cartesian(clip = "off")
}


p_heat_tr <- plot_heatmap_single(df_heat_tr, "Diagnostic performance relative to Nom-%ΔVHS (Training cohort)", "A")
p_nri_tr  <- plot_nri_idi_forest(df_nri_tr, "NRI", "NRI_lower", "NRI_upper", "NRI_p",
                                 "Continuous net reclassification improvement (Training cohort)",
                                 "Continuous NRI (vs. Nom-%ΔVHS)", "B", "#F4F7F9")
p_idi_tr  <- plot_nri_idi_forest(df_nri_tr, "IDI", "IDI_lower", "IDI_upper", "IDI_p",
                                 "Integrated discrimination improvement (Training cohort)",
                                 "IDI (vs. Nom-%ΔVHS)", "C", "#FCF8F5")


p_heat_te <- plot_heatmap_single(df_heat_te, "Diagnostic performance relative to Nom-%ΔVHS (Testing cohort)", "D")
p_nri_te  <- plot_nri_idi_forest(df_nri_te, "NRI", "NRI_lower", "NRI_upper", "NRI_p",
                                 "Continuous net reclassification improvement (Testing cohort)",
                                 "Continuous NRI (vs. Nom-%ΔVHS)", "E", "#F4F7F9")
p_idi_te  <- plot_nri_idi_forest(df_nri_te, "IDI", "IDI_lower", "IDI_upper", "IDI_p",
                                 "Integrated discrimination improvement (Testing cohort)",
                                 "IDI (vs. Nom-%ΔVHS)", "F", "#FCF8F5")


fig_mega_2x3 <- gridExtra::arrangeGrob(
  p_heat_tr, p_nri_tr, p_idi_tr,
  p_heat_te, p_nri_te, p_idi_te,
  nrow = 2, ncol = 3,
  widths = c(1.3, 1, 1),
  padding = grid::unit(1.2, "line")
)


bg_rect <- grid::rectGrob(gp = grid::gpar(fill = "white", col = NA))
fig_mega_final <- gridExtra::arrangeGrob(bg_rect, fig_mega_2x3, heights = grid::unit(c(0, 1), c("null", "null")))

pdf(file.path(FIG_DIR, "Mega_Fig_Heatmap_NRI_IDI_2x3.pdf"), width = 22, height = 12)
grid::grid.draw(fig_mega_final)
dev.off()

png(file.path(FIG_DIR, "Mega_Fig_Heatmap_NRI_IDI_2x3.png"), width = 5200, height = 2800, res = 240)
grid::grid.draw(fig_mega_final)
dev.off()

cat("[OK] Mega 2x3 Figure (Heatmap + NRI + IDI) exported:\n",
    " - ", basename(file.path(FIG_DIR, "Mega_Fig_Heatmap_NRI_IDI_2x3.pdf")), "\n",
    " - ", basename(file.path(FIG_DIR, "Mega_Fig_Heatmap_NRI_IDI_2x3.png")), "\n", sep = "")

cat("\n=== ALL PIPELINE STEPS COMPLETED SUCCESSFULLY ===\n")


cat("\n=== Step 14-16: Build Mega-Figure (Calibration, ROC, DCA in 2x3 Grid) ===\n")


VHS_RELATED <- c("pct_Delta_VHS","pct_Delta_VHS_cat","Delta_VHS","VHS_pre","VHS_post",
                 "pct_Delta_ITH","pct_Delta_ITH_cat","Baseline_ITH","Post_ITH")
vars_clinical_only <- setdiff(sig_vars, intersect(sig_vars, VHS_RELATED))


fit_clin <- stats::glm(as.formula(paste("group ~", paste(vars_clinical_only, collapse = " + "))),
                       data = train_data, family = binomial)

clean_train$score_clin <- as.numeric(predict(fit_clin, newdata = clean_train, type = "response"))
clean_test$score_clin  <- as.numeric(predict(fit_clin, newdata = clean_test, type = "response"))

clean_train$score_dith <- (1 - suppressWarnings(as.numeric(as.character(clean_train$pct_Delta_ITH_cat))))
clean_test$score_dith  <- (1 - suppressWarnings(as.numeric(as.character(clean_test$pct_Delta_ITH_cat))))


roc_nom_tr <- pROC::roc(clean_train$group, clean_train$nom_prob, levels = c(0,1), direction = "<", quiet = TRUE)
roc_nom_te <- pROC::roc(clean_test$group, clean_test$nom_prob, levels = c(0,1), direction = "<", quiet = TRUE)
roc_dit_tr <- pROC::roc(clean_train$group, clean_train$score_dith, levels = c(0,1), direction = "<", quiet = TRUE)
roc_dit_te <- pROC::roc(clean_test$group,  clean_test$score_dith, levels = c(0,1), direction = "<", quiet = TRUE)
roc_cln_tr <- pROC::roc(clean_train$group, clean_train$score_clin, levels = c(0,1), direction = "<", quiet = TRUE)
roc_cln_te <- pROC::roc(clean_test$group,  clean_test$score_clin, levels = c(0,1), direction = "<", quiet = TRUE)


ci_nom_tr <- as.numeric(pROC::ci.auc(roc_nom_tr))
ci_dit_tr <- as.numeric(pROC::ci.auc(roc_dit_tr))
ci_cln_tr <- as.numeric(pROC::ci.auc(roc_cln_tr))

ci_nom_te <- as.numeric(pROC::ci.auc(roc_nom_te))
ci_dit_te <- as.numeric(pROC::ci.auc(roc_dit_te))
ci_cln_te <- as.numeric(pROC::ci.auc(roc_cln_te))

extract_roc_df <- function(r, model_name) { data.frame(FPR = 1 - r$specificities, TPR = r$sensitivities, Model = model_name) }


labels_tr <- c(
  sprintf("Nom-%%ΔVHS (AUC=%.3f, 95%% CI: %.3f-%.3f)", ci_nom_tr[2], ci_nom_tr[1], ci_nom_tr[3]),
  sprintf("%%ΔVHS (AUC=%.3f, 95%% CI: %.3f-%.3f)",     ci_dit_tr[2], ci_dit_tr[1], ci_dit_tr[3]),
  sprintf("Clin_pre (AUC=%.3f, 95%% CI: %.3f-%.3f)",   ci_cln_tr[2], ci_cln_tr[1], ci_cln_tr[3])
)

labels_te <- c(
  sprintf("Nom-%%ΔVHS (AUC=%.3f, 95%% CI: %.3f-%.3f)", ci_nom_te[2], ci_nom_te[1], ci_nom_te[3]),
  sprintf("%%ΔVHS (AUC=%.3f, 95%% CI: %.3f-%.3f)",     ci_dit_te[2], ci_dit_te[1], ci_dit_te[3]),
  sprintf("Clin_pre (AUC=%.3f, 95%% CI: %.3f-%.3f)",   ci_cln_te[2], ci_cln_te[1], ci_cln_te[3])
)

df_roc_tr <- rbind(extract_roc_df(roc_nom_tr, labels_tr[1]), extract_roc_df(roc_dit_tr, labels_tr[2]), extract_roc_df(roc_cln_tr, labels_tr[3]))
df_roc_te <- rbind(extract_roc_df(roc_nom_te, labels_te[1]), extract_roc_df(roc_dit_te, labels_te[2]), extract_roc_df(roc_cln_te, labels_te[3]))
df_roc_tr$Model <- factor(df_roc_tr$Model, levels = labels_tr)
df_roc_te$Model <- factor(df_roc_te$Model, levels = labels_te)


clean_train$p_bnd <- pmin(pmax(clean_train$nom_prob, 1e-6), 1 - 1e-6)
clean_test$p_bnd  <- pmin(pmax(clean_test$nom_prob,  1e-6), 1 - 1e-6)

cal_train_nom <- tryCatch(rms::calibrate(final_lrm, method = "boot", B = 1000), error = function(e) rms::calibrate(final_lrm, method = "boot", B = 200))
cal_train_curve_nom <- .extract_rms_calibrate_curve(cal_train_nom, verbose = FALSE)

cal_train_bin_nom  <- .make_binned_calibration(clean_train, pred_col = "nom_prob", n_bins = 10)
cal_test_bin_nom   <- .make_binned_calibration(clean_test,  pred_col = "nom_prob", n_bins = 10)
cal_test_loess_nom <- .make_loess_calibration(clean_test, pred_col = "nom_prob", span = 0.9, grid_n = 200)

res_train_stats_nom <- tryCatch(rms::val.prob(p = clean_train$p_bnd, y = clean_train$group, pl = FALSE, statloc = FALSE), error = function(e) c(Brier=NA))
res_test_stats_nom  <- tryCatch(rms::val.prob(p = clean_test$p_bnd, y = clean_test$group, pl = FALSE, statloc = FALSE), error = function(e) c(Brier=NA))


dca_train <- dcurves::dca(group ~ nom_prob + score_clin + score_dith, data = clean_train, thresholds = seq(0, 0.95, by = 0.01), label = list(nom_prob="Nom-%ΔVHS", score_clin="Clin_pre", score_dith="%ΔVHS"))
dca_test  <- dcurves::dca(group ~ nom_prob + score_clin + score_dith, data = clean_test,  thresholds = seq(0, 0.95, by = 0.01), label = list(nom_prob="Nom-%ΔVHS", score_clin="Clin_pre", score_dith="%ΔVHS"))


plot_calibration_npj <- function(curve_nom, bin_nom, rug_data, stats_nom, title, tag, bg_color = "#FCF8F5") {
  ggplot2::ggplot() +

    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.9, color = "#8895A5") +

    ggplot2::geom_line(data = curve_nom, ggplot2::aes(x = pred, y = obs), color = "#E64B35", linewidth = 1.3) +

    ggplot2::geom_errorbar(data = bin_nom, ggplot2::aes(x = pred, ymin = lower, ymax = upper), width = 0.02, color = "#E64B35", linewidth = 0.8) +
    ggplot2::geom_point(data = bin_nom, ggplot2::aes(x = pred, y = obs), color = "#E64B35", size = 2.8) +
    ggplot2::geom_rug(data = rug_data, ggplot2::aes(x = nom_prob), sides = "b", alpha = 0.25, color = "#222222") +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(x = "Predicted probability", y = "Observed probability", title = title, tag = tag) +
    ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.grid.major = ggplot2::element_line(color = "#FFFFFF", linewidth = 0.9),
      panel.grid.minor = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 18, color = "black"),
      plot.tag.position = c(0.015, 0.98),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13, color = "#111111", margin=ggplot2::margin(b=14)),
      axis.text = ggplot2::element_text(size = 11, color = "#222222", face = "bold"),
      axis.title = ggplot2::element_text(size = 12, color = "#111111", face = "bold"),
      plot.margin = ggplot2::margin(20, 15, 20, 15)
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip="off") +

    ggplot2::annotate("text", x = 0.96, y = 0.04, hjust = 1, vjust = 0, size = 4.0, family = "serif", color = "#222222", fontface="bold",
                      label = sprintf("Brier = %.3f", stats_nom["Brier"]))
}


plot_npj_roc <- function(df, title, tag, bg_color = "#F4F7F9") {
  ggplot2::ggplot(df, ggplot2::aes(x = FPR, y = TPR, color = Model, linetype = Model)) +
    ggplot2::geom_line(linewidth = 1.15) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "#8895A5", linewidth = 0.8) +
    ggplot2::labs(title = title, x = "1 - Specificity", y = "Sensitivity", tag = tag) +
    ggplot2::scale_color_manual(values = c("#E64B35", "#3C5488", "#00A087")) +
    ggplot2::scale_linetype_manual(values = c("solid", "dotdash", "longdash")) +
    ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.grid.major = ggplot2::element_line(color = "#FFFFFF", linewidth = 0.9),
      panel.grid.minor = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 18, color = "black"),
      plot.tag.position = c(0.015, 0.98),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13, color = "#111111", margin=ggplot2::margin(b=14)),
      axis.text = ggplot2::element_text(size = 11, color = "#222222", face = "bold"),
      axis.title = ggplot2::element_text(size = 12, color = "#111111", face = "bold"),
      legend.justification = c(1, 0),
      legend.position = c(0.98, 0.02),
      legend.background = ggplot2::element_rect(fill = ggplot2::alpha(bg_color, 0.8), color = NA),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 9.5, color = "#222222", family="serif", face="bold"),
      plot.margin = ggplot2::margin(20, 15, 20, 15)
    ) + ggplot2::coord_cartesian(xlim=c(0,1), ylim=c(0,1), clip="off")
}


plot_advanced_dca <- function(dca_obj, title, tag, bg_color = "#F2F6F3") {
  df <- as.data.frame(dca_obj$dca)

  max_nb <- max(df$net_benefit, na.rm = TRUE)
  headroom <- max(0.12, max_nb * 0.4)
  y_lower <- -0.05
  y_upper <- max_nb + headroom

  df$label <- factor(df$label, levels = c("Nom-%ΔVHS", "Clin_pre", "%ΔVHS", "Treat All", "Treat None"))

  ggplot2::ggplot(df, ggplot2::aes(x = threshold, y = net_benefit, color = label, linetype = label)) +
    ggplot2::geom_hline(yintercept = 0, color = "#8895A5", linewidth = 0.6) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::scale_color_manual(values = c("Nom-%ΔVHS"="#E64B35", "Clin_pre"="#00A087", "%ΔVHS"="#3C5488", "Treat All"="#A0A8B4", "Treat None"="#222222")) +
    ggplot2::scale_linetype_manual(values = c("Nom-%ΔVHS"="solid", "Clin_pre"="longdash", "%ΔVHS"="dotdash", "Treat All"="dashed", "Treat None"="solid")) +
    ggplot2::coord_cartesian(xlim = c(0, 0.95), ylim = c(y_lower, y_upper), clip = "off") +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(title = title, x = "Threshold probability", y = "Net benefit", tag = tag) +
    ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.grid.major = ggplot2::element_line(color = "#FFFFFF", linewidth = 0.9),
      panel.grid.minor = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 18, color = "black"),
      plot.tag.position = c(0.015, 0.98),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13, color = "#111111", margin=ggplot2::margin(b=14)),
      axis.text = ggplot2::element_text(size = 11, color = "#222222", face = "bold"),
      axis.title = ggplot2::element_text(size = 12, color = "#111111", face = "bold"),
      legend.justification = c(1, 1),
      legend.position = c(0.98, 0.98),
      legend.background = ggplot2::element_rect(fill = ggplot2::alpha(bg_color, 0.85), color = NA),
      legend.key = ggplot2::element_blank(),
      legend.key.width = ggplot2::unit(1.2, "cm"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10, family = "serif", face = "bold"),
      plot.margin = ggplot2::margin(20, 15, 20, 15)
    )
}


cat("\n=== Step 18: Build Boxplot-Jitter Plots for VHS features (2x3 Grid) ===\n")


plot_boxplot_npj <- function(df, y_col, title, y_label, tag, bg_color, min_y = NULL) {


  d <- df[!is.na(df[[y_col]]) & !is.na(df$group), ]
  d$group_fac <- factor(d$group, levels = c(0, 1), labels = c("Response (0)", "SLR (1)"))


  if (!is.null(min_y)) {
    d <- d[d[[y_col]] >= min_y, ]
  }


  pval <- suppressWarnings(wilcox.test(d[[y_col]] ~ d$group_fac)$p.value)
  p_text <- ifelse(pval < 0.001, "P < 0.001", sprintf("P = %.3f", pval))


  y_max <- max(d[[y_col]], na.rm = TRUE)
  y_min <- min(d[[y_col]], na.rm = TRUE)
  y_range <- y_max - y_min
  y_line <- y_max + 0.06 * y_range
  y_star <- y_max + 0.12 * y_range


  p <- ggplot2::ggplot(d, ggplot2::aes(x = group_fac, y = .data[[y_col]])) +


    ggplot2::geom_boxplot(ggplot2::aes(fill = group_fac), width = 0.45, color = "#222222", outlier.shape = NA, alpha = 0.65, linewidth = 0.8) +


    ggplot2::geom_jitter(ggplot2::aes(color = group_fac), width = 0.15, size = 1.8, alpha = 0.8) +


    ggplot2::scale_fill_manual(values = c("Response (0)" = "#3C5488", "SLR (1)" = "#E64B35")) +
    ggplot2::scale_color_manual(values = c("Response (0)" = "#3C5488", "SLR (1)" = "#E64B35")) +


    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.25))) +


    ggplot2::annotate("segment", x = 1, xend = 2, y = y_line, yend = y_line, color = "#222222", linewidth = 0.7) +
    ggplot2::annotate("segment", x = 1, xend = 1, y = y_line, yend = y_line - 0.015 * y_range, color = "#222222", linewidth = 0.7) +
    ggplot2::annotate("segment", x = 2, xend = 2, y = y_line, yend = y_line - 0.015 * y_range, color = "#222222", linewidth = 0.7) +
    ggplot2::annotate("text", x = 1.5, y = y_star, label = p_text, size = 4.2, family = "serif", fontface = "bold", color = "#222222") +


    ggplot2::labs(title = title, x = NULL, y = y_label, tag = tag) +
    ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.grid.major.y = ggplot2::element_line(color = "#FFFFFF", linewidth = 1.2),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 18, color = "black"),
      plot.tag.position = c(0.015, 0.98),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13, color = "#111111", margin = ggplot2::margin(b = 15)),
      axis.text.x = ggplot2::element_text(size = 12, color = "#222222", face = "bold"),
      axis.text.y = ggplot2::element_text(size = 11, color = "#222222", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12, color = "#111111", face = "bold", margin = ggplot2::margin(r = 10)),
      legend.position = "none",
      plot.margin = ggplot2::margin(t = 20, r = 15, b = 15, l = 15)
    )

  return(p)
}


bg_pre   <- "#F4F7F9"
bg_post  <- "#FCF8F5"
bg_delta <- "#F2F6F3"


p_vhs_pre_tr  <- plot_boxplot_npj(train_data, "Baseline_ITH",  "VHS_pre\n(Training cohort)", "Score", "A", bg_pre)
p_vhs_post_tr <- plot_boxplot_npj(train_data, "Post_ITH",      "VHS_post\n(Training cohort)", "Score", "B", bg_post)
p_delta_vhs_tr<- plot_boxplot_npj(train_data, "pct_Delta_ITH", "%ΔVHS\n(Training cohort)", "Percentage change (%)", "C", bg_delta, min_y = -2000)


p_vhs_pre_te  <- plot_boxplot_npj(test_data, "Baseline_ITH",  "VHS_pre\n(Testing cohort)", "Score", "D", bg_pre)
p_vhs_post_te <- plot_boxplot_npj(test_data, "Post_ITH",      "VHS_post\n(Testing cohort)", "Score", "E", bg_post)
p_delta_vhs_te<- plot_boxplot_npj(test_data, "pct_Delta_ITH", "%ΔVHS\n(Testing cohort)", "Percentage change (%)", "F", bg_delta, min_y = -2000)


fig_box_2x3 <- gridExtra::arrangeGrob(
  p_vhs_pre_tr, p_vhs_post_tr, p_delta_vhs_tr,
  p_vhs_pre_te, p_vhs_post_te, p_delta_vhs_te,
  nrow = 2, ncol = 3,
  padding = grid::unit(1.2, "line")
)


bg_rect <- grid::rectGrob(gp = grid::gpar(fill = "white", col = NA))
fig_box_final <- gridExtra::arrangeGrob(bg_rect, fig_box_2x3, heights = grid::unit(c(0, 1), c("null", "null")))

pdf(file.path(FIG_DIR, "Mega_Fig_Boxplot_VHS_2x3.pdf"), width = 16, height = 12.5)
grid::grid.draw(fig_box_final)
dev.off()

png(file.path(FIG_DIR, "Mega_Fig_Boxplot_VHS_2x3.png"), width = 3800, height = 2900, res = 240)
grid::grid.draw(fig_box_final)
dev.off()

cat("[OK] Robust Boxplot-Jitter 2x3 Figure exported:\n",
    " - ", basename(file.path(FIG_DIR, "Mega_Fig_Boxplot_VHS_2x3.pdf")), "\n",
    " - ", basename(file.path(FIG_DIR, "Mega_Fig_Boxplot_VHS_2x3.png")), "\n", sep = "")


.theme_npj_panel <- function(bg){
  ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = bg, color = NA),
      panel.background = ggplot2::element_rect(fill = bg, color = NA),
      panel.grid.major = ggplot2::element_line(color = "#FFFFFF", linewidth = 1.15),
      panel.grid.minor = ggplot2::element_blank(),

      axis.title       = ggplot2::element_text(size = 14, face = "bold", color = "#111111"),
      axis.text        = ggplot2::element_text(size = 12.5, face = "bold", color = "#222222"),

      plot.title       = ggplot2::element_text(
        size = 14, face = "bold", hjust = 0.5, color = "#111111",
        margin = ggplot2::margin(b = 10)
      ),

      plot.tag         = ggplot2::element_text(size = 22, face = "bold", color = "black"),
      plot.tag.position = c(0.02, 0.98),

      legend.background = ggplot2::element_rect(fill = ggplot2::alpha(bg, 0.80), color = NA),
      legend.title     = ggplot2::element_blank(),
      legend.text      = ggplot2::element_text(size = 12, face = "bold"),

      plot.margin      = ggplot2::margin(18, 14, 14, 14)
    )
}


plot_paired_panel <- function(df, id_col, pre_col, post_col, title, tag, bg){
  stopifnot(all(c("group", id_col, pre_col, post_col) %in% colnames(df)))

  d <- df[, c("group", id_col, pre_col, post_col), drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if(nrow(d) < 5){
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title=paste0(title, " (Insufficient data)")))
  }

  d$group_fac <- factor(d$group, levels = c(0,1), labels = c("Response (0)", "SLR (1)"))
  d[[id_col]] <- trimws(as.character(d[[id_col]]))

  long <- rbind(
    data.frame(ID=d[[id_col]], group_fac=d$group_fac, Time="Pre",  Value=d[[pre_col]]),
    data.frame(ID=d[[id_col]], group_fac=d$group_fac, Time="Post", Value=d[[post_col]])
  )
  long$Time <- factor(long$Time, levels = c("Pre","Post"))

  mean_df <- long %>%
    dplyr::group_by(group_fac, Time) %>%
    dplyr::summarise(
      mean = mean(Value, na.rm = TRUE),
      sdv  = stats::sd(Value, na.rm = TRUE),
      n_id = dplyr::n_distinct(ID),
      .groups="drop"
    )
  mean_df$se <- mean_df$sdv / sqrt(pmax(mean_df$n_id, 1))

  d$Delta <- d[[post_col]] - d[[pre_col]]
  p_between <- tryCatch(stats::wilcox.test(Delta ~ factor(group), data=d)$p.value, error=function(e) NA_real_)
  ptxt <- .p_text(p_between)

  y_max <- max(long$Value, na.rm = TRUE)
  y_min <- min(long$Value, na.rm = TRUE)
  y_rng <- y_max - y_min; if(!is.finite(y_rng) || y_rng == 0) y_rng <- 1
  y_annot <- y_max + 0.20*y_rng

  ggplot2::ggplot(long, ggplot2::aes(x=Time, y=Value, group=ID, color=group_fac)) +
    ggplot2::geom_line(alpha=0.22, linewidth=0.85) +
    ggplot2::geom_point(alpha=0.58, size=2.0) +
    ggplot2::geom_line(
      data=mean_df,
      ggplot2::aes(x=Time, y=mean, group=group_fac, color=group_fac),
      inherit.aes=FALSE,
      linewidth=3.0
    ) +
    ggplot2::geom_errorbar(
      data=mean_df,
      ggplot2::aes(x=Time, ymin=mean-se, ymax=mean+se, color=group_fac),
      inherit.aes=FALSE,
      width=0.12,
      linewidth=1.1
    ) +
    ggplot2::scale_color_manual(values = COL_MAP) +
    ggplot2::labs(title=title, x=NULL, y="VHS score", tag=tag) +
    ggplot2::annotate(
      "text", x=1.5, y=y_annot,
      label=paste0("Δ(Post−Pre) between groups: ", ptxt),
      size=5.0, family="serif", fontface="bold", color="#111111"
    ) +
    .theme_npj_panel(bg) +
    ggplot2::theme(legend.position = c(0.86, 0.22)) +
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0.10, 0.33)))
}


plot_box_panel <- function(df, y_col, title, ylab, tag, bg, min_y=NULL){
  stopifnot(all(c("group", y_col) %in% colnames(df)))

  d <- df[!is.na(df[[y_col]]) & !is.na(df$group), , drop=FALSE]
  if(!is.null(min_y)) d <- d[d[[y_col]] >= min_y, , drop=FALSE]
  d$group_fac <- factor(d$group, levels=c(0,1), labels=c("Response (0)", "SLR (1)"))

  pval <- tryCatch(stats::wilcox.test(d[[y_col]] ~ d$group_fac)$p.value, error=function(e) NA_real_)
  ptxt <- .p_text(pval)

  y_max <- max(d[[y_col]], na.rm=TRUE)
  y_min <- min(d[[y_col]], na.rm=TRUE)
  y_rng <- y_max - y_min; if(!is.finite(y_rng) || y_rng==0) y_rng <- 1
  y_line <- y_max + 0.06*y_rng
  y_star <- y_max + 0.13*y_rng

  ggplot2::ggplot(d, ggplot2::aes(x=group_fac, y=.data[[y_col]])) +
    ggplot2::geom_boxplot(ggplot2::aes(fill=group_fac),
                          width=0.52, alpha=0.72,
                          color="#111111", linewidth=1.0,
                          outlier.shape=NA) +
    ggplot2::geom_jitter(ggplot2::aes(color=group_fac),
                         width=0.14, size=2.4, alpha=0.88) +
    ggplot2::scale_fill_manual(values=COL_MAP) +
    ggplot2::scale_color_manual(values=COL_MAP) +
    ggplot2::annotate("segment", x=1, xend=2, y=y_line, yend=y_line, linewidth=0.9, color="#111111") +
    ggplot2::annotate("text", x=1.5, y=y_star, label=ptxt,
                      size=5.0, family="serif", fontface="bold", color="#111111") +
    ggplot2::labs(title=title, x=NULL, y=ylab, tag=tag) +
    .theme_npj_panel(bg) +
    ggplot2::theme(legend.position="none") +
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0.10, 0.30)))
}


pkgs2 <- c("dplyr","gtsummary","tibble","rlang","forcats","flextable","officer")
for(p in pkgs2){
  if(!require(p, character.only = TRUE)) stop(paste("Missing required R package:", p))
  library(p, character.only = TRUE)
}


train_data$group <- as.numeric(train_data$group)


if(!("pct_Delta_ITH_cat" %in% colnames(train_data))) stop("pct_Delta_ITH_cat not found in train_data.")
train_data$pct_Delta_ITH_cat <- as.factor(train_data$pct_Delta_ITH_cat)
train_data$pct_Delta_ITH_cat <- forcats::fct_drop(train_data$pct_Delta_ITH_cat)
train_data$pct_Delta_ITH_cat <- stats::relevel(train_data$pct_Delta_ITH_cat, ref = "1")


var_labels <- c(
  "pct_Delta_ITH_cat" = "%ΔVHS category (0 vs 1; reference = 1)",
  "Age"               = "Age (years)",
  "Sex"               = "Sex",
  "BMI"               = "Body mass index (kg/m^2)",
  "Disease_location"  = "Disease location (Montreal classification)",
  "Disease_behavior"  = "Disease behavior (Montreal classification)",
  "CDAI"              = "Crohn’s Disease Activity Index (CDAI)",
  "SES_CD"            = "Simple Endoscopic Score for Crohn’s Disease (SES-CD)",
  "CRP"               = "C-reactive protein (mg/L)",
  "ESR"               = "Erythrocyte sedimentation rate (mm/h)",
  "Platelet"          = "Platelet count",
  "NLR"               = "Neutrophil-to-lymphocyte ratio (NLR)",
  "Albumin"           = "Serum albumin (g/L)",
  "Creatinine"        = "Serum creatinine",
  "SII"               = "Systemic immune-inflammation index (SII)",
  "Smoking"           = "Smoking status"
)


candidate_vars <- c("pct_Delta_ITH_cat","Age","Sex","BMI","Disease_location","Disease_behavior",
                    "CDAI","SES_CD","CRP","ESR","Platelet","NLR","Albumin","Creatinine","SII")
if("Smoking" %in% colnames(train_data)) candidate_vars <- c(candidate_vars, "Smoking")
candidate_vars <- candidate_vars[candidate_vars %in% colnames(train_data)]


cat_vars <- intersect(c("Sex","Disease_location","Disease_behavior","Smoking","pct_Delta_ITH_cat"), candidate_vars)
for(v in cat_vars){
  train_data[[v]] <- factor(trimws(as.character(train_data[[v]])))
}

train_data$pct_Delta_ITH_cat <- stats::relevel(train_data$pct_Delta_ITH_cat, ref = "1")


df_uv <- train_data %>%
  dplyr::select(group, dplyr::all_of(candidate_vars)) %>%
  as.data.frame()


tbl_uv <- gtsummary::tbl_uvregression(
  data = df_uv,
  method = glm,
  y = group,
  method.args = list(family = binomial),
  include = all_of(candidate_vars),
  exponentiate = TRUE
) %>%
  gtsummary::modify_header(label = "**Variable**") %>%
  gtsummary::modify_spanning_header(c("estimate","ci","p.value") ~ "**Univariable logistic regression**") %>%
  gtsummary::bold_p(t = 0.05)


tbl_uv <- tbl_uv %>% gtsummary::modify_table_body(
  ~ .x %>% dplyr::mutate(label = dplyr::recode(label, !!!var_labels))
)


get_lrt_p <- function(var_name) {
  dfit <- df_uv[, c("group", var_name), drop = FALSE]
  dfit <- dfit[stats::complete.cases(dfit), , drop = FALSE]
  if(nrow(dfit) < 5) return(NA_real_)
  fit1 <- glm(as.formula(paste("group ~", var_name)), data = dfit, family = binomial)
  fit0 <- glm(group ~ 1, data = dfit, family = binomial)
  anova(fit0, fit1, test = "Chisq")$`Pr(>Chi)`[2]
}

p_lrt <- sapply(candidate_vars, get_lrt_p)
p_lrt <- p_lrt[is.finite(p_lrt)]
sel_vars <- names(p_lrt[p_lrt < 0.20])


if(!("pct_Delta_ITH_cat" %in% sel_vars)) sel_vars <- c("pct_Delta_ITH_cat", sel_vars)
sel_vars <- unique(sel_vars)

df_mv <- train_data %>%
  dplyr::select(group, dplyr::all_of(sel_vars)) %>%
  stats::na.omit() %>%
  as.data.frame()


full_model <- glm(as.formula(paste("group ~", paste(sel_vars, collapse=" + "))),
                  data = df_mv, family = binomial)
lower_model <- glm(group ~ pct_Delta_ITH_cat, data = df_mv, family = binomial)

step_model <- step(full_model,
                   scope = list(lower = formula(lower_model), upper = formula(full_model)),
                   direction = "backward",
                   trace = 0)

tbl_mv <- gtsummary::tbl_regression(step_model, exponentiate = TRUE) %>%
  gtsummary::modify_header(label = "**Variable**") %>%
  gtsummary::modify_spanning_header(c("estimate","ci","p.value") ~ "**Multivariable logistic regression (AIC backward)**") %>%
  gtsummary::bold_p(t = 0.05)


tbl_mv <- tbl_mv %>% gtsummary::modify_table_body(
  ~ .x %>% dplyr::mutate(label = dplyr::recode(label, !!!var_labels))
)


tbl_merged <- gtsummary::tbl_merge(
  tbls = list(tbl_uv, tbl_mv),
  tab_spanner = c("**Univariable**", "**Multivariable**")
) %>%
  gtsummary::modify_caption("**Table. Univariable and Multivariable Logistic Regression for Secondary Loss of Response (SLR)**") %>%
  gtsummary::modify_footnote(
    all_stat_cols() ~
      "Odds ratios (ORs) are reported with 95% confidence intervals. Outcome: SLR = 1. For %ΔVHS category, the reference level is 1."
  )


out_docx <- file.path(TABLE_DIR, "Logistic_Univ_Multiv_Results.docx")

ft <- gtsummary::as_flex_table(tbl_merged)
ft <- flextable::autofit(ft)
ft <- flextable::fontsize(ft, part = "all", size = 10.5)
ft <- flextable::set_table_properties(ft, width = 1, layout = "autofit")

doc <- officer::read_docx()
doc <- officer::body_add_par(doc,
                             "Univariable and Multivariable Logistic Regression Results",
                             style = "heading 1")
doc <- officer::body_add_par(doc,
                             "Outcome: secondary loss of response (SLR, coded as 1). Reference category for %ΔVHS is level 1.",
                             style = "Normal")


doc <- flextable::body_add_flextable(doc, value = ft)

print(doc, target = out_docx)
cat("[OK] Logistic regression Word exported:\n - ", basename(out_docx), "\n", sep = "")


cat("\n=== Step 18 (2x4, NPJ): Paired + Pre/Post/Delta distributions ===\n")


.detect_id_col <- function(df){
  cand <- c("ID","id","Case","case","case_id","CaseID","PatientID","patient_id","SubjectID","subject_id")
  hit <- cand[cand %in% colnames(df)]
  if(length(hit) > 0) return(hit[1])
  return(NULL)
}


.p_text <- function(p){
  if(!is.finite(p)) return("P = NA")
  if(p < 0.001) return("P < 0.001")
  sprintf("P = %.3f", p)
}


BG_PAIR  <- "#F4F7F9"
BG_PRE   <- "#FCF8F5"
BG_POST  <- "#F2F6F3"
BG_DELTA <- "#F3F1F7"

COL_MAP <- c("Response (0)"="#3C5488", "SLR (1)"="#E64B35")


.theme_npj_panel <- function(bg){
  ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      plot.background  = ggplot2::element_rect(fill = bg, color = NA),
      panel.background = ggplot2::element_rect(fill = bg, color = NA),
      panel.grid.major = ggplot2::element_line(color = "#FFFFFF", linewidth = 1.15),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title       = ggplot2::element_text(size = 15, face = "bold", color = "#111111"),
      axis.text        = ggplot2::element_text(size = 13, face = "bold", color = "#222222"),
      plot.title       = ggplot2::element_text(size = 15, face = "bold", hjust = 0.5, color = "#111111",
                                               margin = ggplot2::margin(b = 12)),
      plot.tag         = ggplot2::element_text(size = 22, face = "bold", color = "black"),
      plot.tag.position = c(0.02, 0.98),
      legend.background = ggplot2::element_rect(fill = ggplot2::alpha(bg, 0.80), color = NA),
      legend.title     = ggplot2::element_blank(),
      legend.text      = ggplot2::element_text(size = 12, face = "bold"),
      plot.margin      = ggplot2::margin(18, 14, 14, 14)
    )
}


plot_paired_panel <- function(df, id_col, pre_col, post_col, title, tag, bg){
  stopifnot(all(c("group", id_col, pre_col, post_col) %in% colnames(df)))

  d <- df[, c("group", id_col, pre_col, post_col), drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if(nrow(d) < 5){
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::labs(title=paste0(title, " (Insufficient data)")))
  }

  d$group_fac <- factor(d$group, levels = c(0,1), labels = c("Response (0)", "SLR (1)"))
  d[[id_col]] <- trimws(as.character(d[[id_col]]))

  long <- rbind(
    data.frame(ID=d[[id_col]], group_fac=d$group_fac, Time="Pre",  Value=d[[pre_col]]),
    data.frame(ID=d[[id_col]], group_fac=d$group_fac, Time="Post", Value=d[[post_col]])
  )
  long$Time <- factor(long$Time, levels = c("Pre","Post"))

  mean_df <- long %>%
    dplyr::group_by(group_fac, Time) %>%
    dplyr::summarise(
      mean = mean(Value, na.rm = TRUE),
      sdv  = stats::sd(Value, na.rm = TRUE),
      n_id = dplyr::n_distinct(ID),
      .groups="drop"
    )
  mean_df$se <- mean_df$sdv / sqrt(pmax(mean_df$n_id, 1))

  d$Delta <- d[[post_col]] - d[[pre_col]]
  p_between <- tryCatch(stats::wilcox.test(Delta ~ factor(group), data=d)$p.value, error=function(e) NA_real_)
  ptxt <- .p_text(p_between)

  y_max <- max(long$Value, na.rm = TRUE)
  y_min <- min(long$Value, na.rm = TRUE)
  y_rng <- y_max - y_min; if(!is.finite(y_rng) || y_rng == 0) y_rng <- 1
  y_annot <- y_max + 0.18*y_rng

  ggplot2::ggplot(long, ggplot2::aes(x=Time, y=Value, group=ID, color=group_fac)) +
    ggplot2::geom_line(alpha=0.22, linewidth=0.75) +
    ggplot2::geom_point(alpha=0.55, size=1.7) +
    ggplot2::geom_line(
      data=mean_df,
      ggplot2::aes(x=Time, y=mean, group=group_fac, color=group_fac),
      inherit.aes=FALSE,
      linewidth=2.7
    ) +
    ggplot2::geom_errorbar(
      data=mean_df,
      ggplot2::aes(x=Time, ymin=mean-se, ymax=mean+se, color=group_fac),
      inherit.aes=FALSE,
      width=0.10,
      linewidth=1.0
    ) +
    ggplot2::scale_color_manual(values = COL_MAP) +
    ggplot2::labs(title=title, x=NULL, y="VHS score", tag=tag) +
    ggplot2::annotate(
      "text", x=1.5, y=y_annot,
      label=paste0(expression(Delta), "(Post−Pre) between groups: ", ptxt),
      size=5.5, family="serif", fontface="bold", color="#111111"
    ) +
    .theme_npj_panel(bg) +
    ggplot2::theme(legend.position = c(0.86, 0.22)) +
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0.10, 0.30)))
}


plot_box_panel <- function(df, y_col, title, ylab, tag, bg, min_y=NULL){
  stopifnot(all(c("group", y_col) %in% colnames(df)))

  d <- df[!is.na(df[[y_col]]) & !is.na(df$group), , drop=FALSE]
  if(!is.null(min_y)) d <- d[d[[y_col]] >= min_y, , drop=FALSE]
  d$group_fac <- factor(d$group, levels=c(0,1), labels=c("Response (0)", "SLR (1)"))

  pval <- tryCatch(stats::wilcox.test(d[[y_col]] ~ d$group_fac)$p.value, error=function(e) NA_real_)
  ptxt <- .p_text(pval)

  y_max <- max(d[[y_col]], na.rm=TRUE)
  y_min <- min(d[[y_col]], na.rm=TRUE)
  y_rng <- y_max - y_min; if(!is.finite(y_rng) || y_rng==0) y_rng <- 1
  y_line <- y_max + 0.06*y_rng
  y_star <- y_max + 0.12*y_rng

  ggplot2::ggplot(d, ggplot2::aes(x=group_fac, y=.data[[y_col]])) +
    ggplot2::geom_boxplot(ggplot2::aes(fill=group_fac),
                          width=0.50, alpha=0.70,
                          color="#111111", linewidth=0.9,
                          outlier.shape=NA) +
    ggplot2::geom_jitter(ggplot2::aes(color=group_fac),
                         width=0.14, size=2.0, alpha=0.85) +
    ggplot2::scale_fill_manual(values=COL_MAP) +
    ggplot2::scale_color_manual(values=COL_MAP) +
    ggplot2::annotate("segment", x=1, xend=2, y=y_line, yend=y_line, linewidth=0.8, color="#111111") +
    ggplot2::annotate("text", x=1.5, y=y_star, label=ptxt,
                      size=5.5, family="serif", fontface="bold", color="#111111") +
    ggplot2::labs(title=title, x=NULL, y=ylab, tag=tag) +
    .theme_npj_panel(bg) +
    ggplot2::theme(legend.position="none") +
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0.10, 0.28)))
}


id_tr <- .detect_id_col(train_data)
id_te <- .detect_id_col(test_data)
if(is.null(id_tr) || is.null(id_te)) stop("[ERROR] Cannot find ID column in train/test. Add ID column.")


pA <- plot_paired_panel(
  train_data, id_tr, "Baseline_ITH", "Post_ITH",
  "Paired trajectories of VHS (Pre → Post), Training cohort", "A", BG_PAIR
)

pB <- plot_box_panel(
  train_data, "Baseline_ITH",
  "Group-wise distribution of VHS_pre, Training cohort", "VHS score", "B", BG_PRE
)

pC <- plot_box_panel(
  train_data, "Post_ITH",
  "Group-wise distribution of VHS_post, Training cohort", "VHS score", "C", BG_POST
)

pD <- plot_box_panel(
  train_data, "pct_Delta_ITH",
  "Group-wise distribution of %ΔVHS, Training cohort", "Percentage change (%)", "D", BG_DELTA,
  min_y = -2000
)


pE <- plot_paired_panel(
  test_data, id_te, "Baseline_ITH", "Post_ITH",
  "Paired trajectories of VHS (Pre → Post), Testing cohort", "E", BG_PAIR
)

pF <- plot_box_panel(
  test_data, "Baseline_ITH",
  "Group-wise distribution of VHS_pre, Testing cohort", "VHS score", "F", BG_PRE
)

pG <- plot_box_panel(
  test_data, "Post_ITH",
  "Group-wise distribution of VHS_post, Testing cohort", "VHS score", "G", BG_POST
)

pH <- plot_box_panel(
  test_data, "pct_Delta_ITH",
  "Group-wise distribution of %ΔVHS, Testing cohort", "Percentage change (%)", "H", BG_DELTA,
  min_y = -2000
)


fig_2x4 <- gridExtra::arrangeGrob(
  pA, pB, pC, pD,
  pE, pF, pG, pH,
  nrow=2, ncol=4,
  widths = c(1.25, 1, 1, 1),
  padding = grid::unit(1.2, "line")
)

bg_rect <- grid::rectGrob(gp = grid::gpar(fill = "white", col = NA))
fig_final <- gridExtra::arrangeGrob(bg_rect, fig_2x4, heights = grid::unit(c(0, 1), c("null", "null")))

pdf(file.path(FIG_DIR, "Mega_Fig_Paired_PrePost_Delta_2x4_NPJ.pdf"), width=24, height=12.5)
grid::grid.draw(fig_final)
dev.off()

png(file.path(FIG_DIR, "Mega_Fig_Paired_PrePost_Delta_2x4_NPJ.png"), width=6400, height=3200, res=240)
grid::grid.draw(fig_final)
dev.off()

cat("[OK] Mega 2x4 (NPJ backgrounds) exported:\n",
    " - ", basename(file.path(FIG_DIR, "Mega_Fig_Paired_PrePost_Delta_2x4_NPJ.pdf")), "\n",
    " - ", basename(file.path(FIG_DIR, "Mega_Fig_Paired_PrePost_Delta_2x4_NPJ.png")), "\n", sep="")


cat("\n=== Step 19: Build Mega-Figure (Single Model: Calib -> ROC -> DCA in 2x3 Grid) ===\n")


clean_train <- train_data[!is.na(train_data$group), ] %>% as.data.frame()
clean_test  <- test_data[!is.na(test_data$group), ] %>% as.data.frame()

aligned <- align_factor_levels_for_rms(clean_train, clean_test, final_train_data, sig_vars)
clean_train <- aligned$train
clean_test  <- aligned$test


clean_train$nom_prob <- as.numeric(predict(final_lrm, newdata = clean_train, type = "fitted"))
clean_test$nom_prob  <- as.numeric(predict(final_lrm, newdata = clean_test, type = "fitted"))


roc_nom_tr <- pROC::roc(clean_train$group, clean_train$nom_prob, levels = c(0,1), direction = "<", quiet = TRUE)
roc_nom_te <- pROC::roc(clean_test$group, clean_test$nom_prob, levels = c(0,1), direction = "<", quiet = TRUE)

ci_nom_tr <- as.numeric(pROC::ci.auc(roc_nom_tr))
ci_nom_te <- as.numeric(pROC::ci.auc(roc_nom_te))

label_tr <- sprintf("Nom-%%ΔVHS\n(AUC = %.3f, 95%% CI: %.3f-%.3f)", ci_nom_tr[2], ci_nom_tr[1], ci_nom_tr[3])
label_te <- sprintf("Nom-%%ΔVHS\n(AUC = %.3f, 95%% CI: %.3f-%.3f)", ci_nom_te[2], ci_nom_te[1], ci_nom_te[3])

df_roc_tr <- data.frame(FPR = 1 - roc_nom_tr$specificities, TPR = roc_nom_tr$sensitivities, Label = label_tr)
df_roc_te <- data.frame(FPR = 1 - roc_nom_te$specificities, TPR = roc_nom_te$sensitivities, Label = label_te)


clean_train$p_bnd <- pmin(pmax(clean_train$nom_prob, 1e-6), 1 - 1e-6)
clean_test$p_bnd  <- pmin(pmax(clean_test$nom_prob,  1e-6), 1 - 1e-6)

cal_train_nom <- tryCatch(rms::calibrate(final_lrm, method = "boot", B = 1000), error = function(e) rms::calibrate(final_lrm, method = "boot", B = 200))
cal_train_curve_nom <- .extract_rms_calibrate_curve(cal_train_nom, verbose = FALSE)

cal_train_bin_nom  <- .make_binned_calibration(clean_train, pred_col = "nom_prob", n_bins = 10)
cal_test_bin_nom   <- .make_binned_calibration(clean_test,  pred_col = "nom_prob", n_bins = 10)
cal_test_loess_nom <- .make_loess_calibration(clean_test, pred_col = "nom_prob", span = 0.9, grid_n = 200)

res_train_stats_nom <- tryCatch(rms::val.prob(p = clean_train$p_bnd, y = clean_train$group, pl = FALSE, statloc = FALSE), error = function(e) c(Brier=NA))
res_test_stats_nom  <- tryCatch(rms::val.prob(p = clean_test$p_bnd,  y = clean_test$group,  pl = FALSE, statloc = FALSE), error = function(e) c(Brier=NA))


dca_train <- dcurves::dca(group ~ nom_prob, data = clean_train, thresholds = seq(0, 0.95, by = 0.01), label = list(nom_prob="Nom-%ΔVHS"))
dca_test  <- dcurves::dca(group ~ nom_prob, data = clean_test,  thresholds = seq(0, 0.95, by = 0.01), label = list(nom_prob="Nom-%ΔVHS"))


plot_calibration_npj <- function(curve_nom, bin_nom, rug_data, stats_nom, title, tag, bg_color = "#FCF8F5") {
  ggplot2::ggplot() +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.9, color = "#8895A5") +
    ggplot2::geom_line(data = curve_nom, ggplot2::aes(x = pred, y = obs), color = "#E64B35", linewidth = 1.5) +
    ggplot2::geom_errorbar(data = bin_nom, ggplot2::aes(x = pred, ymin = lower, ymax = upper), width = 0.02, color = "#E64B35", linewidth = 0.9) +
    ggplot2::geom_point(data = bin_nom, ggplot2::aes(x = pred, y = obs), color = "#E64B35", size = 4.0) +
    ggplot2::geom_rug(data = rug_data, ggplot2::aes(x = nom_prob), sides = "b", alpha = 0.25, color = "#222222") +
    ggplot2::scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(x = "Predicted probability", y = "Observed probability", title = title, tag = tag) +
    ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.grid.major = ggplot2::element_line(color = "#FFFFFF", linewidth = 1.1),
      panel.grid.minor = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 24, color = "black"),
      plot.tag.position = c(0.015, 0.98),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16, color = "#111111", margin=ggplot2::margin(b=14)),
      axis.text = ggplot2::element_text(size = 14, color = "#222222", face = "bold"),
      axis.title = ggplot2::element_text(size = 15, color = "#111111", face = "bold"),
      plot.margin = ggplot2::margin(20, 15, 20, 15)
    ) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip="off") +
    ggplot2::annotate("text", x = 0.96, y = 0.04, hjust = 1, vjust = 0, size = 5.5, family = "serif", color = "#222222", fontface="bold",
                      label = sprintf("Brier = %.3f", stats_nom["Brier"]))
}


plot_npj_roc_single <- function(df, title, tag, bg_color = "#F4F7F9") {
  label_text <- unique(df$Label)

  ggplot2::ggplot(df, ggplot2::aes(x = FPR, y = TPR)) +

    ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = TPR), fill = "#E64B35", alpha = 0.25) +
    ggplot2::geom_line(color = "#E64B35", linewidth = 1.6) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#8895A5", linewidth = 1.0) +
    ggplot2::labs(title = title, x = "1 - Specificity", y = "Sensitivity", tag = tag) +
    ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.grid.major = ggplot2::element_line(color = "#FFFFFF", linewidth = 1.1),
      panel.grid.minor = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 24, color = "black"),
      plot.tag.position = c(0.015, 0.98),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16, color = "#111111", margin=ggplot2::margin(b=14)),
      axis.text = ggplot2::element_text(size = 14, color = "#222222", face = "bold"),
      axis.title = ggplot2::element_text(size = 15, color = "#111111", face = "bold"),
      plot.margin = ggplot2::margin(20, 15, 20, 15)
    ) +
    ggplot2::coord_cartesian(xlim=c(0,1), ylim=c(0,1), clip="off") +
    ggplot2::annotate("text", x = 0.5, y = 0.15, label = label_text,
                      size = 5.5, family = "serif", fontface = "bold", color = "#E64B35", hjust = 0.5)
}


plot_advanced_dca_single <- function(dca_obj, title, tag, bg_color = "#F2F6F3") {
  df <- as.data.frame(dca_obj$dca)
  max_nb <- max(df$net_benefit, na.rm = TRUE)
  headroom <- max(0.12, max_nb * 0.3)
  y_lower <- -max_nb * 0.1
  y_upper <- max_nb + headroom

  df$label <- factor(df$label, levels = c("Nom-%ΔVHS", "Treat All", "Treat None"))


  df_model <- df[df$label == "Nom-%ΔVHS", c("threshold", "net_benefit")]
  df_all <- df[df$label == "Treat All", c("threshold", "net_benefit")]
  df_ribbon <- merge(df_model, df_all, by = "threshold", all.x = TRUE)
  names(df_ribbon) <- c("threshold", "nb_model", "nb_all")


  df_ribbon$ymin <- pmax(df_ribbon$nb_all, 0)

  df_ribbon <- df_ribbon[df_ribbon$nb_model > df_ribbon$ymin, ]

  ggplot2::ggplot(df, ggplot2::aes(x = threshold, y = net_benefit, color = label, linetype = label)) +


    ggplot2::geom_ribbon(data = df_ribbon,
                         ggplot2::aes(x = threshold, ymin = ymin, ymax = nb_model),
                         fill = "#E64B35", alpha = 0.15, inherit.aes = FALSE) +

    ggplot2::geom_line(linewidth = 1.6) +
    ggplot2::scale_color_manual(values = c("Nom-%ΔVHS"="#E64B35", "Treat All"="#8895A5", "Treat None"="#222222")) +
    ggplot2::scale_linetype_manual(values = c("Nom-%ΔVHS"="solid", "Treat All"="dashed", "Treat None"="solid")) +
    ggplot2::coord_cartesian(xlim = c(0, 0.95), ylim = c(y_lower, y_upper), clip = "off") +
    ggplot2::scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = seq(0, 1, 0.2)) +
    ggplot2::labs(title = title, x = "Threshold probability", y = "Net benefit", tag = tag) +
    ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.background = ggplot2::element_rect(fill = bg_color, color = NA),
      panel.grid.major = ggplot2::element_line(color = "#FFFFFF", linewidth = 1.1),
      panel.grid.minor = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 24, color = "black"),
      plot.tag.position = c(0.015, 0.98),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16, color = "#111111", margin=ggplot2::margin(b=14)),
      axis.text = ggplot2::element_text(size = 14, color = "#222222", face = "bold"),
      axis.title = ggplot2::element_text(size = 15, color = "#111111", face = "bold"),
      legend.justification = c(1, 1),
      legend.position = c(0.98, 0.98),
      legend.background = ggplot2::element_rect(fill = ggplot2::alpha(bg_color, 0.85), color = NA),
      legend.key = ggplot2::element_blank(),
      legend.key.width = ggplot2::unit(1.8, "cm"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 13, family = "serif", face = "bold"),
      plot.margin = ggplot2::margin(20, 15, 20, 15)
    )
}


p_cal_tr <- plot_calibration_npj(cal_train_curve_nom, cal_train_bin_nom, clean_train, res_train_stats_nom, "Calibration (Training cohort)", "A", "#FCF8F5")
p_roc_tr <- plot_npj_roc_single(df_roc_tr, "ROC Curve (Training cohort)", "B", "#F4F7F9")
p_dca_tr <- plot_advanced_dca_single(dca_train, "Decision Curve (Training cohort)", "C", "#F2F6F3")


p_cal_te <- plot_calibration_npj(cal_test_loess_nom, cal_test_bin_nom, clean_test, res_test_stats_nom, "Calibration (Testing cohort)", "D", "#FCF8F5")
p_roc_te <- plot_npj_roc_single(df_roc_te, "ROC Curve (Testing cohort)", "E", "#F4F7F9")
p_dca_te <- plot_advanced_dca_single(dca_test,  "Decision Curve (Testing cohort)", "F", "#F2F6F3")


fig_cal_roc_dca_2x3 <- gridExtra::arrangeGrob(
  p_cal_tr, p_roc_tr, p_dca_tr,
  p_cal_te, p_roc_te, p_dca_te,
  nrow = 2, ncol = 3,
  padding = grid::unit(1.5, "line")
)


bg_rect <- grid::rectGrob(gp = grid::gpar(fill = "white", col = NA))
fig_final_2x3 <- gridExtra::arrangeGrob(bg_rect, fig_cal_roc_dca_2x3, heights = grid::unit(c(0, 1), c("null", "null")))


pdf(file.path(FIG_DIR, "Mega_Fig_SingleModel_Calib_ROC_DCA_2x3.pdf"), width = 21, height = 13.5)
grid::grid.draw(fig_final_2x3)
dev.off()

png(file.path(FIG_DIR, "Mega_Fig_SingleModel_Calib_ROC_DCA_2x3.png"), width = 5800, height = 3600, res = 280)
grid::grid.draw(fig_final_2x3)
dev.off()

cat("[OK] Mega 2x3 Figure (Calibration -> ROC -> DCA for Single Model) exported:\n",
    " - ", basename(file.path(FIG_DIR, "Mega_Fig_SingleModel_Calib_ROC_DCA_2x3.pdf")), "\n",
    " - ", basename(file.path(FIG_DIR, "Mega_Fig_SingleModel_Calib_ROC_DCA_2x3.png")), "\n", sep = "")


cat("\n=== Step 20: Machine Learning Comparison Experiment (2x2 Grid) ===\n")


ml_pkgs <- c("randomForest", "e1071", "rpart", "xgboost", "nnet", "tidyr", "writexl")
for(p in ml_pkgs){
  if(!require(p, character.only = TRUE)) stop(paste("Missing required R package:", p))
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}


ml_train <- clean_train[, c("group", sig_vars)]
ml_test  <- clean_test[, c("group", sig_vars)]

ml_train <- na.omit(ml_train)
ml_test  <- na.omit(ml_test)


ml_train$group_fac <- factor(ml_train$group, levels = c(0, 1), labels = c("Class0", "Class1"))
ml_test$group_fac  <- factor(ml_test$group, levels = c(0, 1), labels = c("Class0", "Class1"))


predictions_tr <- list()
predictions_te <- list()


predictions_tr[["Nom-% ΔVHS"]] <- as.numeric(predict(final_lrm, newdata = ml_train, type = "fitted"))
predictions_te[["Nom-% ΔVHS"]] <- as.numeric(predict(final_lrm, newdata = ml_test,  type = "fitted"))


set.seed(123)
fit_rf <- randomForest::randomForest(group_fac ~ . - group, data = ml_train,
                                     ntree = 15, maxnodes = 4, nodesize = 20, probability = TRUE)
predictions_tr[["Random Forest"]] <- predict(fit_rf, newdata = ml_train, type = "prob")[, "Class1"]
predictions_te[["Random Forest"]] <- predict(fit_rf, newdata = ml_test,  type = "prob")[, "Class1"]


set.seed(123)
fit_svm <- e1071::svm(group_fac ~ . - group, data = ml_train,
                      kernel = "sigmoid", cost = 0.01, gamma = 0.001, probability = TRUE)
predictions_tr[["SVM"]] <- attr(predict(fit_svm, newdata = ml_train, probability = TRUE), "probabilities")[, "Class1"]
predictions_te[["SVM"]] <- attr(predict(fit_svm, newdata = ml_test,  probability = TRUE), "probabilities")[, "Class1"]


set.seed(123)
fit_dt <- rpart::rpart(group_fac ~ . - group, data = ml_train, method = "class",
                       control = rpart::rpart.control(cp = 0.08, maxdepth = 3))
predictions_tr[["Decision Tree"]] <- predict(fit_dt, newdata = ml_train, type = "prob")[, "Class1"]
predictions_te[["Decision Tree"]] <- predict(fit_dt, newdata = ml_test,  type = "prob")[, "Class1"]


train_mat <- model.matrix(group ~ . - group_fac - 1, data = ml_train)
test_mat  <- model.matrix(group ~ . - group_fac - 1, data = ml_test)

set.seed(123)
fit_xgb <- tryCatch({
  suppressWarnings(
    xgboost::xgboost(
      x = train_mat,
      y = as.logical(ml_train$group),
      max_depth = 2,
      learning_rate = 0.05,
      nrounds = 20,
      min_child_weight = 5,
      gamma = 0.5,
      objective = "binary:logistic",
      verbosity = 0
    )
  )
}, error = function(e) {
  suppressWarnings(
    xgboost::xgboost(
      data = train_mat,
      label = as.numeric(ml_train$group),
      max_depth = 2,
      eta = 0.05,
      nrounds = 20,
      min_child_weight = 5,
      gamma = 0.5,
      objective = "binary:logistic",
      verbose = 0
    )
  )
})
predictions_tr[["XGBoost"]] <- as.numeric(predict(fit_xgb, train_mat))
predictions_te[["XGBoost"]] <- as.numeric(predict(fit_xgb, test_mat))


set.seed(123)

fit_nn <- nnet::nnet(group_fac ~ . - group, data = ml_train,
                     size = 1, decay = 1.5, maxit = 80, trace = FALSE)
predictions_tr[["Neural Network"]] <- as.numeric(predict(fit_nn, newdata = ml_train, type = "raw"))
predictions_te[["Neural Network"]] <- as.numeric(predict(fit_nn, newdata = ml_test,  type = "raw"))


calc_ml_metrics <- function(y_true, y_prob, model_name, cohort, thr = NULL) {
  roc_obj <- pROC::roc(y_true, y_prob, levels = c(0,1), direction = "<", quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))


  if(is.null(thr)){
    thr <- as.numeric(pROC::coords(roc_obj, "best", ret = "threshold", best.method = "youden")[1])
  }
  y_pred <- ifelse(y_prob >= thr, 1, 0)

  TP <- sum(y_pred == 1 & y_true == 1)
  TN <- sum(y_pred == 0 & y_true == 0)
  FP <- sum(y_pred == 1 & y_true == 0)
  FN <- sum(y_pred == 0 & y_true == 1)

  acc <- (TP + TN) / (TP + TN + FP + FN)
  prec <- ifelse((TP + FP) == 0, 0, TP / (TP + FP))
  rec <- ifelse((TP + FN) == 0, 0, TP / (TP + FN))
  f1 <- ifelse((prec + rec) == 0, 0, 2 * prec * rec / (prec + rec))

  data.frame(
    Model = model_name,
    Cohort = cohort,
    Threshold = round(thr, 4),
    AUC = round(auc_val, 4),
    Accuracy = round(acc, 4),
    Precision = round(prec, 4),
    Recall = round(rec, 4),
    F1_Score = round(f1, 4),
    stringsAsFactors = FALSE
  )
}

res_list <- list()
model_names <- names(predictions_tr)

for(m in model_names){

  roc_tr_m <- pROC::roc(ml_train$group, predictions_tr[[m]], levels = c(0,1), direction = "<", quiet = TRUE)
  thr_m   <- as.numeric(pROC::coords(roc_tr_m, "best", ret = "threshold", best.method = "youden")[1])

  res_list[[length(res_list) + 1]] <- calc_ml_metrics(ml_train$group, predictions_tr[[m]], m, "Training", thr = thr_m)
  res_list[[length(res_list) + 1]] <- calc_ml_metrics(ml_test$group,  predictions_te[[m]], m, "Testing",  thr = thr_m)
}
ml_performance <- do.call(rbind, res_list)


out_ml_excel <- file.path(TABLE_DIR, "Machine_Learning_Comparison_Performance.xlsx")
writexl::write_xlsx(ml_performance, out_ml_excel)
cat("[OK] ML Comparison Performance Excel saved to:\n - ", basename(out_ml_excel), "\n")


ml_colors <- c(
  "Nom-% ΔVHS" = "#E64B35",
  "Random Forest" = "#4DBBD5",
  "SVM" = "#00A087",
  "Decision Tree" = "#F39B7F",
  "XGBoost" = "#3C5488",
  "Neural Network" = "#8491B4"
)

ml_line_types <- c(
  "Nom-% ΔVHS" = "solid",
  "Random Forest" = "dashed",
  "SVM" = "dotted",
  "Decision Tree" = "dotdash",
  "XGBoost" = "longdash",
  "Neural Network" = "twodash"
)


build_roc_data <- function(y_true, pred_list, models) {
  df_list <- list()
  auc_labs <- character(length(models))
  names(auc_labs) <- models
  for(m in models){
    r <- pROC::roc(y_true, pred_list[[m]], levels = c(0,1), direction = "<", quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(r))
    auc_labs[m] <- sprintf("%s (AUC=%.3f)", m, auc_val)
    df_list[[length(df_list) + 1]] <- data.frame(FPR = 1 - r$specificities, TPR = r$sensitivities, Model = m)
  }
  df_all <- do.call(rbind, df_list)
  list(df = df_all, labels = auc_labs)
}

roc_tr_obj <- build_roc_data(ml_train$group, predictions_tr, model_names)
roc_te_obj <- build_roc_data(ml_test$group,  predictions_te, model_names)


model_factor_levels <- c("Nom-% ΔVHS", "Random Forest", "SVM", "Decision Tree", "XGBoost", "Neural Network")
roc_tr_obj$df$Model <- factor(roc_tr_obj$df$Model, levels = model_factor_levels)
roc_te_obj$df$Model <- factor(roc_te_obj$df$Model, levels = model_factor_levels)


plot_ml_roc <- function(roc_data, auc_labels, title, tag) {
  ggplot2::ggplot(roc_data, ggplot2::aes(x = FPR, y = TPR, color = Model, linetype = Model)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "#8895A5", linewidth = 0.8) +
    ggplot2::scale_color_manual(values = ml_colors, labels = auc_labels) +
    ggplot2::scale_linetype_manual(values = ml_line_types, labels = auc_labels) +
    ggplot2::labs(title = title, x = "1 - Specificity", y = "Sensitivity", tag = tag) +
    ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "#F4F7F9", color = NA),
      panel.background = ggplot2::element_rect(fill = "#F4F7F9", color = NA),
      panel.grid.major = ggplot2::element_line(color = "#FFFFFF", linewidth = 1.0),
      panel.grid.minor = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 22, color = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15, color = "#111111", margin=ggplot2::margin(b=12)),
      axis.text = ggplot2::element_text(size = 13, color = "#222222", face = "bold"),
      axis.title = ggplot2::element_text(size = 14, color = "#111111", face = "bold"),
      legend.position = c(0.65, 0.20),
      legend.background = ggplot2::element_rect(fill = ggplot2::alpha("#F4F7F9", 0.7), color = NA),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11, face = "bold", family="serif")
    ) + ggplot2::coord_cartesian(xlim=c(0,1), ylim=c(0,1), expand=FALSE) +
    ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 20, 20))
}


build_bar_data <- function(perf_df, cohort_name) {
  sub_df <- perf_df[perf_df$Cohort == cohort_name, ]
  long_df <- tidyr::pivot_longer(sub_df, cols = c("AUC", "Accuracy", "Precision", "Recall", "F1_Score"), names_to = "Metric", values_to = "Value")
  long_df$Metric <- factor(long_df$Metric, levels = c("AUC", "Accuracy", "Precision", "Recall", "F1_Score"))
  long_df$Model <- factor(long_df$Model, levels = model_factor_levels)
  long_df
}

bar_tr_data <- build_bar_data(ml_performance, "Training")
bar_te_data <- build_bar_data(ml_performance, "Testing")


plot_ml_bar <- function(bar_data, title, tag) {
  ggplot2::ggplot(bar_data, ggplot2::aes(x = Metric, y = Value, fill = Model)) +
    ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.8), width = 0.7, color = "#222222", linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = ml_colors) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, 0.2), expand = ggplot2::expansion(mult = c(0, 0.1))) +
    ggplot2::labs(title = title, x = NULL, y = "Value", tag = tag) +
    ggplot2::theme_minimal(base_family = "serif") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "#FCF8F5", color = NA),
      panel.background = ggplot2::element_rect(fill = "#FCF8F5", color = NA),
      panel.grid.major.y = ggplot2::element_line(color = "#FFFFFF", linewidth = 1.0),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.tag = ggplot2::element_text(face = "bold", size = 22, color = "black"),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 15, color = "#111111", margin=ggplot2::margin(b=12)),
      axis.text.x = ggplot2::element_text(size = 12, color = "#222222", face = "bold"),
      axis.text.y = ggplot2::element_text(size = 12, color = "#222222", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 14, color = "#111111", face = "bold"),
      legend.position = "right",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11, face = "bold", family="serif")
    ) + ggplot2::theme(plot.margin = ggplot2::margin(20, 20, 20, 20))
}


p_roc_tr <- plot_ml_roc(roc_tr_obj$df, roc_tr_obj$labels, "ROC Comparison (Training cohort)", "A")
p_bar_tr <- plot_ml_bar(bar_tr_data, "Diagnostic Metrics (Training cohort)", "B")

p_roc_te <- plot_ml_roc(roc_te_obj$df, roc_te_obj$labels, "ROC Comparison (Testing cohort)", "C")
p_bar_te <- plot_ml_bar(bar_te_data, "Diagnostic Metrics (Testing cohort)", "D")


fig_ml_2x2 <- gridExtra::arrangeGrob(
  p_roc_tr, p_bar_tr,
  p_roc_te, p_bar_te,
  ncol = 2, nrow = 2,
  widths = c(1, 1.25),
  padding = grid::unit(1.5, "line")
)

bg_rect_ml <- grid::rectGrob(gp = grid::gpar(fill = "white", col = NA))
fig_ml_final <- gridExtra::arrangeGrob(bg_rect_ml, fig_ml_2x2, heights = grid::unit(c(0, 1), c("null", "null")))

pdf(file.path(FIG_DIR, "Mega_Fig_MachineLearning_Comparison_2x2.pdf"), width = 17, height = 13.5)
grid::grid.draw(fig_ml_final)
dev.off()

png(file.path(FIG_DIR, "Mega_Fig_MachineLearning_Comparison_2x2.png"), width = 4500, height = 3600, res = 260)
grid::grid.draw(fig_ml_final)
dev.off()

cat("[OK] ML Comparison Mega 2x2 Figure exported:\n",
    " - ", basename(file.path(FIG_DIR, "Mega_Fig_MachineLearning_Comparison_2x2.pdf")), "\n",
    " - ", basename(file.path(FIG_DIR, "Mega_Fig_MachineLearning_Comparison_2x2.png")), "\n", sep = "")

cat("\n=== Machine Learning Comparison Successfully Completed! ===\n")


