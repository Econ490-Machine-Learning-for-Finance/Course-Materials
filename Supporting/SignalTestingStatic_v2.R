# =============================================================================
# SignalTesting-Static.R (tidied, fully integrated)
# Purpose: Evaluate a static trading signal on a single cross-section
# Authors: Mike Aguilar and Ziming Huang
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pROC)
  library(moments)
})

# -----------------------------------------------------------------------------
# 0) Utilities
# -----------------------------------------------------------------------------

assert_cols <- function(df, cols, df_name = deparse(substitute(df))) {
  missing <- setdiff(cols, names(df))
  if (length(missing) > 0) {
    stop(sprintf("%s is missing columns: %s", df_name, paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}

safe_slug <- function(x) gsub("[^A-Za-z0-9_\\-]+", "_", x)

standardize_returns <- function(x, return_type = c("simple", "log")) {
  # Convert returns to SIMPLE returns for all classification / confusion / portfolio math
  return_type <- match.arg(return_type)
  x <- as.numeric(x)
  if (return_type == "log") expm1(x) else x
}

classify_outcome <- function(ret, return_threshold = 0,
                             mode = c("direction", "threshold")) {
  # direction: Positive iff ret > 0, Negative otherwise
  # threshold: Positive iff ret >= return_threshold, Negative otherwise
  mode <- match.arg(mode)
  if (mode == "direction") ifelse(ret > 0, "Positive", "Negative")
  else ifelse(ret >= return_threshold, "Positive", "Negative")
}

# -----------------------------------------------------------------------------
# 1) Distribution summaries
# -----------------------------------------------------------------------------

calculate_sample_moments <- function(ret) {
  ret <- ret[is.finite(ret)]
  n <- length(ret)
  if (n < 2) return(c(Mean = NA, Sigma = NA, Skew = NA, Kurt = NA, N = n))
  
  mu <- mean(ret)
  sig <- sd(ret)
  
  skew <- moments::skewness(ret)
  kurt_excess <- moments::kurtosis(ret) - 3
  
  c(Mean = mu, Sigma = sig, Skew = skew, Kurt = kurt_excess, N = n)
}

calculate_sample_quantiles <- function(ret) {
  ret <- ret[is.finite(ret)]
  if (length(ret) < 1) return(rep(NA_real_, 11))
  
  probs <- c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1)
  q <- quantile(ret, probs = probs, na.rm = TRUE)
  names(q) <- c("Min", ".01", ".05", ".1", ".25", ".5", ".75", ".9", ".95", ".99", "Max")
  q
}

compare_density_distributions <- function(ret1, ret2) {
  ret1 <- ret1[is.finite(ret1)]
  ret2 <- ret2[is.finite(ret2)]
  if (length(ret1) < 2 || length(ret2) < 2) return(c(diff_mean_pval = NA, ks_pval = NA))
  
  t_p <- tryCatch(t.test(ret1, ret2)$p.value, error = function(e) NA_real_)
  ks_p <- tryCatch(ks.test(ret1, ret2)$p.value, error = function(e) NA_real_)
  c(diff_mean_pval = t_p, ks_pval = ks_p)
}

# -----------------------------------------------------------------------------
# 2) Information coefficient (two variants)
# -----------------------------------------------------------------------------

calculate_information_coefficient <- function(signal_position, ret) {
  valid <- is.finite(signal_position) & is.finite(ret)
  s <- as.numeric(signal_position[valid])
  r <- as.numeric(ret[valid])
  n <- length(s)
  
  if (n < 3) {
    return(list(ic = NA, t_stat = NA, p_value = NA, n = n,
                mean_s = NA, mean_r = NA, cov_sr = NA, sigma_s = NA, sigma_r = NA))
  }
  
  mean_s <- mean(s); mean_r <- mean(r)
  cov_sr <- sum((s - mean_s) * (r - mean_r)) / n
  sigma_s <- sqrt(sum((s - mean_s)^2) / n)
  sigma_r <- sqrt(sum((r - mean_r)^2) / n)
  
  ic <- if (sigma_s > 0 && sigma_r > 0) cov_sr / (sigma_s * sigma_r) else NA_real_
  
  if (!is.na(ic) && abs(ic) < 1) {
    t_stat <- ic * sqrt(n - 2) / sqrt(1 - ic^2)
    p_val  <- 2 * pt(-abs(t_stat), df = n - 2)
  } else {
    t_stat <- NA_real_
    p_val  <- NA_real_
  }
  
  list(ic = ic, t_stat = t_stat, p_value = p_val, n = n,
       mean_s = mean_s, mean_r = mean_r, cov_sr = cov_sr,
       sigma_s = sigma_s, sigma_r = sigma_r)
}

calculate_ic_acted_only <- function(signal_position, ret) {
  valid <- is.finite(signal_position) & is.finite(ret)
  s <- as.numeric(signal_position[valid])
  r <- as.numeric(ret[valid])
  
  acted <- (s != 0)
  if (sum(acted) < 3) {
    return(list(ic = NA, t_stat = NA, p_value = NA,
                n_acted = sum(acted), n_plus = 0, n_minus = 0,
                mean_plus = NA, mean_minus = NA, sd_acted = NA))
  }
  
  s_a <- s[acted]
  r_a <- r[acted]
  
  n_plus  <- sum(s_a ==  1)
  n_minus <- sum(s_a == -1)
  n_acted <- n_plus + n_minus
  p <- n_plus / n_acted
  
  mean_plus  <- mean(r_a[s_a ==  1], na.rm = TRUE)
  mean_minus <- mean(r_a[s_a == -1], na.rm = TRUE)
  sd_acted   <- sd(r_a, na.rm = TRUE)
  
  ic <- if (sd_acted > 0 && p > 0 && p < 1) {
    (mean_plus - mean_minus) / sd_acted * sqrt(p * (1 - p))
  } else NA_real_
  
  if (!is.na(ic) && abs(ic) < 1 && n_acted > 2) {
    t_stat <- ic * sqrt(n_acted - 2) / sqrt(1 - ic^2)
    p_val  <- 2 * pt(-abs(t_stat), df = n_acted - 2)
  } else {
    t_stat <- NA_real_
    p_val  <- NA_real_
  }
  
  list(ic = ic, t_stat = t_stat, p_value = p_val,
       n_acted = n_acted, n_plus = n_plus, n_minus = n_minus,
       mean_plus = mean_plus, mean_minus = mean_minus, sd_acted = sd_acted)
}

# -----------------------------------------------------------------------------
# 3) Hit rate (null-based SE; matches slide-2)
# -----------------------------------------------------------------------------

calculate_hit_rate <- function(signal_action, ret, p0 = 0.5, one_sided = FALSE) {
  valid <- is.finite(signal_action) & is.finite(ret)
  s <- as.numeric(signal_action[valid])
  r <- as.numeric(ret[valid])
  n <- length(s)
  
  if (n == 0) {
    return(list(hit_rate_overall = NA, overall_hits = NA, overall_n = 0,
                z_overall = NA, p_overall = NA,
                hit_rate_acted = NA, acted_hits = NA, acted_n = 0,
                z_acted = NA, p_acted = NA, coverage = NA))
  }
  
  r_sign <- ifelse(r > 0, 1, -1)
  hits_overall <- sum(s * r_sign > 0, na.rm = TRUE)
  p_hat <- hits_overall / n
  
  se0 <- sqrt(p0 * (1 - p0) / n)
  z  <- (p_hat - p0) / se0
  p_val <- if (one_sided) (1 - pnorm(z)) else (2 * pnorm(-abs(z)))
  
  acted <- (s != 0)
  n_acted <- sum(acted)
  if (n_acted > 0) {
    s_a <- s[acted]
    r_a <- r[acted]
    r_a_sign <- ifelse(r_a > 0, 1, -1)
    hits_acted <- sum(s_a * r_a_sign > 0, na.rm = TRUE)
    
    p_hat_a <- hits_acted / n_acted
    se0_a <- sqrt(p0 * (1 - p0) / n_acted)
    z_a <- (p_hat_a - p0) / se0_a
    p_a <- if (one_sided) (1 - pnorm(z_a)) else (2 * pnorm(-abs(z_a)))
  } else {
    hits_acted <- NA; p_hat_a <- NA; z_a <- NA; p_a <- NA
  }
  
  list(
    hit_rate_overall = p_hat, overall_hits = hits_overall, overall_n = n,
    z_overall = z, p_overall = p_val,
    hit_rate_acted = p_hat_a, acted_hits = hits_acted, acted_n = n_acted,
    z_acted = z_a, p_acted = p_a,
    coverage = n_acted / n
  )
}

# -----------------------------------------------------------------------------
# 4) Confusion matrices (directional label by default)
# -----------------------------------------------------------------------------

create_action_confusion_matrix <- function(signal_action, ret,
                                           return_threshold = 0,
                                           outcome_mode = c("direction", "threshold")) {
  outcome_mode <- match.arg(outcome_mode)
  
  valid <- is.finite(signal_action) & is.finite(ret)
  s <- as.numeric(signal_action[valid])
  r <- as.numeric(ret[valid])
  if (length(s) < 2) {
    return(list(confusion_matrix = NULL, metrics = rep(NA, 7), counts = rep(NA, 4)))
  }
  
  prediction <- ifelse(s != 0, "Action", "No_Action")
  outcome <- classify_outcome(r, return_threshold, outcome_mode)
  
  cm <- table(
    Prediction = factor(prediction, levels = c("Action", "No_Action")),
    Outcome    = factor(outcome,    levels = c("Positive", "Negative"))
  )
  
  TP <- cm["Action",    "Positive"]
  FP <- cm["Action",    "Negative"]
  TN <- cm["No_Action", "Negative"]
  FN <- cm["No_Action", "Positive"]
  
  accuracy    <- (TP + TN) / (TP + TN + FP + FN)
  precision   <- if (TP + FP > 0) TP / (TP + FP) else NA_real_
  recall      <- if (TP + FN > 0) TP / (TP + FN) else NA_real_
  specificity <- if (TN + FP > 0) TN / (TN + FP) else NA_real_
  type1_error <- if (TN + FP > 0) FP / (TN + FP) else NA_real_
  type2_error <- if (TP + FN > 0) FN / (TP + FN) else NA_real_
  f1_score    <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0)
    2 * precision * recall / (precision + recall) else NA_real_
  
  list(
    confusion_matrix = cm,
    metrics = c(Accuracy = accuracy, Precision = precision, Recall = recall,
                Specificity = specificity, Type1_Error = type1_error,
                Type2_Error = type2_error, F1_Score = f1_score),
    counts  = c(TP = TP, FP = FP, TN = TN, FN = FN)
  )
}

create_longshort_confusion_matrix <- function(signal_action, ret,
                                              return_threshold = 0,
                                              outcome_mode = c("direction", "threshold")) {
  outcome_mode <- match.arg(outcome_mode)
  
  valid <- is.finite(signal_action) & is.finite(ret)
  s <- as.numeric(signal_action[valid])
  r <- as.numeric(ret[valid])
  
  acted <- (s != 0)
  if (sum(acted) < 2) {
    return(list(confusion_matrix = NULL, metrics = rep(NA, 7), counts = rep(NA, 4)))
  }
  
  s_a <- s[acted]
  r_a <- r[acted]
  
  prediction <- ifelse(s_a == 1, "Long", "Short")
  outcome <- classify_outcome(r_a, return_threshold, outcome_mode)
  
  cm <- table(
    Prediction = factor(prediction, levels = c("Long", "Short")),
    Outcome    = factor(outcome,    levels = c("Positive", "Negative"))
  )
  
  TP <- cm["Long",  "Positive"]
  FP <- cm["Long",  "Negative"]
  TN <- cm["Short", "Negative"]
  FN <- cm["Short", "Positive"]
  
  accuracy    <- (TP + TN) / (TP + TN + FP + FN)
  precision   <- if (TP + FP > 0) TP / (TP + FP) else NA_real_
  recall      <- if (TP + FN > 0) TP / (TP + FN) else NA_real_
  specificity <- if (TN + FP > 0) TN / (TN + FP) else NA_real_
  type1_error <- if (TN + FP > 0) FP / (TN + FP) else NA_real_
  type2_error <- if (TP + FN > 0) FN / (TP + FN) else NA_real_
  f1_score    <- if (!is.na(precision) && !is.na(recall) && (precision + recall) > 0)
    2 * precision * recall / (precision + recall) else NA_real_
  
  list(
    confusion_matrix = cm,
    metrics = c(Accuracy = accuracy, Precision = precision, Recall = recall,
                Specificity = specificity, Type1_Error = type1_error,
                Type2_Error = type2_error, F1_Score = f1_score),
    counts  = c(TP = TP, FP = FP, TN = TN, FN = FN)
  )
}

# -----------------------------------------------------------------------------
# 5) ROC / AUC (connected to confusion-matrix definitions)
# -----------------------------------------------------------------------------
# IMPORTANT:
# - roc_mode="action": prediction rule is "Action" iff |signal_matrix| > tau, so ROC uses abs(signal_matrix)
# - roc_mode="longshort": prediction uses signed signal_matrix; typically restrict to "acted" obs
# - response (label) matches slides: Positive iff ret > 0 (or threshold mode if chosen)

create_roc_analysis <- function(signal_matrix, ret,
                                roc_mode = c("action", "longshort"),
                                return_threshold = 0,  # used only if outcome_mode="threshold"
                                outcome_mode = c("direction", "threshold"),
                                plot_title = "ROC Curve",
                                acted_only = FALSE,
                                signal_threshold = 0,  # used only for longshort acted_only
                                drop_zero_returns = TRUE,
                                direction = c("auto", ">", "<")) {
  
  roc_mode     <- match.arg(roc_mode)
  outcome_mode <- match.arg(outcome_mode)
  direction    <- match.arg(direction)
  
  valid <- is.finite(signal_matrix) & is.finite(ret)
  s0 <- as.numeric(signal_matrix[valid])
  r0 <- as.numeric(ret[valid])
  
  if (drop_zero_returns) {
    keep0 <- r0 != 0
    s0 <- s0[keep0]
    r0 <- r0[keep0]
  }
  
  # Use absolute signal_matrix for Action-vs-NoAction (since Action groups + and - together)
  s <- if (roc_mode == "action") abs(s0) else s0
  r <- r0
  
  # For longshort ROC, typically evaluate only where the model actually “acts”:
  # acted iff |original signed signal_matrix| > signal_threshold
  if (roc_mode == "longshort" && acted_only) {
    keep <- abs(s0) > signal_threshold
    s <- s[keep]
    r <- r[keep]
  }
  
  if (length(s) < 5 || length(unique(s)) < 2) {
    return(list(plot = NULL, auc = NA, best_threshold = NA, youden_index = NA, roc_object = NULL))
  }
  
  response <- classify_outcome(r, return_threshold, outcome_mode)
  response <- factor(response, levels = c("Negative", "Positive"))
  
  roc_obj <- pROC::roc(
    response  = response,
    predictor = s,
    levels    = c("Negative", "Positive"),
    direction = direction,
    quiet     = TRUE
  )
  
  auc <- as.numeric(pROC::auc(roc_obj))
  
  best <- pROC::coords(
    roc_obj, x = "best", best.method = "youden",
    ret = c("threshold", "specificity", "sensitivity"),
    transpose = FALSE
  )
  
  
  
  best_thr <- as.numeric(best["threshold"])
  spec <- as.numeric(best["specificity"])
  sens <- as.numeric(best["sensitivity"])
  J <- sens + spec - 1
  
  curve <- pROC::coords(
    roc_obj, x = "all",
    ret = c("specificity", "sensitivity")
  )
  
  plot_df <- data.frame(
    FPR = 1 - curve$specificity,
    TPR = curve$sensitivity
  ) %>%
    arrange(FPR, TPR) %>%
    distinct(FPR, TPR, .keep_all = TRUE)
  
  p <- ggplot(plot_df, aes(FPR, TPR)) +
    geom_line(linewidth = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5) +
    annotate("point", x = 1 - spec, y = sens, size = 3) +
    labs(
      title = plot_title,
      subtitle = paste("AUC =", round(auc, 3),
                       "| Optimal Threshold =", round(best_thr, 3),
                       "| Youden's J =", round(J, 3)),
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme_minimal() +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  
  list(plot = p, auc = auc, best_threshold = best_thr, youden_index = J, roc_object = roc_obj)
}

# -----------------------------------------------------------------------------
# 6) MAIN
# -----------------------------------------------------------------------------

# SignalTestingStatic()
# Evaluate a static, cross-sectional trading signal against realized returns.
#
# Notation (matches slides):
#   i indexes assets (tickers) in a single cross-section.
#   s_i ∈ {-1, 0, +1} is the DISCRETE action: short / no-action / long.
#   x_i ∈ ℝ is the CONTINUOUS signal_matrix used to rank assets (e.g., momentum).
#   r_i is the realized test-period return for asset i (converted to SIMPLE returns).
#
# Two evaluation tasks:
#   (1) Action vs No Action:
#       Predict "Action" iff s_i ≠ 0. ROC uses |x_i| (strength irrespective of sign).
#   (2) Long vs Short (acted only):
#       Restrict to acted assets (|x_i| > signal_threshold or equivalently s_i ≠ 0),
#       then predict direction; ROC uses x_i (signed).
#
# Outcome definition (rows of confusion matrices; labels for ROC):
#   outcome_mode = "direction": Positive iff r_i > 0 (matches deck)
#   outcome_mode = "threshold": Positive iff r_i ≥ return_threshold (e.g. ≥ 10 bps)
#
SignalTestingStatic <- function(
    signal_position,            # data.frame with columns: ticker, signal_position   (discrete s_i ∈ {-1,0,+1})
    signal_matrix,        # data.frame with columns: ticker, signal_matrix   (continuous x_i ∈ ℝ)
    returns,           # data.frame with columns: ticker, returns       (realized r_i; simple or log)
    Meta,              # list with: assetname, benchmarkname, signalname (used for filenames/labels)
    signal_threshold = 0.001,   # τ for "acted" filter: treat as acted if |x_i| > τ (Long vs Short ROC)
    return_threshold = 0.001,   # κ for outcome_mode="threshold": Positive iff r_i ≥ κ
    outcome_mode = c("direction", "threshold"),
    return_type  = c("simple", "log"),
    one_sided_hit = FALSE,      # if TRUE, test H_A: p > 0.5; else two-sided H_A: p ≠ 0.5
    seed = NULL,                 # optional reproducibility
    output_dir = "SignalTestingStatic"
){
  
  outcome_mode <- match.arg(outcome_mode)
  return_type  <- match.arg(return_type)
  
  if (!is.null(seed)) set.seed(seed)
  
  # --- Validate inputs
  assert_cols(signal_position,     c("ticker", "signal_position"), "signal_position")
  assert_cols(signal_matrix, c("ticker", "signal_matrix"), "signal_matrix")
  assert_cols(returns,    c("ticker", "returns"),     "returns")
  
  stopifnot(
    is.list(Meta),
    all(c("assetname", "benchmarkname", "signalname") %in% names(Meta))
  )
  
  asset_name  <- Meta$assetname
  bench_name  <- Meta$benchmarkname
  signal_name <- Meta$signalname
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  cat("==========================================\n")
  cat("SIGNAL TESTING - STATIC ANALYSIS\n")
  cat("==========================================\n")
  cat("Asset:        ", asset_name, "\n", sep = "")
  cat("Benchmark:    ", bench_name, "\n", sep = "")
  cat("Signal:       ", signal_name, "\n", sep = "")
  cat("Outcome mode: ", outcome_mode, "\n", sep = "")
  cat("Return type:  ", return_type, "\n", sep = "")
  cat("Return thr:   ", return_threshold, "\n", sep = "")
  cat("Signal thr:   ", signal_threshold, " (acted if |signal_matrix| > thr)\n", sep = "")
  cat("Hit test:     ", ifelse(one_sided_hit, "one-sided (p>0.5)", "two-sided (p≠0.5)"), "\n", sep = "")
  cat("==========================================\n")
  
  # --- Merge into a single analysis table
  data <- signal_position %>%
    inner_join(returns,    by = "ticker") %>%
    inner_join(signal_matrix, by = "ticker") %>%
    transmute(
      ticker       = ticker,
      signal_position = as.numeric(signal_position),  # intended discrete {-1,0,1}
      signal_matrix   = as.numeric(signal_matrix),  # continuous signal_matrix
      ret_raw      = as.numeric(returns)
    ) %>%
    filter(is.finite(signal_position), is.finite(signal_matrix), is.finite(ret_raw)) %>%
    mutate(
      # Convert to SIMPLE returns for all downstream evaluation
      ret = standardize_returns(ret_raw, return_type = return_type)
    )
  
  n_obs <- nrow(data)
  cat("Observations after merge/filter: ", n_obs, "\n", sep = "")
  if (n_obs == 0) stop("No usable observations after merging/filtering.")
  
  # --- Discrete signal groups (drives confusion matrices + group stats)
  data <- data %>%
    mutate(
      signal_group = case_when(
        signal_position ==  1 ~ "Signal_Plus1",
        signal_position == -1 ~ "Signal_Minus1",
        TRUE               ~ "Signal_Zero"
      ),
      signal_group  = factor(signal_group, levels = c("Signal_Plus1", "Signal_Zero", "Signal_Minus1")),
      signal_action = case_when(
        signal_group == "Signal_Plus1"  ~  1,
        signal_group == "Signal_Minus1" ~ -1,
        TRUE                            ~  0
      )
    )
  
  # --- Group summary (console)
  group_counts <- data %>%
    group_by(signal_group) %>%
    summarise(
      n        = dplyr::n(),
      mean_ret = mean(ret, na.rm = TRUE),
      sd_ret   = sd(ret, na.rm = TRUE),
      .groups  = "drop"
    )
  print(group_counts)
  
  # ===========================================================================
  # A) Distribution summaries (CSV outputs)
  # ===========================================================================
  
  moments_list <- lapply(levels(data$signal_group), function(g) {
    calculate_sample_moments(data$ret[data$signal_group == g])
  })
  names(moments_list) <- levels(data$signal_group)
  moments_tbl <- as.data.frame(do.call(cbind, moments_list))
  rownames(moments_tbl) <- c("Mean", "Sigma", "Skew", "Kurt", "N")
  
  write.csv(
    moments_tbl,
    file = file.path(output_dir, paste0(safe_slug(asset_name), "_", safe_slug(signal_name), "_moments_tbl.csv"))
  )
  
  quant_list <- lapply(levels(data$signal_group), function(g) {
    calculate_sample_quantiles(data$ret[data$signal_group == g])
  })
  names(quant_list) <- levels(data$signal_group)
  quant_tbl <- as.data.frame(do.call(cbind, quant_list))
  
  write.csv(
    quant_tbl,
    file = file.path(output_dir, paste0(safe_slug(asset_name), "_", safe_slug(signal_name), "_quantiles_tbl.csv"))
  )
  
  # --- Boxplot (PNG)
  boxplot_gg <- ggplot(data, aes(signal_group, ret, fill = signal_group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.2, alpha = 0.25, size = 1.2) +
    labs(
      title = paste("Return Distribution by Signal Group -", signal_name),
      subtitle = paste("Asset:", asset_name),
      x = "Signal Group",
      y = "Returns (simple space)"
    ) +
    theme_minimal() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(
    file.path(output_dir, paste0(safe_slug(asset_name), "_", safe_slug(signal_name), "_boxplots.png")),
    boxplot_gg, width = 10, height = 6, dpi = 300
  )
  
  # --- Density tests (CSV)
  comps <- list(
    "+1 vs -1" = c("Signal_Plus1", "Signal_Minus1"),
    "+1 vs 0"  = c("Signal_Plus1", "Signal_Zero"),
    "0 vs -1"  = c("Signal_Zero",  "Signal_Minus1")
  )
  
  dens_tbl <- do.call(rbind, lapply(names(comps), function(nm) {
    g <- comps[[nm]]
    out <- compare_density_distributions(
      data$ret[data$signal_group == g[1]],
      data$ret[data$signal_group == g[2]]
    )
    data.frame(pair = nm,
               diff_mean_pval = out["diff_mean_pval"],
               ks_pval = out["ks_pval"])
  }))
  
  write.csv(
    dens_tbl,
    file = file.path(output_dir, paste0(safe_slug(asset_name), "_", safe_slug(signal_name), "_densitytests_tbl.csv")),
    row.names = FALSE
  )
  
  # ===========================================================================
  # B) Information coefficient (CSV)
  # ===========================================================================
  
  ic_all   <- calculate_information_coefficient(data$signal_position, data$ret)
  ic_acted <- calculate_ic_acted_only(data$signal_position, data$ret)
  
  ic_tbl <- data.frame(
    Metric  = c("IC (All)", "IC (Acted Only)"),
    Value   = c(ic_all$ic, ic_acted$ic),
    t_stat  = c(ic_all$t_stat, ic_acted$t_stat),
    p_value = c(ic_all$p_value, ic_acted$p_value),
    N       = c(ic_all$n, ic_acted$n_acted),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    ic_tbl,
    file = file.path(output_dir, paste0(safe_slug(asset_name), "_", safe_slug(signal_name), "_ic_tbl.csv")),
    row.names = FALSE
  )
  
  # ===========================================================================
  # C) Confusion matrices + ROC/AUC (CSV + PNG outputs)
  # ===========================================================================
  
  # --- Confusion matrices
  action_conf <- create_action_confusion_matrix(
    signal_action   = data$signal_action,
    ret             = data$ret,
    return_threshold = return_threshold,
    outcome_mode    = outcome_mode
  )
  
  longshort_conf <- create_longshort_confusion_matrix(
    signal_action   = data$signal_action,
    ret             = data$ret,
    return_threshold = return_threshold,
    outcome_mode    = outcome_mode
  )
  
  # --- ROC / AUC
  # Pass RAW signed signal_matrix in both calls; ROC function handles abs() internally for action mode.
  action_roc <- create_roc_analysis(
    signal_matrix            = data$signal_matrix,
    ret              = data$ret,
    roc_mode         = "action",
    return_threshold = return_threshold,
    outcome_mode     = outcome_mode,
    plot_title       = paste("ROC Curve - Action vs No Action -", signal_name),
    acted_only       = FALSE
  )
  
  longshort_roc <- create_roc_analysis(
    signal_matrix            = data$signal_matrix,
    ret              = data$ret,
    roc_mode         = "longshort",
    return_threshold = return_threshold,
    outcome_mode     = outcome_mode,
    plot_title       = paste("ROC Curve - Long vs Short -", signal_name),
    acted_only       = TRUE,
    signal_threshold = signal_threshold
  )
  
  # --- Save confusion matrix tables (tidy)
  action_conf_tbl <- as.data.frame(action_conf$confusion_matrix) %>%
    rename(Count = Freq)
  
  longshort_conf_tbl <- as.data.frame(longshort_conf$confusion_matrix) %>%
    rename(Count = Freq)
  
  write.csv(
    action_conf_tbl,
    file = file.path(output_dir, paste0(safe_slug(asset_name), "_", safe_slug(signal_name), "_action_confusion_tbl.csv")),
    row.names = FALSE
  )
  
  write.csv(
    longshort_conf_tbl,
    file = file.path(output_dir, paste0(safe_slug(asset_name), "_", safe_slug(signal_name), "_longshort_confusion_tbl.csv")),
    row.names = FALSE
  )
  
  # --- Save confusion summary tables (metrics + AUC)
  metrics_order <- c("Accuracy", "Precision", "Recall", "Specificity", "Type1_Error", "Type2_Error")
  
  action_metrics <- action_conf$metrics[c("Accuracy", "Precision", "Recall", "Specificity", "Type1_Error", "Type2_Error")]
  longshort_metrics <- longshort_conf$metrics[c("Accuracy", "Precision", "Recall", "Specificity", "Type1_Error", "Type2_Error")]
  
  action_summary_tbl <- data.frame(
    Metric = c(metrics_order, "AUC"),
    Action_vs_NoAction = c(as.numeric(action_metrics), action_roc$auc),
    stringsAsFactors = FALSE
  )
  
  longshort_summary_tbl <- data.frame(
    Metric = c(metrics_order, "AUC"),
    Long_vs_Short = c(as.numeric(longshort_metrics), longshort_roc$auc),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    action_summary_tbl,
    file = file.path(output_dir, paste0(safe_slug(asset_name), "_", safe_slug(signal_name), "_action_confusionsummary_tbl.csv")),
    row.names = FALSE
  )
  
  write.csv(
    longshort_summary_tbl,
    file = file.path(output_dir, paste0(safe_slug(asset_name), "_", safe_slug(signal_name), "_longshort_confusionsummary_tbl.csv")),
    row.names = FALSE
  )
  
  # --- Save ROC plots (PNG)
  if (!is.null(action_roc$plot)) {
    ggsave(
      file.path(output_dir, paste0(safe_slug(asset_name), "_", safe_slug(signal_name), "_action_roc.png")),
      action_roc$plot, width = 8, height = 6, dpi = 300
    )
  }
  
  if (!is.null(longshort_roc$plot)) {
    ggsave(
      file.path(output_dir, paste0(safe_slug(asset_name), "_", safe_slug(signal_name), "_longshort_roc.png")),
      longshort_roc$plot, width = 8, height = 6, dpi = 300
    )
  }
  
  # ===========================================================================
  # D) Hit rate (CSV)
  # ===========================================================================
  
  hit <- calculate_hit_rate(
    signal_action = data$signal_action,
    ret           = data$ret,
    p0            = 0.5,
    one_sided     = one_sided_hit
  )
  
  hit_tbl <- data.frame(
    Metric  = c("Hit Rate (Overall)", "Hit Rate (Acted Only)", "Coverage"),
    Value   = c(hit$hit_rate_overall, hit$hit_rate_acted, hit$coverage),
    z_stat  = c(hit$z_overall, hit$z_acted, NA),
    p_value = c(hit$p_overall, hit$p_acted, NA),
    N       = c(hit$overall_n, hit$acted_n, hit$overall_n),
    Hits    = c(hit$overall_hits, hit$acted_hits, NA),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    hit_tbl,
    file = file.path(output_dir, paste0(safe_slug(asset_name), "_", safe_slug(signal_name), "_hit_tbl.csv")),
    row.names = FALSE
  )
  
  # ===========================================================================
  # Return outputs
  # ===========================================================================
  
  list(
    data                        = data,
    group_counts                = group_counts,
    
    moments_tbl                 = moments_tbl,
    quantiles_tbl               = quant_tbl,
    density_tests_tbl           = dens_tbl,
    
    ic_tbl                      = ic_tbl,
    
    action_confusion            = action_conf,
    longshort_confusion         = longshort_conf,
    action_confusion_tbl        = action_conf_tbl,
    longshort_confusion_tbl     = longshort_conf_tbl,
    action_confusion_summary_tbl   = action_summary_tbl,
    longshort_confusion_summary_tbl = longshort_summary_tbl,
    
    action_roc                  = action_roc,
    longshort_roc               = longshort_roc,
    
    hit_tbl                     = hit_tbl,
    output_dir                  = output_dir
  )
  
}



# -----------------------------
# Example: run the full pipeline
# -----------------------------
# 
# set.seed(123)
# 
# # 1) Create a toy cross-section of 300 assets
# n <- 300
# tickers <- paste0("A", sprintf("%03d", 1:n))
# 
# # continuous "momentum" signal (raw)
# signal_matrix_vec <- rnorm(n, mean = 0, sd = 1)
# 
# # discretize into {-1,0,+1} using a threshold tau
# tau <- 0.05
# signal_position_vec <- ifelse(signal_matrix_vec >  tau,  1,
#                           ifelse(signal_matrix_vec < -tau, -1, 0))
# 
# # create test-period returns that are mildly related to the signal
# eps <- rnorm(n, 0, 0.05)
# ret_simple <- 0.01 * signal_matrix_vec + eps   # simple returns
# 
# # 2) Build the three required input data frames
# signal_position <- data.frame(
#   ticker = tickers,
#   signal_position = signal_position_vec
# )
# 
# signal_matrix <- data.frame(
#   ticker = tickers,
#   signal_matrix = signal_matrix_vec
# )
# 
# returns <- data.frame(
#   ticker = tickers,
#   returns = ret_simple   # can also be log returns if set return_type="log"
# )
# 
# Meta <- list(
#   assetname = "SP500",
#   benchmarkname = "SP500",
#   signalname = "Momentum"
# )
# 
# # 3) Source functions file (or paste them above) then run:
# res <- SignalTestingStatic(
#   signal_position = signal_position,
#   signal_matrix = signal_matrix,
#   returns = returns,
#   Meta = Meta,
#   signal_threshold = tau,          # acted iff |signal_matrix| > tau
#   return_threshold = 0,            # for outcome_mode="direction" it doesn't matter
#   outcome_mode = "direction",
#   return_type = "simple",
#   one_sided_hit = FALSE,
#   seed = 123,
#   output_dir = "../Output/SignalTestingStatic"
# )
# 
# # 4) Inspect outputs
# names(res)
# 
# # Confusion matrices
# res$action_confusion$confusion_matrix
# res$longshort_confusion$confusion_matrix
# 
# # AUCs
# res$action_roc$auc
# res$longshort_roc$auc
# 
# # Show ROC plots in RStudio
# print(res$action_roc$plot)
# print(res$longshort_roc$plot)

