# PerformanceStats.R
# Purpose: Compute performance + risk diagnostics for an asset relative to a benchmark.
# Note: If computing sample statistics for a portfolio, the portfolio must be created before feeding into this function.
#
#
# Inputs:
#   price_df: data.frame with columns [Date, Asset, Benchmark] (levels)
#   meta:     list(asset_name = "...", benchmark_name = "...")
#   rf_annual: annual risk-free rate (decimal, e.g. 0.02). Default 0.
#   output_dir: directory for saved outputs; created if missing.
#
# Output:
#   A list of plots/tables; also saves .png and .csv files under output_dir.
#
# Authors:
#   Mike Aguilar (mike.aguilar@duke.edu), Ziming Huang

source("../Supporting/PackageLoads.R")

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

# Validate + standardize column names
standardize_price_df <- function(price_df) {
  stopifnot(is.data.frame(price_df))
  stopifnot(ncol(price_df) >= 3)
  
  out <- price_df[, 1:3]
  names(out) <- c("Date", "Asset", "Benchmark")
  
  out$Date <- as.Date(out$Date)
  out <- out %>% dplyr::arrange(Date)
  
  out
}

ensure_output_dir <- function(output_dir) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }
  invisible(output_dir)
}

# Simple returns from price levels
compute_simple_returns <- function(df_levels) {
  df_levels %>%
    dplyr::mutate(
      Asset_Return     = Asset / dplyr::lag(Asset) - 1,
      Benchmark_Return = Benchmark / dplyr::lag(Benchmark) - 1
    ) %>%
    tidyr::drop_na(Asset_Return, Benchmark_Return)
}

# Aggregate to frequency and (optionally) annualize the aggregated return
# freq: "month", "quarter", "year", or "day"
aggregate_returns <- function(returns, dates, freq = "month", annualize = TRUE) {
  df <- data.frame(Date = as.Date(dates), Return = returns)
  df$Period <- lubridate::floor_date(df$Date, unit = freq)
  
  agg <- df %>%
    dplyr::group_by(Period) %>%
    dplyr::summarise(
      Return = prod(1 + Return) - 1,
      .groups = "drop"
    ) %>%
    dplyr::filter(is.finite(Return))
  
  if (!annualize) return(agg$Return)
  
  ann_factor <- dplyr::case_when(
    freq == "day"     ~ 252,
    freq == "month"   ~ 12,
    freq == "quarter" ~ 4,
    freq == "year"    ~ 1,
    TRUE              ~ 252
  )
  
  (1 + agg$Return)^ann_factor - 1
}

# Omega ratio at threshold (per-period threshold, NOT annual)
omega_ratio <- function(returns, threshold = 0) {
  returns <- returns[is.finite(returns)]
  if (length(returns) == 0) return(NA_real_)
  
  gains  <- returns[returns >  threshold] - threshold
  losses <- threshold - returns[returns <= threshold]
  
  sg <- sum(gains)
  sl <- sum(losses)
  
  if (sl == 0) return(NA_real_)
  sg / sl
}

# Risk stats assuming DAILY simple returns by default
risk_stats_daily <- function(returns_daily, rf_annual = 0) {
  r <- returns_daily[is.finite(returns_daily)]
  if (length(r) == 0) return(rep(NA_real_, 13))
  
  rf_daily <- rf_annual / 252
  
  ann_mean <- mean(r) * 252
  ann_sd   <- sd(r) * sqrt(252)
  
  sharpe <- ifelse(ann_sd > 0, (ann_mean - rf_annual) / ann_sd, NA_real_)
  
  downside <- r[r < 0]
  downside_dev <- ifelse(length(downside) > 1, sd(downside) * sqrt(252), NA_real_)
  sortino <- ifelse(is.finite(downside_dev) && downside_dev > 0,
                    (ann_mean - rf_annual) / downside_dev, NA_real_)
  
  omega <- omega_ratio(r, threshold = rf_daily)
  
  skew <- moments::skewness(r, na.rm = TRUE)
  kurt <- moments::kurtosis(r, na.rm = TRUE) - 3  # excess
  
  quants <- stats::quantile(r, probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = TRUE)
  
  c(Sigma = ann_sd, Sharpe = sharpe, Sortino = sortino, Omega = omega,
    Skew = skew, Kurt = kurt, quants)
}

save_plot <- function(p, path, w = 12, h = 6) {
  ggplot2::ggsave(filename = path, plot = p, width = w, height = h)
  invisible(path)
}

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------

PerformanceStats <- function(price_df,
                             meta,
                             rf_annual = 0,
                             output_dir = "PerformanceStat") {
  
  stopifnot(is.list(meta))
  stopifnot(all(c("asset_name", "benchmark_name") %in% names(meta)))
  
  ensure_output_dir(output_dir)
  
  asset_name     <- meta$asset_name
  benchmark_name <- meta$benchmark_name
  
  message("==========================================")
  message("Performance analysis starts")
  message("Asset:     ", asset_name)
  message("Benchmark: ", benchmark_name)
  message("RF (annual): ", rf_annual)
  message("Output dir:  ", output_dir)
  message("==========================================")
  
  # 0) Standardize + compute returns
  df_levels  <- standardize_price_df(price_df)
  df_returns <- compute_simple_returns(df_levels)
  
  asset_xts <- xts::xts(df_returns$Asset_Return, order.by = df_returns$Date)
  bench_xts <- xts::xts(df_returns$Benchmark_Return, order.by = df_returns$Date)
  
  # 1) Price plot + MAs (levels)
  message("Step 1/9: Price plot...")
  df_plot_price <- df_levels %>%
    dplyr::mutate(
      Asset_MA50  = zoo::rollmean(Asset, 50,  fill = NA, align = "right"),
      Asset_MA250 = zoo::rollmean(Asset, 250, fill = NA, align = "right")
    )
  
  price_plot <- ggplot2::ggplot(df_plot_price, ggplot2::aes(x = Date)) +
    ggplot2::geom_line(ggplot2::aes(y = Asset, color = asset_name),
                       linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = Asset_MA50, color = paste(asset_name, "50-day MA")),
                       linetype = "dashed", linewidth = 0.6, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = Asset_MA250, color = paste(asset_name, "250-day MA")),
                       linetype = "dashed", linewidth = 0.6, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = Benchmark, color = benchmark_name),
                       linewidth = 0.8, alpha = 0.7, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c("black", "blue", "red", "gray")) +
    ggplot2::labs(
      title = paste(asset_name, "vs", benchmark_name, "- Price Series"),
      x = "Date", y = "Value", color = "Series"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  
  save_plot(price_plot, file.path(output_dir, paste0(asset_name, "_priceplot.png")))
  
  # 2) Rolling 5-day returns (daily data, 5-day compounded)
  message("Step 2/9: Return plot...")
  df_plot_ret <- df_returns %>%
    dplyr::mutate(
      Asset_5D = zoo::rollapply(Asset_Return, 5, function(x) prod(1 + x) - 1,
                                fill = NA, align = "right"),
      Benchmark_5D = zoo::rollapply(Benchmark_Return, 5, function(x) prod(1 + x) - 1,
                                    fill = NA, align = "right")
    )
  
  return_plot <- ggplot2::ggplot(df_plot_ret, ggplot2::aes(x = Date)) +
    ggplot2::geom_line(ggplot2::aes(y = Asset_5D, color = asset_name),
                       linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = Benchmark_5D, color = benchmark_name),
                       linewidth = 0.8, alpha = 0.7, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c("blue", "gray")) +
    ggplot2::labs(title = "Rolling 5-Day Returns", x = "Date", y = "5-Day Return", color = "Series") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  
  save_plot(return_plot, file.path(output_dir, paste0(asset_name, "_returnsplot.png")))
  
  # 3) Rolling 6-month sigma (126 trading days), annualized
  message("Step 3/9: Sigma plot...")
  df_plot_sig <- df_returns %>%
    dplyr::mutate(
      Asset_Sigma = zoo::rollapply(Asset_Return, 126, function(x) sd(x) * sqrt(252),
                                   fill = NA, align = "right"),
      Benchmark_Sigma = zoo::rollapply(Benchmark_Return, 126, function(x) sd(x) * sqrt(252),
                                       fill = NA, align = "right")
    )
  
  sigma_plot <- ggplot2::ggplot(df_plot_sig, ggplot2::aes(x = Date)) +
    ggplot2::geom_line(ggplot2::aes(y = Asset_Sigma, color = asset_name),
                       linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = Benchmark_Sigma, color = benchmark_name),
                       linewidth = 0.8, alpha = 0.7, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c("blue", "gray")) +
    ggplot2::labs(
      title = "Rolling 6-Month Daily Sigma (Annualized)",
      x = "Date", y = "Annualized Volatility", color = "Series"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  
  save_plot(sigma_plot, file.path(output_dir, paste0(asset_name, "_sigmaplot.png")))
  
  # 4) Sample stats table (levels-based total return + CAGR over trailing windows)
  message("Step 4/9: Sample stats...")
  windows <- c("5D" = 5, "1M" = 21, "3M" = 63, "6M" = 126, "1Y" = 252, "3Y" = 756, "5Y" = 1260)
  
  sample_stats_data <- lapply(names(windows), function(lbl) {
    days <- windows[[lbl]]
    if (nrow(df_levels) < days) return(rep(NA_real_, 4))
    
    sub <- tail(df_levels, days)
    start_a <- sub$Asset[1];     end_a <- sub$Asset[nrow(sub)]
    start_b <- sub$Benchmark[1]; end_b <- sub$Benchmark[nrow(sub)]
    
    total_a <- (end_a / start_a - 1) * 100
    total_b <- (end_b / start_b - 1) * 100
    
    years <- days / 252
    cagr_a <- ((end_a / start_a)^(1 / years) - 1) * 100
    cagr_b <- ((end_b / start_b)^(1 / years) - 1) * 100
    
    c(total_a, cagr_a, total_b, cagr_b)
  })
  names(sample_stats_data) <- names(windows)
  
  sample_stats_tbl <- as.data.frame(do.call(cbind, sample_stats_data))
  rownames(sample_stats_tbl) <- c("Asset Total %", "Asset CAGR %", "Benchmark Total %", "Benchmark CAGR %")
  
  utils::write.csv(sample_stats_tbl,
                   file = file.path(output_dir, paste0(asset_name, "_samplestats_tbl.csv")))
  
  # Risk-return table (daily returns; annualized mean/sd; higher moments; quantiles)
  message("  Risk-return table...")
  asset_stats <- risk_stats_daily(df_returns$Asset_Return, rf_annual = rf_annual)
  bench_stats <- risk_stats_daily(df_returns$Benchmark_Return, rf_annual = rf_annual)
  
  risk_return_tbl <- data.frame(
    Metric    = c("Sigma", "Sharpe", "Sortino", "Omega", "Skew", "Excess Kurtosis",
                  "Min", "Q.05", "Q.25", "Q.50", "Q.75", "Q.95", "Max"),
    Asset     = as.numeric(asset_stats),
    Benchmark = as.numeric(bench_stats)
  )
  
  utils::write.csv(risk_return_tbl,
                   file = file.path(output_dir, paste0(asset_name, "_riskreturn_tbl.csv")),
                   row.names = FALSE)
  
  # Boxplot data: annualized returns by frequency (consistent definition)
  message("  Return distribution boxplot...")
  daily_ann     <- (1 + df_returns$Asset_Return)^252 - 1
  monthly_ann   <- aggregate_returns(df_returns$Asset_Return, df_returns$Date, "month",   annualize = TRUE)
  quarterly_ann <- aggregate_returns(df_returns$Asset_Return, df_returns$Date, "quarter", annualize = TRUE)
  annual_ann    <- aggregate_returns(df_returns$Asset_Return, df_returns$Date, "year",    annualize = TRUE)
  
  box_data <- data.frame(
    Return = c(daily_ann, monthly_ann, quarterly_ann, annual_ann),
    Frequency = factor(
      rep(c("Daily", "Monthly", "Quarterly", "Annual"),
          c(length(daily_ann), length(monthly_ann), length(quarterly_ann), length(annual_ann))),
      levels = c("Daily", "Monthly", "Quarterly", "Annual")
    )
  )
  
  return_boxplot <- ggplot2::ggplot(box_data, ggplot2::aes(x = Frequency, y = Return)) +
    ggplot2::geom_boxplot(outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.2, alpha = 0.25, size = 0.8) +
    ggplot2::labs(title = "Annualized Return Distributions by Frequency",
                  x = "Frequency", y = "Annualized Return") +
    ggplot2::theme_minimal()
  
  save_plot(return_boxplot, file.path(output_dir, paste0(asset_name, "_returnboxplot.png")), w = 10, h = 6)
  
  # 5) Relative performance (daily return diff)
  message("Step 5/9: Relative performance...")
  df_rel <- df_returns %>% dplyr::mutate(Relative_Perf = Asset_Return - Benchmark_Return)
  
  rel_perf_plot <- ggplot2::ggplot(df_rel, ggplot2::aes(x = Date, y = Relative_Perf)) +
    ggplot2::geom_bar(stat = "identity", width = 1, fill = "steelblue", alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    ggplot2::labs(title = "Daily Relative Performance (Asset - Benchmark)",
                  x = "Date", y = "Return Difference") +
    ggplot2::theme_minimal()
  
  save_plot(rel_perf_plot, file.path(output_dir, paste0(asset_name, "_relperf.png")))
  
  # 6) Seasonality (monthly heatmap last 5 years; weekday x month heatmap)
  message("Step 6/9: Seasonality...")
  df_seas <- df_returns %>%
    dplyr::mutate(
      Year    = lubridate::year(Date),
      Month   = lubridate::month(Date, label = TRUE, abbr = TRUE),
      Weekday = lubridate::wday(Date, label = TRUE, abbr = TRUE, week_start = 1)
    )
  
  recent_years <- seq(max(df_seas$Year), max(df_seas$Year) - 4, by = -1)
  
  monthly_data <- df_seas %>%
    dplyr::filter(Year %in% recent_years) %>%
    dplyr::group_by(Year, Month) %>%
    dplyr::summarise(Avg_Return = mean(Asset_Return, na.rm = TRUE) * 100, .groups = "drop")
  
  monthly_heatmap <- ggplot2::ggplot(monthly_data,
                                     ggplot2::aes(x = factor(Year), y = Month, fill = Avg_Return)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "red", mid = "white", high = "green",
                                  midpoint = 0, name = "Avg Return (%)") +
    ggplot2::labs(title = "Monthly Average Returns (Last 5 Years)",
                  x = "Year", y = "Month") +
    ggplot2::theme_minimal()
  
  save_plot(monthly_heatmap, file.path(output_dir, paste0(asset_name, "_monthlyheatmap.png")), w = 8, h = 6)
  
  daily_data <- df_seas %>%
    dplyr::group_by(Month, Weekday) %>%
    dplyr::summarise(Avg_Return = mean(Asset_Return, na.rm = TRUE) * 100, .groups = "drop")
  
  daily_heatmap <- ggplot2::ggplot(daily_data, ggplot2::aes(x = Month, y = Weekday, fill = Avg_Return)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::scale_fill_gradient2(low = "red", mid = "white", high = "green",
                                  midpoint = 0, name = "Avg Return (%)") +
    ggplot2::labs(title = "Average Returns by Day of Week and Month",
                  x = "Month", y = "Day of Week") +
    ggplot2::theme_minimal()
  
  save_plot(daily_heatmap, file.path(output_dir, paste0(asset_name, "_dailyheatmap.png")), w = 10, h = 6)
  
  # 7) Comovement (rolling correlation)
  message("Step 7/9: Comovement...")
  df_corr <- df_returns %>%
    dplyr::mutate(
      Corr_3M = zoo::rollapply(cbind(Asset_Return, Benchmark_Return), 63,
                               function(x) cor(x[, 1], x[, 2], use = "complete.obs"),
                               by.column = FALSE, fill = NA, align = "right"),
      Corr_1Y = zoo::rollapply(cbind(Asset_Return, Benchmark_Return), 252,
                               function(x) cor(x[, 1], x[, 2], use = "complete.obs"),
                               by.column = FALSE, fill = NA, align = "right")
    )
  
  corr_plot <- ggplot2::ggplot(df_corr, ggplot2::aes(x = Date)) +
    ggplot2::geom_line(ggplot2::aes(y = Corr_3M, color = "3-Month"), linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_line(ggplot2::aes(y = Corr_1Y, color = "1-Year"),  linewidth = 0.8, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    ggplot2::scale_color_manual(values = c("3-Month" = "blue", "1-Year" = "red")) +
    ggplot2::labs(title = "Rolling Correlation with Benchmark", x = "Date", y = "Correlation", color = "Window") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
  
  save_plot(corr_plot, file.path(output_dir, paste0(asset_name, "_benchmarkcor.png")))
  
  # 8) Drawdowns (underwater plot + table)
  message("Step 8/9: Drawdowns...")
  df_dd <- df_returns %>%
    dplyr::mutate(
      Asset_Cum      = cumprod(1 + Asset_Return),
      Benchmark_Cum  = cumprod(1 + Benchmark_Return),
      Asset_Max      = cummax(Asset_Cum),
      Benchmark_Max  = cummax(Benchmark_Cum),
      Asset_Drawdown = (Asset_Cum - Asset_Max) / Asset_Max,
      Bench_Drawdown = (Benchmark_Cum - Benchmark_Max) / Benchmark_Max
    )
  
  underwater_plot <- ggplot2::ggplot(df_dd, ggplot2::aes(x = Date)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = Asset_Drawdown, ymax = 0), fill = "red", alpha = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = Asset_Drawdown), color = "darkred", linewidth = 0.8) +
    ggplot2::scale_y_continuous(labels = scales::percent_format()) +
    ggplot2::labs(title = "Underwater Plot - Asset Drawdowns", x = "Date", y = "Drawdown") +
    ggplot2::theme_minimal()
  
  save_plot(underwater_plot, file.path(output_dir, paste0(asset_name, "_underwater.png")))
  
  drawdowns_tbl <- PerformanceAnalytics::table.Drawdowns(asset_xts, top = 5)
  utils::write.csv(drawdowns_tbl, file.path(output_dir, paste0(asset_name, "_drawdowns_tbl.csv")))
  
  # 9) Persistence (ACF/CCF)
  message("Step 9/9: Persistence (ACF/CCF)...")
  acf_daily <- stats::acf(df_returns$Asset_Return, lag.max = 50, plot = FALSE)
  acf_daily_df <- data.frame(Lag = as.numeric(acf_daily$lag), ACF = as.numeric(acf_daily$acf))
  
  weekly_asset <- df_returns %>%
    dplyr::mutate(Week = lubridate::floor_date(Date, "week")) %>%
    dplyr::group_by(Week) %>%
    dplyr::summarise(Return = prod(1 + Asset_Return) - 1, .groups = "drop")
  
  weekly_bench <- df_returns %>%
    dplyr::mutate(Week = lubridate::floor_date(Date, "week")) %>%
    dplyr::group_by(Week) %>%
    dplyr::summarise(Return = prod(1 + Benchmark_Return) - 1, .groups = "drop")
  
  acf_weekly <- stats::acf(weekly_asset$Return, lag.max = 20, plot = FALSE)
  acf_weekly_df <- data.frame(Lag = as.numeric(acf_weekly$lag), ACF = as.numeric(acf_weekly$acf))
  
  p1 <- ggplot2::ggplot(acf_daily_df, ggplot2::aes(x = Lag, y = ACF)) +
    ggplot2::geom_bar(stat = "identity", width = 0.5, fill = "steelblue") +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::labs(title = "ACF - Daily Returns", x = "Lag", y = "ACF") +
    ggplot2::theme_minimal()
  
  p2 <- ggplot2::ggplot(acf_weekly_df, ggplot2::aes(x = Lag, y = ACF)) +
    ggplot2::geom_bar(stat = "identity", width = 0.5, fill = "steelblue") +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::labs(title = "ACF - Weekly Returns", x = "Lag", y = "ACF") +
    ggplot2::theme_minimal()
  
  acf_combined <- gridExtra::grid.arrange(p1, p2, nrow = 2)
  ggplot2::ggsave(file.path(output_dir, paste0(asset_name, "_acf.png")),
                  acf_combined, width = 10, height = 8)
  
  ccf_daily <- stats::ccf(df_returns$Asset_Return, df_returns$Benchmark_Return,
                          lag.max = 50, plot = FALSE)
  ccf_daily_df <- data.frame(Lag = as.numeric(ccf_daily$lag), CCF = as.numeric(ccf_daily$acf))
  
  ccf_daily_plot <- ggplot2::ggplot(ccf_daily_df, ggplot2::aes(x = Lag, y = CCF)) +
    ggplot2::geom_bar(stat = "identity", width = 0.5, fill = "steelblue") +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::labs(title = "CCF - Daily Returns (Asset vs Benchmark)", x = "Lag", y = "Cross-Correlation") +
    ggplot2::theme_minimal()
  
  save_plot(ccf_daily_plot, file.path(output_dir, paste0(asset_name, "_ccf_daily.png")), w = 10, h = 6)
  
  ccf_weekly <- stats::ccf(weekly_asset$Return, weekly_bench$Return, lag.max = 20, plot = FALSE)
  ccf_weekly_df <- data.frame(Lag = as.numeric(ccf_weekly$lag), CCF = as.numeric(ccf_weekly$acf))
  
  ccf_weekly_plot <- ggplot2::ggplot(ccf_weekly_df, ggplot2::aes(x = Lag, y = CCF)) +
    ggplot2::geom_bar(stat = "identity", width = 0.5, fill = "steelblue") +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::labs(title = "CCF - Weekly Returns (Asset vs Benchmark)", x = "Lag", y = "Cross-Correlation") +
    ggplot2::theme_minimal()
  
  save_plot(ccf_weekly_plot, file.path(output_dir, paste0(asset_name, "_ccf_weekly.png")), w = 10, h = 6)
  
  message("==========================================")
  message("Performance analysis complete!")
  message("Saved outputs under: ", output_dir)
  message("==========================================")
  
  # Return outputs (keep naming consistent)
  list(
    price_plot       = price_plot,
    return_plot      = return_plot,
    sigma_plot       = sigma_plot,
    sample_stats_tbl = sample_stats_tbl,
    risk_return_tbl  = risk_return_tbl,
    return_boxplot   = return_boxplot,
    rel_perf_plot    = rel_perf_plot,
    monthly_heatmap  = monthly_heatmap,
    daily_heatmap    = daily_heatmap,
    corr_plot        = corr_plot,
    underwater_plot  = underwater_plot,
    drawdowns_tbl    = drawdowns_tbl,
    acf_combined     = acf_combined,
    ccf_daily_plot   = ccf_daily_plot,
    ccf_weekly_plot  = ccf_weekly_plot
  )
}

# Example usage (fixed lengths)
# set.seed(124)
# n <- 2000
# Price <- data.frame(
#   Date      = seq.Date(from = as.Date("2018-01-01"), by = "day", length.out = n),
#   Asset     = cumprod(1 + rnorm(n, 0.0005, 0.02)),
#   Benchmark = cumprod(1 + rnorm(n, 0.0003, 0.015))
# )
# 
# Meta <- list(asset_name = "MyPortfolio", benchmark_name = "SP500")
# results <- PerformanceStats(Price, Meta, rf_annual = 0, output_dir = "../Output/PerformanceStat")
