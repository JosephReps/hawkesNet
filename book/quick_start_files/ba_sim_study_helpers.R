# quick_start_files/ba_sim_study_helpers.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(readr)
  library(ggplot2)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
as_num <- function(x) suppressWarnings(as.numeric(x))

read_ba_rep <- function(file) {
  obj <- tryCatch(readRDS(file), error = function(e) NULL)
  if (is.null(obj) || !is.list(obj)) return(NULL)

  params_true <- obj$params_true %||% list()

  tibble(
    file = file,
    run_id = basename(dirname(file)),
    rep_id = as.integer(obj$rep_id %||% NA),
    seed = as_num(obj$seed %||% NA),
    T_end = as_num(obj$T_end %||% NA),
    status = as.character(obj$status %||% NA),

    sim_seconds = as_num(obj$sim_seconds %||% NA),
    fit_seconds = as_num(obj$fit_seconds %||% NA),

    n_events = as.integer(obj$n_events %||% NA),
    n_nodes  = as.integer(obj$n_nodes %||% NA),
    n_edges  = as.integer(obj$n_edges %||% NA),

    mean_degree = as_num(obj$mean_degree %||% NA),
    max_degree  = as_num(obj$max_degree %||% NA),
    density     = as_num(obj$density %||% NA),

    dt_mean   = as_num(obj$dt_mean %||% NA),
    dt_median = as_num(obj$dt_median %||% NA),
    lambda_mean = as_num(obj$lambda_mean %||% NA),
    lambda_max  = as_num(obj$lambda_max %||% NA),

    convergence = as.integer(obj$convergence %||% NA),
    converged   = isTRUE(obj$converged %||% FALSE),
    loglik      = as_num(obj$loglik %||% NA),
    fit_message = as.character(obj$fit_message %||% NA),

    mu_true         = as_num(params_true$mu %||% NA),
    K_true          = as_num(params_true$K %||% NA),
    beta_true       = as_num(params_true$beta %||% NA),
    beta_edges_true = as_num(params_true$beta_edges %||% NA),

    mu_hat         = as_num(obj$mu_hat %||% NA),
    K_hat          = as_num(obj$K_hat %||% NA),
    beta_hat       = as_num(obj$beta_hat %||% NA),
    beta_edges_hat = as_num(obj$beta_edges_hat %||% NA)
  ) %>%
    mutate(
      ok_status = is.na(status) | status %in% c("ok", "OK", "success"),
      keep = ok_status & converged & !is.na(convergence) & convergence == 0
    ) %>%
    mutate(
      ok_status = is.na(status) | status %in% c("ok", "OK", "success"),
      ok_conv   = converged & !is.na(convergence) & convergence == 0,

      # cap filter (your request)
      ok_caps   = (is.na(K_hat) | K_hat <= 20) & (is.na(beta_hat) | beta_hat <= 20),

      keep = ok_status & ok_conv & ok_caps,

      drop_reason = case_when(
        !(ok_status) ~ "bad_status",
        !(ok_conv)   ~ "not_converged",
        !(ok_caps)   ~ "cap_exceeded",
        TRUE         ~ NA_character_
      )
    )
}

# Choose latest ba_run_* folder (by folder mtime), unless run_id specified
find_ba_run_dir <- function(base_dir = "sim_study_results/ba", run_id = NULL) {
  if (!dir.exists(base_dir)) stop("Directory not found: ", base_dir)

  if (!is.null(run_id)) {
    cand <- file.path(base_dir, run_id)
    if (!dir.exists(cand)) stop("Run dir not found: ", cand)
    return(cand)
  }

  runs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  runs <- runs[grepl("ba_run_", basename(runs))]
  if (length(runs) == 0) stop("No ba_run_* folders found under: ", base_dir)

  info <- file.info(runs)
  runs[which.max(info$mtime)]
}

load_ba_run <- function(base_dir = "sim_study_results/ba", run_id = NULL) {
  run_dir <- find_ba_run_dir(base_dir, run_id)
  files <- list.files(run_dir, pattern = "^ba_rep_\\d+\\.rds$", full.names = TRUE)
  if (length(files) == 0) stop("No ba_rep_###.rds found in: ", run_dir)

  df <- map_dfr(files, read_ba_rep) %>%
    filter(!is.na(T_end)) %>%
    mutate(T_end = as.numeric(T_end))

  list(run_dir = run_dir, df = df)
}

make_param_long <- function(df) {
  df %>%
    select(run_id, rep_id, T_end, keep,
           mu_hat, K_hat, beta_hat, beta_edges_hat,
           mu_true, K_true, beta_true, beta_edges_true) %>%
    pivot_longer(cols = ends_with("_hat"), names_to = "param_hat", values_to = "estimate") %>%
    mutate(
      param = str_remove(param_hat, "_hat$"),
      true = case_when(
        param == "mu" ~ mu_true,
        param == "K" ~ K_true,
        param == "beta" ~ beta_true,
        param == "beta_edges" ~ beta_edges_true,
        TRUE ~ NA_real_
      ),
      T_end_f = factor(T_end, levels = sort(unique(T_end)))
    )
}

plot_estimates_vs_T <- function(df, keep_only = FALSE, log_scale = FALSE) {

  long <- make_param_long(df) %>%
    filter(!is.na(estimate), !is.na(T_end))

  if (keep_only) {
    long <- long %>% filter(keep)
  }

  # If using log scale, remove non-positive estimates / truths
  if (log_scale) {
    long <- long %>%
      filter(estimate > 0, true > 0)
  }

  p <- ggplot(long, aes(x = T_end_f, y = estimate)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(alpha = keep), width = 0.2, size = 1) +
    geom_hline(aes(yintercept = true), linetype = "dashed") +
    facet_wrap(~param, scales = "free_y") +
    labs(
      title = if (log_scale)
        "Parameter estimates vs T_end (log scale)"
      else
        "Parameter estimates vs T_end",
      subtitle = if (keep_only)
        "Keep-only (converged & ok)"
      else
        "All runs (alpha shows keep)",
      x = "T_end",
      y = if (log_scale) "Estimate (log scale)" else "Estimate",
      alpha = "keep"
    ) +
    theme_minimal()

  if (log_scale) {
    p <- p + scale_y_log10()
  }

  return(p)
}


plot_counts_vs_T <- function(df) {
  long <- df %>%
    select(T_end, keep, n_events, n_nodes, n_edges) %>%
    pivot_longer(cols = c(n_events, n_nodes, n_edges),
                 names_to = "metric", values_to = "value") %>%
    mutate(T_end_f = factor(T_end, levels = sort(unique(T_end))))

  ggplot(long %>% filter(!is.na(T_end), !is.na(value)), aes(x = T_end_f, y = value)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(alpha = keep), width = 0.2, size = 1) +
    facet_wrap(~metric, scales = "free_y") +
    labs(
      title = "Counts vs T_end",
      x = "T_end",
      y = "Value",
      alpha = "keep"
    ) +
    theme_minimal()
}

plot_keep_rate_vs_T <- function(df) {
  sumdf <- df %>%
    group_by(T_end) %>%
    summarise(n = n(), keep_rate = mean(keep, na.rm = TRUE), .groups = "drop") %>%
    arrange(T_end)

  ggplot(sumdf, aes(x = T_end, y = keep_rate)) +
    geom_line() +
    geom_point() +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Keep rate vs T_end", x = "T_end", y = "Keep rate") +
    theme_minimal()
}

plot_rmse_vs_T <- function(df) {
  long <- make_param_long(df) %>%
    filter(keep, !is.na(estimate), !is.na(true), !is.na(T_end))

  sumdf <- long %>%
    group_by(T_end, param) %>%
    summarise(
      bias = mean(estimate - true),
      rmse = sqrt(mean((estimate - true)^2)),
      sd = sd(estimate),
      .groups = "drop"
    )

  ggplot(sumdf, aes(x = T_end, y = rmse)) +
    geom_line() +
    geom_point() +
    facet_wrap(~param, scales = "free_y") +
    labs(title = "RMSE vs T_end (keep-only)", x = "T_end", y = "RMSE") +
    theme_minimal()
}

plot_runtime_vs_T <- function(df) {
  long <- df %>%
    mutate(total_seconds = sim_seconds + fit_seconds) %>%
    select(T_end, keep, sim_seconds, fit_seconds, total_seconds) %>%
    pivot_longer(cols = c(sim_seconds, fit_seconds, total_seconds),
                 names_to = "runtime", values_to = "seconds") %>%
    mutate(T_end_f = factor(T_end, levels = sort(unique(T_end))))

  ggplot(long %>% filter(!is.na(T_end), !is.na(seconds)), aes(x = T_end_f, y = seconds)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(alpha = keep), width = 0.2, size = 1) +
    facet_wrap(~runtime, scales = "free_y") +
    labs(title = "Runtime vs T_end", x = "T_end", y = "Seconds", alpha = "keep") +
    theme_minimal()
}

summary_table <- function(df) {
  # nice little table for the book
  df %>%
    group_by(T_end) %>%
    summarise(
      n = n(),
      keep_rate = mean(keep, na.rm = TRUE),
      n_events_mean = mean(n_events, na.rm = TRUE),
      n_events_sd   = sd(n_events, na.rm = TRUE),
      sim_s_mean = mean(sim_seconds, na.rm = TRUE),
      fit_s_mean = mean(fit_seconds, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(T_end)
}

load_all_ba_runs <- function(base_dir = "sim_study_results/ba") {

  runs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
  runs <- runs[grepl("ba_run_", basename(runs))]

  if (length(runs) == 0) stop("No ba_run_* folders found in: ", base_dir)

  cat("Found", length(runs), "runs:\n")
  print(basename(runs))

  df <- purrr::map_dfr(runs, function(run_dir) {
    files <- list.files(run_dir, pattern = "^ba_rep_\\d+\\.rds$", full.names = TRUE)
    if (length(files) == 0) return(NULL)
    purrr::map_dfr(files, read_ba_rep)
  })

  df <- df %>%
    filter(!is.na(T_end)) %>%
    mutate(T_end = as.numeric(T_end))

  return(df)
}

