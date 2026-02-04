suppressPackageStartupMessages({
  library(parallel)
  library(ernm)
  library(devtools)
})

# --- env helpers (like your old script) ---
as_int <- function(x, default) {
  if (!nzchar(x)) return(as.integer(default))
  suppressWarnings({
    y <- as.integer(x)
    if (is.na(y)) as.integer(default) else y
  })
}

parse_times <- function(x, default = 1:2) {
  if (!nzchar(x)) return(default)
  # allow "1:15" or "c(1,2,3)" or "1,2,3"
  x <- gsub("\\s+", "", x)
  if (grepl(":", x) && !grepl("c\\(", x)) {
    parts <- strsplit(x, ":", fixed = TRUE)[[1]]
    return(seq.int(as.integer(parts[1]), as.integer(parts[2])))
  }
  if (grepl("^c\\(", x)) return(eval(parse(text = x)))
  as.integer(strsplit(x, ",", fixed = TRUE)[[1]])
}

# ---- settings ----
times <- parse_times(Sys.getenv("TIMES", ""), default = 1:15)
reps  <- as_int(Sys.getenv("REPS", ""), 10)

version_label <- Sys.getenv("VERSION_LABEL", "v1")
out_dir <- Sys.getenv("OUT_DIR", "bench_results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# cores: prefer N_CORES (your old style), else slurm cpus
n_cores <- if (nzchar(Sys.getenv("N_CORES", ""))) {
  as_int(Sys.getenv("N_CORES", ""), 1)
} else {
  as_int(Sys.getenv("SLURM_CPUS_PER_TASK", ""), 1)
}
n_cores <- max(1L, n_cores)

# Load package on master exactly once
pkg_path <- normalizePath(getwd())
devtools::load_all(pkg_path, quiet = TRUE)

params_cs_true <- list(
  mu = 5, K = 0.5, beta = 2, beta_edges = 0.5, node_lambda = 4,
  CS_star.2 = -1, CS_star.3 = -3
)

params_cs_init <- list(
  mu = 3, K = 1, beta = 1, beta_edges = 1, node_lambda = 5,
  CS_star.2 = 0, CS_star.3 = 0
)

run_one <- function(T_end, rep_id) {
  set.seed(1000 + 100 * T_end + rep_id)

  sim <- sim_hawkesNet(
    params = params_cs_true,
    T_end = T_end,
    mark_type = "cs",
    formula_rhs = "star(c(2,3))"
  )

  t0 <- proc.time()[["elapsed"]]
  fit <- fit_hawkesNet(
    ev = sim$ev,
    params_init = params_cs_init,
    mark_type = "cs",
    transform = list(CS_star.2 = "none", CS_star.3 = "none"),
    formula_rhs = "star(c(2,3))"
  )
  t1 <- proc.time()[["elapsed"]]

  nn <- if (!is.null(sim$ev$nodes) && is.data.frame(sim$ev$nodes)) nrow(sim$ev$nodes) else NA_integer_
  ne <- if (!is.null(sim$ev$edges) && is.data.frame(sim$ev$edges)) nrow(sim$ev$edges) else NA_integer_

  data.frame(
    T_end = T_end,
    rep = rep_id,
    fit_seconds = (t1 - t0),
    n_nodes = nn,
    n_edges = ne,
    stringsAsFactors = FALSE
  )
}

grid <- expand.grid(T_end = times, rep = seq_len(reps))

# ---- run ----
pt0 <- proc.time()

if (n_cores == 1L || nrow(grid) == 1L) {
  # single-core path
  res_list <- lapply(seq_len(nrow(grid)), function(i) run_one(grid$T_end[i], grid$rep[i]))
} else {
  cl <- parallel::makeCluster(n_cores)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  parallel::clusterSetRNGStream(cl, 1)

  parallel::clusterEvalQ(cl, {
    library(stats)
    library(ernm)
    library(devtools)
    NULL
  })

  parallel::clusterExport(cl, varlist = c("pkg_path"), envir = environment())
  parallel::clusterEvalQ(cl, {
    devtools::load_all(pkg_path, quiet = TRUE)
    NULL
  })

  parallel::clusterExport(
    cl,
    varlist = c("params_cs_true", "params_cs_init", "run_one", "grid"),
    envir = environment()
  )

  res_list <- parallel::parLapply(cl, seq_len(nrow(grid)), function(i) {
    run_one(grid$T_end[i], grid$rep[i])
  })
}

elapsed_sec <- (proc.time() - pt0)[["elapsed"]]

results <- do.call(rbind, res_list)

# aggregate
agg <- do.call(rbind, lapply(split(results, results$T_end), function(d) {
  data.frame(
    T_end = d$T_end[1],
    reps = nrow(d),
    fit_median = median(d$fit_seconds),
    fit_mean   = mean(d$fit_seconds),
    fit_sd     = sd(d$fit_seconds),
    fit_p10    = unname(quantile(d$fit_seconds, 0.10)),
    fit_p90    = unname(quantile(d$fit_seconds, 0.90)),
    nodes_mean = mean(d$n_nodes),
    edges_mean = mean(d$n_edges),
    stringsAsFactors = FALSE
  )
}))

raw_file <- file.path(out_dir, sprintf("bench_fit_times_raw_%s.csv", version_label))
agg_file <- file.path(out_dir, sprintf("bench_fit_times_agg_%s.csv", version_label))

write.csv(results, raw_file, row.names = FALSE)
write.csv(agg, agg_file, row.names = FALSE)

cat("Elapsed_sec:", elapsed_sec, "\n")
cat("Saved:\n", normalizePath(raw_file), "\n", normalizePath(agg_file), "\n")
