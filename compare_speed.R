library(parallel)
library(ernm)
devtools::load_all(".")  # master loads it

params_cs_true <- list(
  mu = 5, K = 0.5, beta = 2, beta_edges = 0.5, node_lambda = 4,
  CS_star.2 = -1, CS_star.3 = -3
)

params_cs_init <- list(
  mu = 3, K = 1, beta = 1, beta_edges = 1, node_lambda = 5,
  CS_star.2 = 0, CS_star.3 = 0
)

times <- 1:2
reps  <- 10

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

  # be defensive: if structure differs, don’t crash the whole run
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

# choose cores
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "2"))
n_cores <- max(1L, n_cores)

# IMPORTANT: workers need to start in the same repo path
repo_dir <- normalizePath(getwd())

cl <- makeCluster(n_cores)

clusterSetRNGStream(cl, 1)

# Load packages on workers, set wd, load your package
clusterEvalQ(cl, {
  library(ernm)
  library(devtools)
  NULL
})

clusterExport(cl, "repo_dir", envir = environment())
clusterEvalQ(cl, {
  setwd(repo_dir)
  devtools::load_all(".")
  NULL
})

# export objects + function
clusterExport(
  cl,
  varlist = c("params_cs_true", "params_cs_init", "run_one", "grid"),
  envir = environment()
)

# run tasks (note: we pass T_end/rep directly, not via indexing issues)
res_list <- parLapply(cl, seq_len(nrow(grid)), function(i) {
  run_one(grid$T_end[i], grid$rep[i])
})

stopCluster(cl)

results <- do.call(rbind, res_list)
print(results)

# ---- save results ----

version_label <- Sys.getenv("VERSION_LABEL", "v1")

out_dir <- "bench_results"
dir.create(out_dir, showWarnings = FALSE)

raw_file <- file.path(out_dir, sprintf("bench_fit_times_raw_%s.csv", version_label))
agg_file <- file.path(out_dir, sprintf("bench_fit_times_agg_%s.csv", version_label))

# aggregate per T_end (median etc.)
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

write.csv(results, raw_file, row.names = FALSE)
write.csv(agg, agg_file, row.names = FALSE)

cat("Saved:\n",
    normalizePath(raw_file), "\n",
    normalizePath(agg_file), "\n")
