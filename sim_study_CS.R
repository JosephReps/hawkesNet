## sim_study_CS.R
## Cluster-friendly simulation study for fixNet CS Hawkes growth:
## simulate -> fit -> save results

suppressPackageStartupMessages({
  library(parallel)
  library(ernm)
})

# ---- Helpers ----
log_msg <- function(..., log_con = NULL) {
  msg <- paste0(...)
  stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", stamp, "] ", msg, "\n")
  cat(line)
  if (!is.null(log_con)) writeLines(line, con = log_con, sep = "")
}

as_int <- function(x, default) {
  if (!nzchar(x)) return(as.integer(default))
  suppressWarnings({
    y <- as.integer(x)
    if (is.na(y)) as.integer(default) else y
  })
}

as_num <- function(x, default) {
  if (!nzchar(x)) return(as.numeric(default))
  suppressWarnings({
    y <- as.numeric(x)
    if (is.na(y)) as.numeric(default) else y
  })
}


as_num_vec <- function(x, default = numeric(0)) {
  if (!nzchar(x)) return(default)
  # allow commas and/or whitespace
  parts <- unlist(strsplit(x, "[,[:space:]]+"))
  parts <- parts[nzchar(parts)]
  suppressWarnings({
    v <- as.numeric(parts)
  })
  if (anyNA(v)) stop("Could not parse numeric vector from: ", x)
  v
}

# ---- Settings (override via env vars) ----
TIME   <- as_num(Sys.getenv("TIME", ""), 20)
T0     <- as_num(Sys.getenv("T0", ""), 0)
N_SIMS <- as_int(Sys.getenv("N_SIMS", ""), 4)
SEED   <- as_int(Sys.getenv("SEED", ""), 1267)

# ---- CS kernel controls (override via env vars) ----
FORMULA_RHS   <- Sys.getenv("FORMULA_RHS", "edges  + triangles() + star(c(2,3))")
TRUNCATION    <- as_int(Sys.getenv("TRUNCATION", ""), 500)
MAX_NODE_TIME <- as_num(Sys.getenv("MAX_NODE_TIME", ""), 10000)

# CS params can be passed as a vector string, e.g.
#   TRUE_CS_PARAMS="-6,0.5,0.3,-0.1"   (commas or spaces both OK)
# Back-compat: if TRUE_CS_PARAMS/INIT_CS_PARAMS not set, fall back to TRUE_CS_1/INIT_CS_1.
CS_true <- if (nzchar(Sys.getenv("TRUE_CS_PARAMS", ""))) {
  as_num_vec(Sys.getenv("TRUE_CS_PARAMS", ""), default = numeric(0))
} else {
  c(as_num(Sys.getenv("TRUE_CS_1", ""), 0.25))
}
CS_true <- c(-6,0.5,0.3,-0.1)

CS_init <- if (nzchar(Sys.getenv("INIT_CS_PARAMS", ""))) {
  as_num_vec(Sys.getenv("INIT_CS_PARAMS", ""), default = numeric(0))
} else {
  c(as_num(Sys.getenv("INIT_CS_1", ""), 0.10))
}
CS_init <- c(-10,0,0,0)

if (length(CS_true) != length(CS_init)) {
  stop(
    "CS_params length mismatch: TRUE has length ", length(CS_true),
    " but INIT has length ", length(CS_init)
  )
}
if (length(CS_true) != length(CS_init)) {
  stop("CS_params length mismatch: TRUE has length ", length(CS_true),
       " but INIT has length ", length(CS_init))
}

# True parameters (CS adds node_lambda + CS_params)
params_true <- list(
  mu = as_num(Sys.getenv("TRUE_MU", ""), 1),
  K = as_num(Sys.getenv("TRUE_K", ""), 0.5),
  beta = as_num(Sys.getenv("TRUE_BETA", ""), 2),
  beta_edges = as_num(Sys.getenv("TRUE_BETA_EDGES", ""), 0.5),
  node_lambda = as_num(Sys.getenv("TRUE_NODE_LAMBDA", ""), 2.0),
  CS_params = CS_true
)

# Initial values for optimizer
params_init <- list(
  mu = as_num(Sys.getenv("INIT_MU", ""), 1),
  K = as_num(Sys.getenv("INIT_K", ""), 0.1),
  beta = as_num(Sys.getenv("INIT_BETA", ""), 1),
  beta_edges = as_num(Sys.getenv("INIT_BETA_EDGES", ""), 1),
  node_lambda = as_num(Sys.getenv("INIT_NODE_LAMBDA", ""), 1.0),
  CS_params = CS_init
)

# mu,   K,     beta,  beta_edges, node_lambda,  CS_params
# log,  log,   log,   log,        log,          none...
TRANSFORM <- c(rep("log", 5), "none")
names(TRANSFORM) <- names(params_true)

# Output directory + run id
out_dir <- Sys.getenv("OUT_DIR", "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
run_id <- Sys.getenv("RUN_ID", "")
if (!nzchar(run_id)) run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")

log_path <- file.path(out_dir, paste0("sim_study_CS_", run_id, ".log"))
log_con <- file(log_path, open = "wt")
on.exit({ try(close(log_con), silent = TRUE) }, add = TRUE)

# Parallel cores (AWS-friendly: set N_CORES, fallback to SLURM_CPUS_PER_TASK)
N_CORES <- if (nzchar(Sys.getenv("N_CORES", ""))) {
  as_int(Sys.getenv("N_CORES", ""), 1)
} else {
  as_int(Sys.getenv("SLURM_CPUS_PER_TASK", ""), 4)
}
N_CORES <- max(1L, N_CORES)

# Support Slurm array (one replicate per job) if present
array_task <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
rep_ids <- if (nzchar(array_task)) {
  as.integer(array_task)
} else {
  seq_len(N_SIMS)
}

log_msg(
  "Starting sim_study_CS | run_id=", run_id,
  " | reps=", length(rep_ids),
  " | cores=", N_CORES,
  " | TIME=", TIME,
  " | T0=", T0,
  " | SEED=", SEED,
  " | FORMULA_RHS=", FORMULA_RHS,
  " | TRUNCATION=", TRUNCATION,
  " | MAX_NODE_TIME=", MAX_NODE_TIME,
  log_con = log_con
)

# ---- Worker function: one replicate ----
one_rep <- function(rep_id) {
  set.seed(SEED + rep_id)

  # 1) simulate (CS)
  sim <- sim_hg(
    params = params_true,
    T_end = TIME,
    T0 = T0,
    return_ev = TRUE,
    return_state = FALSE,
    debug = FALSE,
    sim_mark_fun = sim_mark_cs,
    mark_type = "cs",
    formula_rhs = FORMULA_RHS,
    truncation = TRUNCATION,
    max_node_time = MAX_NODE_TIME
  )

  nodes_n  <- if (!is.null(sim$ev$nodes)) nrow(sim$ev$nodes) else NA_integer_
  edges_n  <- if (!is.null(sim$ev$edges)) nrow(sim$ev$edges) else NA_integer_
  events_n <- if (!is.null(sim$events$t)) length(sim$events$t) else NA_integer_

  # 2) fit (CS)
  fit <- tryCatch(
    fit_hg(
      ev = sim$ev,
      params_init = params_init,
      method = "BFGS",
      hessian = FALSE,
      control = list(maxit = as_int(Sys.getenv("MAXIT", ""), 500)),
      T0 = T0,
      T_end = TIME,
      transform = TRANSFORM,
      mark_logpmf = log_pmf_cs,
      mark_type = "cs",
      formula_rhs = FORMULA_RHS,
      truncation = TRUNCATION,
      max_node_time = MAX_NODE_TIME
    ),
    error = function(e) list(error = e$message)
  )

  if (!is.null(fit$error)) {
    return(list(
      ok = FALSE,
      rep = rep_id,
      error = fit$error,
      nodes_n = nodes_n,
      edges_n = edges_n,
      events_n = events_n
    ))
  }

  list(
    ok = TRUE,
    rep = rep_id,
    par = fit$par,
    loglik = fit$loglik,
    convergence = fit$convergence,
    message = fit$message,
    nodes_n = nodes_n,
    edges_n = edges_n,
    events_n = events_n
  )
}

safe_one_rep <- function(rep_id) {
  out <- tryCatch(
    one_rep(rep_id),
    error = function(e) {
      list(
        ok = FALSE,
        rep = rep_id,
        error = conditionMessage(e),
        nodes_n = NA_integer_,
        edges_n = NA_integer_,
        events_n = NA_integer_
      )
    }
  )

  # Final guard: ensure consistent structure even if something weird comes back
  if (!is.list(out) || is.null(out$ok)) {
    out <- list(
      ok = FALSE,
      rep = rep_id,
      error = "Non-list result returned from one_rep()",
      nodes_n = NA_integer_,
      edges_n = NA_integer_,
      events_n = NA_integer_
    )
  }

  out
}

# ---- Progress helpers (live bar, works with multi-core PSOCK clusters) ----
lapply_progress <- function(x, fun) {
  n <- length(x)
  if (n == 0L) return(list())
  pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
  on.exit(close(pb), add = TRUE)
  out <- vector('list', n)
  for (i in seq_len(n)) {
    out[[i]] <- fun(x[[i]])
    utils::setTxtProgressBar(pb, i)
  }
  out
}

parLapply_progress <- function(cl, x, fun) {
  n <- length(x)
  if (n == 0L) return(list())
  pb <- utils::txtProgressBar(min = 0, max = n, style = 3)
  on.exit(close(pb), add = TRUE)
  p <- length(cl)
  res <- vector('list', n)
  i <- 0L
  completed <- 0L
  # prime the workers
  for (w in seq_len(min(p, n))) {
    i <- i + 1L
    parallel:::sendCall(cl[[w]], fun, list(x[[i]]), tag = i)
  }
  while (completed < n) {
    d <- parallel:::recvOneResult(cl)
    completed <- completed + 1L
    res[[d$tag]] <- d$value
    utils::setTxtProgressBar(pb, completed)
    if (i < n) {
      i <- i + 1L
      parallel:::sendCall(cl[[d$node]], fun, list(x[[i]]), tag = i)
    }
  }
  res
}

# ---- Run sims (parallel) ----
pt0 <- proc.time()

if (N_CORES == 1L || length(rep_ids) == 1L) {
  # Ensure package is available in this R session
  # (Assumes you run from the package repo; otherwise install & library(yourpkg))
  if (!requireNamespace('devtools', quietly = TRUE)) {
    stop('devtools is required for load_all() in this script. Install it, or modify script to library(yourPackage).')
  }
  devtools::load_all(normalizePath(getwd()), quiet = TRUE)

  res <- lapply_progress(as.list(rep_ids), safe_one_rep)

} else {
  cl <- parallel::makeCluster(N_CORES)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  pkg_path <- normalizePath(getwd())

  parallel::clusterEvalQ(cl, {
    library(stats)
    library(devtools)
    library(ernm)
  })

  parallel::clusterExport(cl, varlist = c('pkg_path'), envir = environment())

  parallel::clusterEvalQ(cl, {
    devtools::load_all(pkg_path, quiet = TRUE)
  })

  parallel::clusterExport(
    cl,
    varlist = c(ls(envir = environment())),
    envir = environment()
  )

  # Run with a live progress bar
  res <- parLapply_progress(cl, as.list(rep_ids), safe_one_rep)
}

elapsed_sec <- (proc.time() - pt0)[['elapsed']]

# ---- Collect results ----
ok_res  <- res[vapply(res, function(x) isTRUE(x$ok), logical(1))]
bad_res <- res[vapply(res, function(x) isFALSE(x$ok), logical(1))]

log_msg("Finished | OK=", length(ok_res), " FAIL=", length(bad_res), " | elapsed_sec=", elapsed_sec, log_con = log_con)
if (length(bad_res)) {
  log_msg("Example failure: ", bad_res[[1]]$error, log_con = log_con)
}

# Convert fits to data.frame (base R)
est_df <- if (length(ok_res)) {
  do.call(rbind, lapply(ok_res, function(x) {

    cs <- unlist(x$par$CS_params, use.names = FALSE)
    cs_cols <- setNames(as.list(cs), paste0("CS_params_", seq_along(cs)))

    rest_of_cols <- data.frame(
      rep = x$rep,
      convergence = x$convergence,
      loglik = x$loglik,
      nodes_n = x$nodes_n,
      edges_n = x$edges_n,
      events_n = x$events_n,
      mu = x$par$mu,
      K = x$par$K,
      beta = x$par$beta,
      beta_edges = x$par$beta_edges,
      node_lambda = x$par$node_lambda,
      stringsAsFactors = FALSE
    )

    data.frame(rest_of_cols, cs_cols, check.names = FALSE, stringsAsFactors = FALSE)
  }))
} else {
  data.frame()
}

fail_df <- if (length(bad_res)) {
  do.call(rbind, lapply(bad_res, function(x) {
    data.frame(
      rep = x$rep,
      error = x$error,
      nodes_n = x$nodes_n,
      edges_n = x$edges_n,
      events_n = x$events_n,
      stringsAsFactors = FALSE
    )
  }))
} else {
  data.frame()
}

# ---- Save outputs ----
obj <- list(
  study = "CS",
  run_id = run_id,
  rep_ids = rep_ids,
  timing = list(elapsed_sec = elapsed_sec),
  config = list(
    TIME = TIME, T0 = T0, N_SIMS = N_SIMS, SEED = SEED,
    FORMULA_RHS = FORMULA_RHS,
    TRUNCATION = TRUNCATION,
    MAX_NODE_TIME = MAX_NODE_TIME,
    params_true = params_true,
    params_init = params_init
  ),
  est_df = est_df,
  fail_df = fail_df
)
out_path <- file.path(out_dir, paste0("sim_study_CS_", run_id, ".rds"))
saveRDS(obj, out_path)

log_msg("Saved results: ", out_path, log_con = log_con)
log_msg("Saved log: ", log_path, log_con = log_con)
