# fit_hawkes_growth.R

#' Fit Hawkes growth model by MLE (optim wrapper)
#'
#' Maximizes loglik(ev, params, ...) using optim with optional parameter transforms
#' and support for fixing some parameters.
#'
#' @param ev hg_events object (your canonical event stream).
#' @param params_init named list of initial values. Must include mu, K, beta.
#'        Can also include beta_edges (if your mark pmf uses it), lambda_new (BA-bip),
#'        node_lambda (CS), and others.
#' @param fixed_params character vector of parameter names to hold fixed.
#' @param method optim method; default "BFGS" works well with smooth transforms.
#' @param transform named character vector giving transform per parameter.
#'        Supported: "log" (positive) or "none".
#' @param hessian logical, ask optim for Hessian.
#' @param control list passed to optim.
#' @param ... passed through to loglik(), e.g. T0, T_end, mark_logpmf, mark_type, delta, etc.
#'
#' @return list with components fit (optim output), par (named list), loglik, hessian (if requested)
#' @export
fit_hg <- function(
    ev,
    params_init,
    fixed_params = NULL,
    method = "BFGS",
    transform = NULL,
    hessian = FALSE,
    control = list(),
    ...
) {
  stopifnot(is.list(params_init), length(params_init) > 0)

  start_time <- proc.time()

  # --- defaults: treat these as positive ---
  if (is.null(transform)) {
    transform <- rep("none", length(params_init))
    names(transform) <- names(params_init)
  }
  # ensure all init params have transform entry
  missing_tf <- setdiff(names(params_init), names(transform))
  if (length(missing_tf) > 0) {
    transform[missing_tf] <- "none"
  }

  # Default to log transform for common positive parameters
  for (nm in intersect(
    c("mu", "K", "beta", "beta_edges", "lambda_new", "node_lambda"),
    names(transform)
  )) {
    if (transform[[nm]] == "none") transform[[nm]] <- "log"
  }

  # --- expand any vector-valued parameters (e.g. CS_params) into scalar entries ---
  vec_map <- list()
  nms0 <- names(params_init)

  for (nm in nms0) {
    val <- params_init[[nm]]
    if (is.atomic(val) && length(val) > 1L) {
      new_nms <- paste0(nm, "_", seq_along(val))
      vec_map[[nm]] <- new_nms

      # expand params_init into scalar entries
      params_init[new_nms] <- as.list(as.numeric(val))
      params_init[[nm]] <- NULL

      # expand transform: replicate transform[[nm]] if it exists, else "none"
      tf_nm <- transform[[nm]]
      if (is.null(tf_nm)) tf_nm <- "none"
      transform[new_nms] <- rep(tf_nm, length(new_nms))

      # delete old transform entry safely (vector)
      transform <- transform[setdiff(names(transform), nm)]

      # if user fixed the vector param by base name, expand to the new names later
      if (!is.null(fixed_params) && nm %in% fixed_params) {
        fixed_params <- setdiff(fixed_params, nm)
        fixed_params <- c(fixed_params, new_nms)
      }
    }
  }

  # ensure transform entries exist for ALL (expanded) params
  missing_tf <- setdiff(names(params_init), names(transform))
  if (length(missing_tf) > 0) transform[missing_tf] <- "none"

  # --- split free vs fixed ---
  fixed_params <- fixed_params %||% character(0)
  fixed_params <- intersect(fixed_params, names(params_init))
  free_params  <- setdiff(names(params_init), fixed_params)

  if (length(free_params) == 0L) {
    ll <- loglik(ev = ev, params = params_init, ...)
    return(list(
      fit = NULL,
      par = params_init,
      loglik = ll,
      convergence = 0L,
      message = "No free parameters (all fixed)."
    ))
  }

  # --- helpers: forward/backward transforms ---
  to_working <- function(p_list) {
    x <- numeric(length(free_params))
    names(x) <- free_params
    for (nm in free_params) {
      val <- p_list[[nm]]
      tf  <- transform[[nm]]
      if (tf == "log") {
        if (!is.finite(val) || val <= 0) stop("Initial value for ", nm, " must be > 0 for log-transform.")
        x[[nm]] <- log(val)
      } else if (tf == "none") {
        x[[nm]] <- val
      } else {
        stop("Unknown transform '", tf, "' for param ", nm)
      }
    }
    x
  }

  from_working <- function(x) {
    p <- params_init
    for (nm in free_params) {
      tf <- transform[[nm]]
      if (tf == "log") {
        p[[nm]] <- exp(x[[nm]])
      } else {
        p[[nm]] <- x[[nm]]
      }
    }
    for (nm in fixed_params) p[[nm]] <- params_init[[nm]]

    # collapse any expanded vector params back into vectors (e.g. CS_params)
    if (length(vec_map)) {
      for (base_nm in names(vec_map)) {
        new_nms <- vec_map[[base_nm]]
        p[[base_nm]] <- unname(as.numeric(p[new_nms]))
        p[new_nms] <- NULL
      }
    }

    p
  }

  # --- objective for optim: MINIMIZE negative loglik ---
  fn <- function(x) {
    p <- from_working(x)
    ll <- tryCatch(
      loglik(ev = ev, params = p, ...),
      error = function(e) -Inf
    )
    if (!is.finite(ll)) return(1e100)
    -ll
  }

  x0 <- to_working(params_init)

  control2 <- modifyList(list(fnscale = 1), control)

  fit <- stats::optim(
    par = x0,
    fn = fn,
    method = method,
    control = control2,
    hessian = hessian
  )

  p_hat <- from_working(fit$par)
  ll_hat <- loglik(ev = ev, params = p_hat, ...)

  out <- list(
    fit = fit,
    par = p_hat,
    loglik = ll_hat,
    convergence = fit$convergence,
    message = fit$message
  )
  if (hessian) out$hessian <- fit$hessian

  print(paste0("Fitting took ", round((proc.time()-start_time)[3],2), " seconds"))
  out
}

# tiny helper (so you don't need rlang)
`%||%` <- function(x, y) if (is.null(x)) y else x
