# diagnostics_temporal.R

#' Temporal diagnostics for HawkesNet
#'
#' Computes basic time-series diagnostics for the exponential Hawkes temporal
#' component at observed event times, including inter-arrival times, intensity
#' at events, compensator increments, and time-rescaled residuals.
#'
#' @param events An events object created by `make_events()`. Must contain
#'   `events$times$t` (numeric event times).
#' @param params Named list of model parameters on the natural scale. Must
#'   contain `mu`, `K`, `beta`.
#' @param T_end Optional end time for diagnostics window. Defaults to last event.
#' @param T0 Start time for diagnostics window.
#'
#' @return A list with components:
#' \describe{
#'   \item{t}{Event times.}
#'   \item{dt}{Inter-arrival times (first is `t[1] - T0`).}
#'   \item{lambda}{Temporal intensity evaluated just before each event time.}
#'   \item{comp_inc}{Compensator increments over each gap (same length as `t`).}
#'   \item{tau}{Time-rescaled increments (should be Exp(1) under correct temporal model).}
#'   \item{meta}{List with `T0`, `T_end`, and basic checks.}
#' }
#' @export
diag_temporal <- function(events, params, T_end = NULL, T0 = 0) {
  if (is.null(events) || is.null(events$times) || !"t" %in% names(events$times)) {
    stop("diag_temporal(): `events$times$t` is required.", call. = FALSE)
  }

  t <- as.numeric(events$times$t)
  if (length(t) == 0L) {
    return(list(
      t = numeric(0),
      dt = numeric(0),
      lambda = numeric(0),
      comp_inc = numeric(0),
      tau = numeric(0),
      meta = list(T0 = T0, T_end = T_end, ok = TRUE, note = "No events.")
    ))
  }
  if (any(!is.finite(t))) stop("diag_temporal(): non-finite event times.", call. = FALSE)
  if (is.unsorted(t, strictly = TRUE) && any(diff(t) < 0)) {
    stop("diag_temporal(): event times must be nondecreasing.", call. = FALSE)
  }

  mu <- params$mu
  K <- params$K
  beta <- params$beta
  if (!is.finite(mu) || !is.finite(K) || !is.finite(beta) || mu <= 0 || K < 0 || beta <= 0) {
    return(list(
      t = t,
      dt = c(t[1] - T0, diff(t)),
      lambda = rep(NA_real_, length(t)),
      comp_inc = rep(NA_real_, length(t)),
      tau = rep(NA_real_, length(t)),
      meta = list(T0 = T0, T_end = T_end, ok = FALSE, note = "Invalid mu/K/beta.")
    ))
  }

  if (is.null(T_end)) T_end <- max(t)
  if (!is.finite(T_end) || T_end < max(t)) stop("diag_temporal(): T_end must be >= max(t).", call. = FALSE)

  # Helpers (duplicated here to keep this file standalone)
  cache_init <- function(T0) list(last_t = T0, S = 0.0)
  cache_decay_to <- function(cache, t_new, beta) {
    dt <- t_new - cache$last_t
    if (dt < 0) stop("diag_temporal(): times must be increasing.", call. = FALSE)
    cache$S <- exp(-beta * dt) * cache$S
    cache$last_t <- t_new
    cache
  }
  lambda_from_cache <- function(cache, mu, K) mu + K * cache$S
  compensator_inc <- function(dt, mu, K, beta, S_post_prev) {
    if (dt < 0) stop("diag_temporal(): times must be increasing.", call. = FALSE)
    mu * dt + (K / beta) * S_post_prev * (1 - exp(-beta * dt))
  }

  n <- length(t)
  dt <- c(t[1] - T0, diff(t))
  lambda <- numeric(n)
  comp_inc <- numeric(n)

  cache <- cache_init(T0)
  t_prev <- T0
  S_post_prev <- 0.0

  for (k in seq_len(n)) {
    t_k <- t[k]

    # Compensator over gap (t_prev, t_k)
    comp_inc[k] <- compensator_inc(
      dt = t_k - t_prev, mu = mu, K = K, beta = beta, S_post_prev = S_post_prev
    )

    # Decay cache to just before t_k
    cache <- cache_decay_to(cache, t_k, beta = beta)
    lambda[k] <- lambda_from_cache(cache, mu = mu, K = K)

    # Update post-event: jump by 1 at event time (standard Hawkes)
    cache$S <- cache$S + 1
    t_prev <- t_k
    S_post_prev <- cache$S
  }

  # Time-rescaled increments (tau_k = compensator increment over each gap)
  tau <- comp_inc

  list(
    t = t,
    dt = dt,
    lambda = lambda,
    comp_inc = comp_inc,
    tau = tau,
    meta = list(T0 = T0, T_end = T_end, ok = TRUE)
  )
}


#' Plot temporal diagnostics from `diag_temporal()`
#'
#' @param diag Output of `diag_temporal()`.
#' @param which Character vector of plots to show. Any of:
#'   `"count"`, `"dt"`, `"lambda"`, `"tau_hist"`, `"tau_qq"`.
#' @param ask Logical; if `TRUE`, prompt between plots.
#' @export
plot_temporal_diag <- function(diag,
                               which = c("count", "dt", "lambda", "tau_hist", "tau_qq"),
                               ask = FALSE) {
  if (is.null(diag) || is.null(diag$t)) stop("plot_temporal_diag(): invalid `diag`.", call. = FALSE)

  t <- diag$t
  if (length(t) == 0L) {
    plot.new()
    text(0.5, 0.5, "No events to plot.")
    return(invisible(NULL))
  }

  old_ask <- par(ask = ask)
  on.exit(par(old_ask), add = TRUE)

  which <- unique(which)

  if ("count" %in% which) {
    plot(t, seq_along(t), type = "s", xlab = "time", ylab = "N(t)",
         main = "Cumulative event count")
  }

  if ("dt" %in% which) {
    dt <- diag$dt
    hist(dt, breaks = "FD", xlab = "Δt", main = "Inter-arrival times", col = "grey90", border = "grey40")
  }

  if ("lambda" %in% which) {
    lam <- diag$lambda
    plot(t, lam, type = "l", xlab = "time", ylab = expression(lambda(t[k]^-1)),
         main = "Intensity at event times", lty = 1)
    # lines(t, lam, lty = 1)
  }

  if ("tau_hist" %in% which) {
    tau <- diag$tau
    hist(tau, breaks = "FD", probability = TRUE, xlab = expression(tau[k]),
         main = "Time-rescaled increments", col = "grey90", border = "grey40")
    xx <- seq(0, max(tau, na.rm = TRUE), length.out = 200)
    lines(xx, dexp(xx), lty = 1)
    legend("topright", legend = c("hist", "Exp(1) density"), bty = "n", lty = c(NA, 1), pch = c(15, NA))
  }

  if ("tau_qq" %in% which) {
    tau <- diag$tau
    tau <- tau[is.finite(tau) & tau >= 0]
    if (length(tau) < 5) {
      plot.new()
      text(0.5, 0.5, "Not enough tau values for QQ plot.")
    } else {
      n <- length(tau)
      theo <- stats::qexp(ppoints(n))
      samp <- sort(tau)
      plot(theo, samp, xlab = "Exp(1) theoretical quantiles", ylab = "Sample quantiles",
           main = "QQ plot: tau vs Exp(1)")
      abline(0, 1, lty = 2)
    }
  }

  invisible(NULL)
}
