# loglik.R

# Internal: initialize recursion cache for an exponential Hawkes kernel.
time_cache_init <- function(T0) {
  list(last_t = T0, S = 0.0, n = 0L)
}

time_cache_step <- function(cache, t, beta) {
  dt <- t - cache$last_t
  if (dt < 0) stop("times must be increasing")
  cache$S <- exp(-beta * dt) * cache$S
  cache$last_t <- t
  cache
}

temporal_intensity_from_cache <- function(cache, mu, K, beta, kernel_density = FALSE) {
  scale <- if (kernel_density) beta else 1.0
  mu + K * scale * cache$S
}

compensator_inc_exp <- function(dt, mu, K, beta, S_post_prev) {
  if (dt < 0) stop("times must be increasing")
  mu * dt + (K / beta) * S_post_prev * (1 - exp(-beta * dt))
}

#' Log-likelihood for a HawkesNet model
#'
#' Computes the log-likelihood for the observed event sequence under an
#' exponential temporal Hawkes kernel and a chosen mark model (`mark_type`).
#'
#' @inheritParams hawkesnet-lik-args
#' @param params Named list of model parameters on the natural scale. Must
#'   contain `mu`, `K`, and `beta`, plus any mark parameters required by
#'   `mark_type` (e.g. `beta_edges` for BA).
#' @param net_init_fun Function to initialize the pre-event network if `net0` is
#'   `NULL`. Defaults to `net_init`.
#' @param net0 Optional initial network object to start recursion from.
#' @param net_step_fun Function used to update the network after each event.
#'   Defaults to `net_add_event`.
#'
#' @return A numeric scalar log-likelihood. Returns `-Inf` for invalid parameter
#'   values or if the likelihood is undefined.
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' params_true <- list(mu = 0.5, K = 0.5, beta = 0.5, beta_edges = 2)
#' sim <- sim_hawkesNet(params = params_true, T_end = 10, mark_type = "ba")
#' loglik(sim$ev, params_true, mark_type = "ba")
loglik <- function(
    events, params, T_end = NULL, T0 = 0,
    mark_type = "ba",
    debug = FALSE,
    net_init_fun = net_init, net0 = NULL,
    net_step_fun = net_add_event,
    ... # passed to mark kernels
) {
  # Pull out common parameters from params object
  mu <- params$mu
  K <- params$K
  beta <- params$beta

  if (!is.finite(mu) || !is.finite(K) || !is.finite(beta) || mu <= 0 || K < 0 || beta <= 0)
    return(-Inf)

  event_times <- events$times$t
  if (length(event_times) == 0L) return(0)

  if (is.null(T_end)) T_end <- max(event_times)
  if (T_end < max(event_times)) stop("T_end must be >= max event time")

  # Mark kernel dispatch
  if (mark_type == "ba") {
    mark_logpmf <- log_pmf_ba
  } else if (mark_type == "ba_bip") {
    mark_logpmf <- log_pmf_ba_bip
  } else if (mark_type == "cs") {
    mark_logpmf <- log_pmf_cs
  } else {
    stop("mark_type must be one of 'ba', 'cs', 'ba_bip'")
  }

  arrivals <- nodes_by_event(events)
  edges <- edges_by_event(events)

  ll <- 0.0
  net <- if (is.null(net0)) net_init_fun() else net0

  time_cache <- time_cache_init(T0 = T0)

  # For compensator
  t_prev <- T0
  S_post_prev <- 0.0
  compensator <- 0.0

  # CS: build model cache once
  cs_model_cache <- NULL
  if (mark_type == "cs") cs_model_cache <- list(rhs = NULL, model = NULL)

  # Debug: precompute layout from full network derived from events
  debug_layout <- NULL
  if (debug == TRUE) {
    full_net <- create_net_from_events(events)
    debug_layout <- .debug_prepare_full_layout(full_net)
  }

  for (idx in seq_along(event_times)) {
    t_k <- event_times[idx]
    new_nodes <- arrivals[[idx]]
    new_edges <- edges[[idx]]

    # (1) compensator over gap
    compensator <- compensator + compensator_inc_exp(
      dt = t_k - t_prev, mu = mu, K = K, beta = beta, S_post_prev = S_post_prev
    )

    # (2) decay cache to t_k
    time_cache <- time_cache_step(time_cache, t_k, beta = beta)

    # (3) temporal intensity
    lambda_g <- temporal_intensity_from_cache(time_cache, mu = mu, K = K, beta = beta)
    if (!is.finite(lambda_g) || lambda_g <= 0) return(-Inf)

    # (4) mark log-pmf using pre-event net
    edge_probs <- NULL
    if (mark_type == "cs") {
      res <- mark_logpmf(net, new_nodes, new_edges, t_k, params, model_cache = cs_model_cache, ...)
      mark_ll <- res$logp
      edge_probs <- res$edge_probs
      cs_model_cache <- res$model_cache
    } else {
      res <- mark_logpmf(net, new_nodes, new_edges, t_k, params, ...)
      mark_ll <- res$logp
      edge_probs <- res$edge_probs
    }
    if (!is.finite(mark_ll)) return(-Inf)

    # (5) update loglik
    ll <- ll + mark_ll + log(lambda_g)

    # (6) update net state
    net <- net_step_fun(net, new_nodes, new_edges, t_k)

    # Debug plot (two panels)
    if (isTRUE(debug)) {
      .debug_plot_step(net, debug_layout, t_k, mark_ll, lambda_g, ll, new_nodes,
                       new_edges, event_id = events$times$event_id[idx], edge_probs = edge_probs)
      readline("Press Enter / Return to continue to next event...")
    }

    # (7) S_post update
    time_cache$S <- time_cache$S + 1
    t_prev <- t_k
    S_post_prev <- time_cache$S
  }

  compensator <- compensator + compensator_inc_exp(
    dt = T_end - t_prev, mu = mu, K = K, beta = beta, S_post_prev = S_post_prev
  )

  ll - compensator
}
