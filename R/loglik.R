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
#' @param formula_rhs String, the CS kernel formula
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
# Internal: build a compiled log-likelihood function that reuses mark caches
# (critical for optimizers that evaluate the objective many times).
.build_loglik_fn <- function(
    events, T_end = NULL, T0 = 0,
    mark_type = "ba",
    debug = FALSE,
    net_init_fun = net_init, net0 = NULL,
    net_step_fun = net_add_event,
    formula_rhs = NULL,
    truncation = Inf,
    ...
) {
  stopifnot(inherits(events, "events"))
  event_times <- events$times$t

  # In the case of no events
  if (length(event_times) == 0L) {
    return(function(params) 0.0)
  }

  # If the user does not specify T_end, assume the event window lasts exactly
  # up until the last event
  if (is.null(T_end)) T_end <- max(event_times)
  if (T_end < max(event_times)) stop("T_end must be >= max event time", call. = FALSE)

  # Grab the lists of node and edge arrivals by event
  arrivals <- nodes_by_event(events)
  edges <- edges_by_event(events)

  # Compile the mark kernel function.
  # This allow us to calculate the mark-likelihood contribution as a simple
  # function of the parameters.
  # The reason this works is that the network changes are fixed regardless of
  # the parameter values.
  mark_logpmf <- switch(
    mark_type,
    "ba" = compile_logpmf_ba(
      events = events,
      net_init_fun = net_init_fun,
      net0 = net0,
      net_step_fun = net_step_fun,
      delta = 0.001,
      ...
    ),
    "ba_bip" = compile_logpmf_ba_bip(
      events = events,
      net_init_fun = net_init_fun,
      net0 = net0,
      net_step_fun = net_step_fun,
      delta = 0.001,
      ...
    ),
    "cs" = compile_logpmf_cs(
      events = events,
      formula_rhs = formula_rhs,
      truncation = truncation,
      net_init_fun = net_init_fun,
      net0 = net0,
      net_step_fun = net_step_fun,
      model_cache = NULL,
      ...
    ),
    "cs_bip" = compile_logpmf_cs_bip(
      events = events,
      formula_rhs = formula_rhs,
      truncation = truncation,
      net_init_fun = net_init_fun,
      net0 = net0,
      net_step_fun = net_step_fun,
      model_cache = NULL,
      ...
    ),
    stop("mark_type must be one of 'ba', 'cs', 'ba_bip', 'cs_bip'", call. = FALSE)
  )



  # Debug: precompute layout from full network derived from events
  debug_layout <- NULL
  if (isTRUE(debug)) {
    full_net <- create_net_from_events(events)
    debug_layout <- .debug_prepare_full_layout(full_net)
  }

  # Return the likelihood function over params (natural scale)
  function(params) {
    mu <- params$mu
    K <- params$K
    beta <- params$beta

    if (!is.finite(mu) || !is.finite(K) || !is.finite(beta) || mu <= 0 || K < 0 || beta <= 0)
      return(-Inf)

    # mark loglik ONCE (it already sums over events)
    mark_ll_total <- mark_logpmf(params)
    if (!is.finite(mark_ll_total)) return(-Inf)

    ll_time <- 0.0
    net <- if (is.null(net0)) net_init_fun() else net0
    time_cache <- time_cache_init(T0 = T0)

    t_prev <- T0
    S_post_prev <- 0.0
    compensator <- 0.0

    for (idx in seq_along(event_times)) {
      t_k <- event_times[idx]
      new_nodes <- arrivals[[idx]]
      new_edges <- edges[[idx]]

      compensator <- compensator + compensator_inc_exp(
        dt = t_k - t_prev, mu = mu, K = K, beta = beta, S_post_prev = S_post_prev
      )

      time_cache <- time_cache_step(time_cache, t_k, beta = beta)
      lambda_g <- temporal_intensity_from_cache(time_cache, mu = mu, K = K, beta = beta)
      if (!is.finite(lambda_g) || lambda_g <= 0) return(-Inf)

      ll_time <- ll_time + log(lambda_g)

      # Debug plot (two panels)
      if (isTRUE(debug)) {
        net <- net_step_fun(net, new_nodes, new_edges, t_k)

        .debug_plot_step(net, debug_layout, t_k, 0, ll_time, ll_time + mark_ll_total, new_nodes,
                         new_edges, event_id = events$times$event_id[idx], edge_probs = NULL)
        readline("Press Enter / Return to continue to next event...")
      }

      time_cache$S <- time_cache$S + 1
      t_prev <- t_k
      S_post_prev <- time_cache$S
    }

    compensator <- compensator + compensator_inc_exp(
      dt = T_end - t_prev, mu = mu, K = K, beta = beta, S_post_prev = S_post_prev
    )

    (ll_time - compensator) + mark_ll_total
  }
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
#' @param formula_rhs String, the CS kernel formula
#'
#' @return A numeric scalar log-likelihood. Returns `-Inf` for invalid parameter
#'   values or if the likelihood is undefined.
#'
#' @export
loglik <- function(
    events, params, T_end = NULL, T0 = 0,
    mark_type = "ba",
    debug = FALSE,
    net_init_fun = net_init, net0 = NULL,
    net_step_fun = net_add_event,
    formula_rhs = NULL,
    truncation = Inf,
    ...
) {
  ll_fn <- .build_loglik_fn(
    events = events, T_end = T_end, T0 = T0,
    mark_type = mark_type, debug = debug,
    net_init_fun = net_init_fun, net0 = net0,
    net_step_fun = net_step_fun,
    formula_rhs = formula_rhs,
    truncation = truncation,
    ...
  )
  ll_fn(params)
}
