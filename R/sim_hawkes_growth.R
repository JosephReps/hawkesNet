# sim_hawkes_growth.R

#' Title
#'
#' @param params
#' @param T_end
#' @param T0
#' @param mu_multiplier
#' @param state_init_fun
#' @param state_step_fun
#' @param sim_mark_fun
#' @param return_ev
#' @param return_state
#' @param debug
#'
#' @return
#' @export
#'
#' @examples
sim_hg_batched <- function(
    params,
    T_end,
    T0 = 0,
    mu_multiplier = 10,
    state_init_fun = state_init,
    state_step_fun = state_step,
    sim_mark_fun = sim_mark_ba,
    delta = 0.5,
    # outputs
    return_ev = TRUE,
    return_state = FALSE,
    debug = FALSE,
    ...
) {
  start_time <- proc.time()
  mu <- params$mu
  K <- params$K
  beta <- params$beta

  if (!is.finite(mu) || !is.finite(K) || !is.finite(beta) || mu <= 0 || K < 0 || beta <= 0)
    stop("Invalid params (mu>0, K>=0, beta>0 required)")

  lambda_upper <- mu_multiplier * mu
  if (!is.finite(lambda_upper) || lambda_upper <= 0) stop("mu_multiplier*mu must be > 0")

  # Determine mark type (strict switch)
  mark_type <- NULL
  if (exists("sim_mark_ba_bip", mode = "function") && identical(sim_mark_fun, sim_mark_ba_bip)) {
    mark_type <- "ba_bip"
  } else if (exists("sim_mark_cs", mode = "function") && identical(sim_mark_fun, sim_mark_cs)) {
    mark_type <- "cs"
  } else {
    mark_type <- "ba"
  }
  mark_type <- match.arg(mark_type, c("ba", "ba_bip", "cs"))

  # CS: model cache for ernm
  cs_model_cache <- NULL
  if (mark_type == "cs") {
    cs_model_cache <- new.env(parent = emptyenv())
    cs_model_cache$rhs <- NULL
    cs_model_cache$model <- NULL
  }

  # 1) Propose candidate times from homogeneous Poisson upper bound
  n_bg <- stats::rpois(1, lambda_upper * (T_end - T0))
  cand_times <- sort(stats::runif(n_bg, min = T0, max = T_end))

  # 2) Initialize state
  state <- state_init_fun()
  cache <- time_cache_init(T0 = T0)

  events_t <- numeric(0)
  accept_probs <- numeric(0)

  # For building hg_events at the end
  nodes_rows <- list()
  edges_rows <- list()

  # 3) Thinning loop
  for (t_k in cand_times) {
    # Decay cache to candidate time
    cache <- time_cache_step(cache, t_k, beta = beta)
    lambda_g <- temporal_intensity_from_cache(cache, mu = mu, K = K, beta = beta)

    if (!is.finite(lambda_g) || lambda_g <= 0) next

    a <- lambda_g / lambda_upper
    accept_probs <- c(accept_probs, a)

    if (debug) message("t_k=", t_k, "  S_pre=", cache$S, "  lambda_g=", lambda_g, "  a=", a)

    if (stats::runif(1) < a) {
      # 4) Accept: sample mark given current state
      if (mark_type == "cs") {
        mk <- sim_mark_fun(state = state, t_k = t_k, params = params, model_cache = cs_model_cache, ...)
      } else {
        mk <- sim_mark_fun(state = state, t_k = t_k, params = params, delta = delta, ...)
      }

      arrivals_k <- mk$arrivals
      edges_k <- mk$edges

      # BA-bip roles: event first then perps
      role_k <- NULL
      if (mark_type == "ba_bip") {
        if (length(arrivals_k) < 1L) stop("BA-bip strict: sim_mark_fun produced no arrivals.")
        role_k <- c("event", rep("perp", length(arrivals_k) - 1L))
      }

      # 5) Record to tables for make_events()
      if (length(arrivals_k) > 0) {
        df_nodes <- data.frame(
          id = arrivals_k,
          time = rep(t_k, length(arrivals_k)),
          stringsAsFactors = FALSE
        )
        if (!is.null(role_k)) df_nodes$role <- role_k
        nodes_rows[[length(nodes_rows) + 1L]] <- df_nodes
      }
      if (!is.null(edges_k) && nrow(edges_k) > 0) {
        edges_rows[[length(edges_rows) + 1L]] <- data.frame(
          i = as.character(edges_k$i),
          j = as.character(edges_k$j),
          time = rep(t_k, nrow(edges_k)),
          stringsAsFactors = FALSE
        )
      }

      # 6) Update state
      if (is.null(role_k)) {
        state <- state_step_fun(state, arrivals_k, edges_k, t_k)
      } else {
        state <- state_step_fun(state, arrivals_k, edges_k, t_k, role_k = role_k)
      }

      # 7) Update temporal cache: add event contribution (S_post = S_pre + 1)
      cache$S <- cache$S + 1

      # 8) Save accepted time
      events_t <- c(events_t, t_k)
    }
  }

  out <- list(
    events = list(n = length(events_t), t_k = events_t),
    accept_probs = accept_probs
  )

  # Build ev if requested
  if (return_ev) {
    nodes <- if (length(nodes_rows)) do.call(rbind, nodes_rows) else data.frame(id=character(0), time=numeric(0))
    edges <- if (length(edges_rows)) do.call(rbind, edges_rows) else data.frame(i=character(0), j=character(0), time=numeric(0))

    out$ev <- make_events(nodes = nodes, edges = edges)
  }

  if (return_state) out$state <- state

  print(paste0("Simulation took ", round((proc.time()-start_time)[3],2), " seconds"))

  out
}


#' Title
#'
#' @param params
#' @param T_end
#' @param T0
#' @param state_init_fun
#' @param state_step_fun
#' @param sim_mark_fun
#' @param delta
#' @param mark_type
#' @param return_ev
#' @param return_state
#' @param debug
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
sim_hg <- function(
    params, T_end, T0 = 0,
    state_init_fun = state_init,
    state_step_fun = state_step,
    sim_mark_fun = sim_mark_ba,
    delta = 0.001,
    mark_type = NULL,
    return_ev = TRUE,
    return_state = FALSE,
    debug = FALSE,
    ...
) {
  start_time <- proc.time()

  mu   <- params$mu
  K    <- params$K
  beta <- params$beta
  if (!is.finite(mu) || mu <= 0) stop("mu>0 required")
  if (!is.finite(K)  || K  < 0) stop("K>=0 required")
  if (!is.finite(beta) || beta <= 0) stop("beta>0 required")

  # Determine mark type (strict switch)
  if (is.null(mark_type)) {
    if (exists("sim_mark_ba_bip", mode = "function") && identical(sim_mark_fun, sim_mark_ba_bip)) {
      mark_type <- "ba_bip"
    } else if (exists("sim_mark_cs", mode = "function") && identical(sim_mark_fun, sim_mark_cs)) {
      mark_type <- "cs"
    } else {
      mark_type <- "ba"
    }
  }
  mark_type <- match.arg(mark_type, c("ba", "ba_bip", "cs"))

  # CS: model cache for ernm
  cs_model_cache <- NULL
  if (mark_type == "cs") {
    cs_model_cache <- new.env(parent = emptyenv())
    cs_model_cache$rhs <- NULL
    cs_model_cache$model <- NULL
  }

  state <- state_init_fun()
  cache <- time_cache_init(T0 = T0)

  t <- T0

  events_t <- numeric(0)
  accept_probs <- numeric(0)
  nodes_rows <- list()
  edges_rows <- list()

  # compute initial bound M = lambda(t)
  cache <- time_cache_step(cache, t, beta = beta)
  lambda_now <- temporal_intensity_from_cache(cache, mu = mu, K = K, beta = beta)
  M <- lambda_now

  while (t < T_end) {
    if (!is.finite(M) || M <= 0) break

    # propose next candidate time
    dt <- stats::rexp(1, rate = M)
    t_prop <- t + dt
    if (t_prop > T_end) break

    # move cache to candidate time and compute true intensity there
    cache_prop <- time_cache_step(cache, t_prop, beta = beta)
    lambda_prop <- temporal_intensity_from_cache(cache_prop, mu = mu, K = K, beta = beta)

    a <- lambda_prop / M
    accept_probs <- c(accept_probs, a)
    if (debug) message("t=", t_prop, " M=", M, " lambda=", lambda_prop, " a=", a)

    if (stats::runif(1) < a) {
      # accepted: sample mark given *pre-event* state at t_prop
      if (mark_type == "cs") {
        mk <- sim_mark_fun(state = state, t_k = t_prop, params = params, model_cache = cs_model_cache, ...)
      } else {
        mk <- sim_mark_fun(state = state, t_k = t_prop, params = params, delta = delta, ...)
      }

      arrivals_k <- mk$arrivals
      edges_k <- mk$edges

      # Strict BA-bip roles at this accepted time (event first, then perps)
      role_k <- NULL
      if (mark_type == "ba_bip") {
        if (length(arrivals_k) < 1L) stop("BA-bip strict: sim_mark_fun produced no arrivals.")
        role_k <- c("event", rep("perp", length(arrivals_k) - 1L))
      }

      if (length(arrivals_k) > 0) {
        df_nodes <- data.frame(
          id = arrivals_k,
          time = rep(t_prop, length(arrivals_k)),
          stringsAsFactors = FALSE
        )
        if (!is.null(role_k)) df_nodes$role <- role_k
        nodes_rows[[length(nodes_rows) + 1L]] <- df_nodes
      }

      if (!is.null(edges_k) && nrow(edges_k) > 0) {
        edges_rows[[length(edges_rows) + 1L]] <- data.frame(
          i = as.character(edges_k$i),
          j = as.character(edges_k$j),
          time = rep(t_prop, nrow(edges_k)),
          stringsAsFactors = FALSE
        )
      }

      if (is.null(role_k)) {
        state <- state_step_fun(state, arrivals_k, edges_k, t_prop)
      } else {
        state <- state_step_fun(state, arrivals_k, edges_k, t_prop, role_k = role_k)
      }

      # Hawkes jump at an accepted event (event-driven ground process)
      cache_prop$S <- cache_prop$S + 1
      events_t <- c(events_t, t_prop)
    }

    # advance time regardless; update cache + bound for next iteration
    t <- t_prop
    cache <- cache_prop
    M <- lambda_prop
  }

  out <- list(
    events = list(n = length(events_t), t = events_t),
    accept_probs = accept_probs
  )

  if (return_ev) {
    nodes <- if (length(nodes_rows)) do.call(rbind, nodes_rows) else data.frame(id=character(0), time=numeric(0))
    edges <- if (length(edges_rows)) do.call(rbind, edges_rows) else data.frame(i=character(0), j=character(0), time=numeric(0))
    out$ev <- make_events(
      nodes = nodes,
      edges = edges,
      allow_implicit_birth = (mark_type != "ba_bip")
    )
  }
  if (return_state) out$state <- state

  print(paste0("Simulation took ", round((proc.time()-start_time)[3],2), " seconds"))

  out
}
