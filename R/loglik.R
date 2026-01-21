# loglik.R

# Internal: initialize recursion cache for an exponential Hawkes kernel.
# (unchanged)
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
  mu * dt + (K / beta) * S_post_prev * (1 - exp(-beta*dt))
}

# Internal: split roles per event, aligned with nodes_by_event ordering
roles_by_event <- function(events) {
  validate_events(events)
  n <- nrow(events$times)
  out <- replicate(n, character(0), simplify = FALSE)

  nd <- events$nodes
  if (nrow(nd) == 0L) return(out)

  if (!("role" %in% names(nd))) stop("BA-bip strict: ev$nodes must include a `role` column.", call. = FALSE)

  spl <- split(seq_len(nrow(nd)), nd$event_id)
  idx <- as.integer(names(spl))

  out[idx] <- lapply(spl, function(ii) {
    id <- as.character(nd$id[ii])
    rr <- as.character(nd$role[ii])

    keep <- !is.na(id) & nzchar(id)
    id <- id[keep]
    rr <- rr[keep]

    if (length(id) != length(rr)) stop("Internal error: role/id length mismatch after filtering.", call. = FALSE)
    if (anyNA(rr)) stop("BA-bip strict: node roles cannot be NA.", call. = FALSE)

    ok <- rr %in% c("event", "perp")
    if (!all(ok)) stop("BA-bip strict: roles must be 'event' or 'perp'.", call. = FALSE)

    rr
  })

  out
}

loglik <- function(
    ev, params, T_end = NULL, T0 = 0,
    mark_logpmf = log_pmf_ba,
    mark_type = NULL,
    state_init_fun = state_init, state0 = NULL,
    state_step_fun = state_step,
    debug = FALSE,
    ...
) {

  mu <- params$mu
  K <- params$K
  beta <- params$beta

  if (!is.finite(mu) || !is.finite(K) || !is.finite(beta) || mu <= 0 || K < 0 || beta <= 0)
    return(-Inf)

  event_times <- ev$times$t
  if (length(event_times) == 0L) return(0)

  # Define window end
  if (is.null(T_end)) T_end <- max(event_times)
  if (T_end < max(event_times)) stop("T_end must be >= max event time")

  # Determine mark type (strict switch)
  if (is.null(mark_type)) {
    if (exists("log_pmf_ba_bip", mode = "function") && identical(mark_logpmf, log_pmf_ba_bip)) {
      mark_type <- "ba_bip"
    } else if (exists("log_pmf_cs", mode = "function") && identical(mark_logpmf, log_pmf_cs)) {
      mark_type <- "cs"
    } else {
      mark_type <- "ba"
    }
  }
  mark_type <- match.arg(mark_type, c("ba", "ba_bip", "cs"))

  arrivals <- nodes_by_event(ev)
  edges <- edges_by_event(ev)

  roles <- NULL
  if (mark_type == "ba_bip") {
    roles <- roles_by_event(ev)

    # Extra strict sanity: every event must have exactly one 'event' arrival
    for (k in seq_along(arrivals)) {
      if (length(arrivals[[k]]) == 0L) stop("BA-bip strict: event ", k, " has zero node arrivals.")
      rk <- roles[[k]]
      if (length(rk) != length(arrivals[[k]])) stop("BA-bip strict: role_k length mismatch at event ", k, ".")
      if (sum(rk == "event") != 1L) stop("BA-bip strict: event ", k, " must have exactly one 'event' arrival.")
      if (!all(rk %in% c("event", "perp"))) stop("BA-bip strict: invalid roles at event ", k, ".")
    }
  }

  ll <- 0.0

  state <- if (is.null(state0)) state_init_fun() else state0

  time_cache <- time_cache_init(T0 = T0)

  # For compensator
  t_prev <- T0
  S_post_prev <- 0.0
  compensator <- 0.0

  # CS: build a model cache once (ernm model object can be reused)
  cs_model_cache <- NULL
  if (mark_type == "cs") {
    cs_model_cache <- new.env(parent = emptyenv())
    cs_model_cache$rhs <- NULL
    cs_model_cache$model <- NULL
  }

  # Iterate over events
  for (idx in seq_along(event_times)) {
    t_k <- event_times[idx]

    arrivals_k <- arrivals[[idx]]
    edges_k <- edges[[idx]]

    role_k <- NULL
    if (mark_type == "ba_bip") role_k <- roles[[idx]]

    # (1) Compensator over (t_prev, t_k)
    compensator <- compensator + compensator_inc_exp(
      dt = t_k - t_prev, mu = mu, K = K, beta = beta, S_post_prev = S_post_prev
    )

    # (2) Decay cache to t_k => S_pre at t_k
    time_cache <- time_cache_step(time_cache, t_k, beta = beta)

    # (3) Temporal intensity (strictly past events)
    lambda_g <- temporal_intensity_from_cache(time_cache, mu = mu, K = K, beta = beta)
    if (!is.finite(lambda_g) || lambda_g <= 0) return(-Inf)

    # (4) Mark log-pmf using pre-event network state
    if (mark_type == "ba_bip") {
      mark_ll <- mark_logpmf(state, arrivals_k, edges_k, t_k, params, role_k = role_k, ...)
    } else if (mark_type == "cs") {
      # pass model cache so ernm model is built once
      mark_ll <- mark_logpmf(state, arrivals_k, edges_k, t_k, params, model_cache = cs_model_cache, ...)
    } else {
      mark_ll <- mark_logpmf(state, arrivals_k, edges_k, t_k, params, ...)
    }
    if (!is.finite(mark_ll)) return(-Inf)

    # (5) Update log-likelihood
    ll <- ll + mark_ll + log(lambda_g)

    if (debug) {
      message("Event at time: t = ", t_k, "\nlog q = ",
              mark_ll, "\nlog lambda_g = ", log(lambda_g), "\nll = ", ll, "\n")
    }

    # (6) Update state
    if (is.null(role_k)) {
      state <- state_step_fun(state, arrivals_k, edges_k, t_k)
    } else {
      state <- state_step_fun(state, arrivals_k, edges_k, t_k, role_k = role_k)
    }

    # (7) Record event for future, S_post = S_pre + 1
    time_cache$S <- time_cache$S + 1

    # Update trackers for next compensator gap
    t_prev <- t_k
    S_post_prev <- time_cache$S
  }

  # Tail compensator (t_last, T_end]
  compensator <- compensator + compensator_inc_exp(
    dt = T_end - t_prev, mu = mu, K = K, beta = beta, S_post_prev = S_post_prev
  )

  ll - compensator
}
