#' Title
#'
#' @param params
#' @param T_end
#' @param T0
#' @param delta
#' @param mark_type
#' @param return_ev
#' @param debug
#' @param net_init_fun
#' @param net_step_fun
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
sim_hawkesNet <- function(
    params, T_end, T0 = 0,
    delta = 0.001,
    mark_type = "ba",
    return_ev = TRUE,
    debug = FALSE,
    net_init_fun = net_init,
    net_step_fun = net_add_event,
    ...
) {
  start_time <- proc.time()

  mu <- params$mu
  K <- params$K
  beta <- params$beta
  if (!is.finite(mu) || mu <= 0) stop("mu>0 required")
  if (!is.finite(K)  || K  < 0) stop("K>=0 required")
  if (!is.finite(beta) || beta <= 0) stop("beta>0 required")

  # Mark kernel dispatch
  if (mark_type == "ba") {
    sim_mark_fun <- sim_mark_ba
  } else if (mark_type == "ba_bip") {
    sim_mark_fun <- sim_mark_ba_bip
  } else if (mark_type == "cs") {
    sim_mark_fun <- sim_mark_cs
  } else if (mark_type == "cs_bip") {
    sim_mark_fun <- sim_mark_cs_bip
  } else {
    stop("mark_type must be one of 'ba', 'cs', 'ba_bip', 'cs_bip'")
  }


  # CS: model cache for ernm
  cs_model_cache <- NULL
  if (mark_type %in% c("cs", "cs_bip")) {
    cs_model_cache <- list(rhs = NULL, model = NULL)
  }

  net <- net_init_fun()
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

    accepted <- FALSE
    a <- lambda_prop / M
    accept_probs <- c(accept_probs, a)
    if (debug) message("t=", t_prop, " M=", M, " lambda=", lambda_prop, " a=", a)

    if (stats::runif(1) < a) {
      accepted <- TRUE
      # accepted: sample mark given *pre-event* state at t_prop
      if (mark_type %in% c("cs", "cs_bip")) {
        mk <- sim_mark_fun(net = net, t_k = t_prop, params = params, model_cache = cs_model_cache, ...)
      } else {
        mk <- sim_mark_fun(net = net, t_k = t_prop, params = params, delta = delta, ...)
      }

      # CS: carry forward ERNM model cache (to avoid rebuilding each event)
      if (mark_type %in% c("cs", "cs_bip") && !is.null(mk$model_cache)) {
        cs_model_cache <- mk$model_cache
      }

      arrivals_k <- mk$arrivals
      edges_k <- mk$edges

      # Normalise arrivals to a data.frame with at least `id`
      if (is.null(arrivals_k)) {
        arrivals_k <- data.frame(id = character(0), stringsAsFactors = FALSE)
      } else if (is.character(arrivals_k)) {
        arrivals_k <- data.frame(id = as.character(arrivals_k), stringsAsFactors = FALSE)
      }
      if (!is.data.frame(arrivals_k) || !("id" %in% names(arrivals_k))) {
        stop("sim_mark_fun must return `arrivals` as a data.frame with column `id`.")
      }

      # BA-bip: ensure roles exist (event first, then perps) if not provided by sim_mark_fun
      if (mark_type == "ba_bip" && !("role" %in% names(arrivals_k))) {
        if (nrow(arrivals_k) < 1L) stop("BA-bip strict: sim_mark_fun produced no arrivals.")
        arrivals_k$role <- c("event", rep("perp", nrow(arrivals_k) - 1L))
      }

      if (nrow(arrivals_k) > 0) {
        df_nodes <- arrivals_k
        df_nodes$time <- rep(t_prop, nrow(arrivals_k))
        nodes_rows[[length(nodes_rows) + 1L]] <- df_nodes
      }

      if (!is.null(edges_k) && nrow(edges_k) > 0) {
        df_edges <- data.frame(
          i = as.character(edges_k$i),
          j = as.character(edges_k$j),
          time = rep(t_prop, nrow(edges_k)),
          stringsAsFactors = FALSE
        )
        edges_rows[[length(edges_rows) + 1L]] <- df_edges
      }

      net <- net_step_fun(net, arrivals_k, edges_k, t_prop)

      # Hawkes jump at an accepted event (event-driven ground process)
      cache_prop$S <- cache_prop$S + 1
      events_t <- c(events_t, t_prop)

      lambda_after <- temporal_intensity_from_cache(cache_prop, mu = mu, K = K, beta = beta)
    }

    # advance time regardless; update cache + bound for next iteration
    t <- t_prop
    cache <- cache_prop
    if (accepted) {
      M <- lambda_after
    } else {
      M <- lambda_prop
    }
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

  out$net <- net

  print(paste0("Simulation took ", round((proc.time()-start_time)[3],2), " seconds"))

  out
}
