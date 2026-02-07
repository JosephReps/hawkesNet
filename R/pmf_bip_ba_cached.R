# R/pmf_bip_ba_cached.R
#
# Compile BA-bip mark log-likelihood into a cheap function(params) evaluator.
#
# Matches log_pmf_ba_bip() semantics:
# - Exactly one new EVENT node per event
# - K_new new PERP nodes per event; Poisson(K_new | lambda_new)
# - Deterministic edges: new_event -- each new_perp (prob 1)
# - Bernoulli-product edges: new_event -- each OLD perp, with
#     p_old ∝ (deg_old + delta) * exp(-beta_edges * age_old)
#
# After compilation, evaluation is arithmetic only (no network calls).

compile_logpmf_ba_bip <- function(
    events,
    T0 = 0,
    net_init_fun = net_init,
    net0 = NULL,
    net_step_fun = net_add_event,
    delta = 0.001,
    eps = 1e-12
) {
  stopifnot(is.numeric(delta), length(delta) == 1L, is.finite(delta), delta >= 0)
  stopifnot(is.numeric(eps), length(eps) == 1L, is.finite(eps), eps > 0, eps < 0.5)

  arrivals <- nodes_by_event(events)
  edges_by <- edges_by_event(events)
  event_times <- .get_event_times(events)

  net <- if (is.null(net0)) net_init_fun() else net0
  cache <- vector("list", length(event_times))

  for (k in seq_along(event_times)) {
    t_k <- event_times[k]
    new_nodes <- arrivals[[k]]
    new_edges <- edges_by[[k]]

    # Normalise new_nodes / new_edges
    if (is.null(new_nodes)) new_nodes <- data.frame(id = character(0), stringsAsFactors = FALSE)
    stopifnot(is.data.frame(new_nodes), "id" %in% names(new_nodes))

    if (is.null(new_edges)) {
      edges_k <- data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
    } else {
      edges_k <- new_edges
      stopifnot(is.data.frame(edges_k), all(c("i","j") %in% names(edges_k)))
      if (nrow(edges_k) == 0L) edges_k <- data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
    }

    if (!("role" %in% names(new_nodes))) {
      stop("BA-bip compile: expected `role` column in new_nodes.", call. = FALSE)
    }
    new_nodes$role <- as.character(new_nodes$role)

    is_event <- !is.na(new_nodes$role) & new_nodes$role == "event"
    is_perp  <- !is.na(new_nodes$role) & new_nodes$role == "perp"

    if (sum(is_event) != 1L) {
      cache[[k]] <- list(type = "invalid", reason = "not_exactly_one_event")
      net <- net_step_fun(net, new_nodes, new_edges, t_k)
      next
    }

    new_event <- as.character(new_nodes$id[which(is_event)[1]])
    new_perps <- as.character(new_nodes$id[which(is_perp)])

    K_new <- length(new_perps)

    # Deterministic edges: new_event -- each new_perp must be present if K_new > 0
    if (K_new > 0L) {
      edge_keys <- .edge_keys_from_df_undirected(edges_k)
      need <- .edge_key_undirected(rep.int(new_event, K_new), new_perps)
      if (!all(need %in% edge_keys)) {
        cache[[k]] <- list(type = "invalid", reason = "missing_event_newperp_edges")
        net <- net_step_fun(net, new_nodes, new_edges, t_k)
        next
      }
    }

    # Old perp nodes only
    roles <- network::get.vertex.attribute(net, "role")
    ids <- network::network.vertex.names(net)
    if (is.null(ids)) ids <- as.character(seq_len(net %n% "n"))

    old_perp_idx <- which(!is.na(roles) & roles == "perp")
    if (length(old_perp_idx) == 0L) {
      cache[[k]] <- list(type = "no_old_perps", K_new = K_new)
      net <- net_step_fun(net, new_nodes, new_edges, t_k)
      next
    }

    born <- network::get.vertex.attribute(net, "time")
    if (is.null(born)) stop("BA-bip compile requires vertex attribute `time`.", call. = FALSE)

    old_perp_ids <- ids[old_perp_idx]
    age <- t_k - born[old_perp_idx]
    deg <- sna::degree(net)[old_perp_idx]

    # Which old perps are connected to the new event in this mark?
    edge_keys <- .edge_keys_from_df_undirected(edges_k)
    is_obs <- .edge_key_undirected(rep.int(new_event, length(old_perp_ids)), old_perp_ids) %in% edge_keys

    cache[[k]] <- list(
      type = "ba_bip",
      K_new = K_new,
      old_perp_ids = old_perp_ids,
      deg = as.numeric(deg),
      age = as.numeric(age),
      is_obs = is_obs
    )

    net <- net_step_fun(net, new_nodes, new_edges, t_k)
  }

  force(cache); force(delta); force(eps)

  function(params) {
    beta_edges <- params[["beta_edges"]]
    lambda_new <- params[["lambda_new"]]

    if (!(is.numeric(beta_edges) && length(beta_edges) == 1L && is.finite(beta_edges) && beta_edges >= 0)) return(-Inf)
    if (!(is.numeric(lambda_new) && length(lambda_new) == 1L && is.finite(lambda_new) && lambda_new >= 0)) return(-Inf)

    ll <- 0.0
    for (k in seq_along(cache)) {
      c <- cache[[k]]
      if (identical(c$type, "invalid")) return(-Inf)

      # Poisson term for number of new perps
      ll <- ll + stats::dpois(c$K_new, lambda = lambda_new, log = TRUE)

      if (!identical(c$type, "ba_bip")) next

      w <- (c$deg + delta) * exp(-beta_edges * c$age)
      W <- sum(w)
      if (!is.finite(W) || W <= 0) return(-Inf)

      p <- w / W
      p <- .clamp_prob(p, eps = eps)

      ll <- ll + .bern_ll(p, c$is_obs)
    }

    as.numeric(ll)
  }
}
