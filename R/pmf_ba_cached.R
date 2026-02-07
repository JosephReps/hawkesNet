# Compile BA mark log-likelihood into a cheap function(params) evaluator.
#
# This mirrors your current log_pmf_ba() behaviour:
# - pre-event old nodes have attachment probs p_old
# - for the observed event, "present" old nodes are those that appear in new_edges
# - logp = sum_{present} log(p_old) + sum_{absent} log(1 - p_old)
#
# After compilation, evaluation is all vectorised arithmetic (no network calls).

compile_logpmf_ba <- function(
    events,
    T0 = 0,
    net_init_fun = net_init,
    net0 = NULL,
    net_step_fun = net_add_event,
    delta = 0.001,
    eps = 1e-12
) {
  stopifnot(is.numeric(delta), length(delta) == 1L, is.finite(delta), delta >= 0)
  stopifnot(is.numeric(eps), length(eps) == 1L, is.finite(eps), eps > 0)

  arrivals <- nodes_by_event(events)
  edges_by <- edges_by_event(events)
  event_times <- events$times$t

  # One pass through the data (network manipulation happens ONCE here)
  net <- if (is.null(net0)) net_init_fun() else net0

  cache <- vector("list", length(event_times))

  for (idx in seq_along(event_times)) {
    t_k <- event_times[idx]
    new_nodes <- arrivals[[idx]]
    new_edges <- edges_by[[idx]]

    # Your BA kernel requires exactly 1 arriving node
    if (!is.null(new_nodes) && nrow(new_nodes) != 1L) {
      stop("BA compile: expected exactly 1 arriving node per event.", call. = FALSE)
    }

    old_ids <- network::network.vertex.names(net)
    if (is.null(old_ids)) old_ids <- character(0)

    if (length(old_ids) == 0L) {
      # Must match your current behaviour:
      # - if no existing nodes and no edges: logp = 0
      # - if edges exist: invalid
      if (!is.null(new_edges) && nrow(new_edges) > 0L) {
        stop("BA event has no eligible existing nodes, but edges were observed.", call. = FALSE)
      }
      cache[[idx]] <- list(type = "empty")  # contributes 0
    } else {
      # Precompute quantities that DO NOT depend on beta_edges:
      # degrees + ages at t_k for all old nodes.
      age <- t_k - network::get.vertex.attribute(net, "time")
      deg <- sna::degree(net)

      # Determine which old nodes were attached-to in this event.
      # Your current BA PMF uses:
      # present <- intersect(existing_nodes, c(new_edges$i, new_edges$j))
      present_ids <- character(0)
      if (!is.null(new_edges) && nrow(new_edges) > 0L) {
        present_ids <- intersect(old_ids, c(new_edges$i, new_edges$j))
      }
      present <- old_ids %in% present_ids

      cache[[idx]] <- list(
        type = "ba",
        ids = old_ids,
        deg = as.numeric(deg),
        age = as.numeric(age),
        present = present
      )
    }

    # Advance net along observed history
    net <- net_step_fun(net, new_nodes, new_edges, t_k)
  }

  # Return a cheap evaluator: function(params) -> total mark loglik
  force(cache); force(delta); force(eps)
  function(params) {
    beta_edges <- params[["beta_edges"]]
    if (!(is.numeric(beta_edges) && length(beta_edges) == 1L && is.finite(beta_edges) && beta_edges >= 0)) {
      return(-Inf)
    }

    ll <- 0.0
    for (idx in seq_along(cache)) {
      c <- cache[[idx]]
      if (identical(c$type, "empty")) next

      # weights w_i = (deg_i + delta) * exp(-beta_edges * age_i)
      w <- (c$deg + delta) * exp(-beta_edges * c$age)
      W <- sum(w)

      if (!is.finite(W) || W <= 0) return(-Inf)

      p <- w / W

      # Numerical safety (avoid log(0) / log(negative))
      p <- pmin(pmax(p, eps), 1 - eps)

      # logp = sum_{present} log(p) + sum_{absent} log(1-p)
      ll <- ll + sum(ifelse(c$present, log(p), log1p(-p)))
    }

    as.numeric(ll)
  }
}
