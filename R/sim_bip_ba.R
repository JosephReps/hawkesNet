# R/sim_bip_ba.R
#' Simulate BA-bip mark (arrivals + edges) for one event
#'
#' Strict BA-bip generative story:
#' - Exactly 1 new EVENT node arrives at each accepted event time.
#' - Additionally, K new PERP nodes arrive, with K ~ Poisson(lambda_new).
#' - Deterministic edges: new_event -- each new_perp.
#' - Bernoulli edges: new_event -- each OLD perp i independently with probability p_i,
#'   where p_i is proportional to (deg_i + delta) * exp(-beta_edges * age_i),
#'   normalized over OLD perp nodes only.
#'
#' This is the network-based version (net stores vertex attrs `time` and `role`).
#' The return format matches sim_mark_ba(): `arrivals` is a data.frame of new nodes,
#' and `edges` is a data.frame of new edges.
#'
#' @param net Pre-event network state.
#' @param t_k Numeric scalar event time.
#' @param params List containing at least beta_edges and lambda_new.
#' @param delta Nonnegative scalar baseline for attractiveness.
#' @param new_event_id Optional character scalar to force the event node id (tests).
#'
#' @return list(arrivals = <data.frame(id, role)>, edges = <data.frame(i,j)>)
#' @export
sim_mark_ba_bip <- function(net, t_k, params, delta = 0.001, new_event_id = NULL) {
  stopifnot(length(t_k) == 1L, is.numeric(t_k), is.finite(t_k))
  stopifnot(is.list(params))
  stopifnot("beta_edges" %in% names(params), "lambda_new" %in% names(params))

  beta_edges <- params$beta_edges
  lambda_new <- params$lambda_new

  stopifnot(length(beta_edges) == 1L, is.numeric(beta_edges), is.finite(beta_edges), beta_edges >= 0)
  stopifnot(length(lambda_new) == 1L, is.numeric(lambda_new), is.finite(lambda_new), lambda_new >= 0)
  stopifnot(length(delta) == 1L, is.numeric(delta), is.finite(delta), delta >= 0)

  # --- choose/generate new EVENT node id ---
  if (is.null(new_event_id)) {
    new_event_id <- next_node_id(net)
  }
  stopifnot(length(new_event_id) == 1L, nzchar(new_event_id))
  new_event_id <- as.character(new_event_id)

  # --- sample number of new PERP nodes ---
  K <- stats::rpois(1L, lambda_new)
  if (!is.finite(K) || K < 0) stop("Invalid K from rpois().", call. = FALSE)

  # --- generate ids for new PERP nodes (sequentially after the event id) ---
  new_perps <- character(0)
  if (K > 0L) {
    # Prefer integer ids if possible; otherwise fall back to n<idx>
    existing <- network::network.vertex.names(net)
    if (is.null(existing)) existing <- character(0)

    # Try integer id scheme
    suppressWarnings(ints <- as.integer(c(existing, new_event_id)))
    ints <- ints[is.finite(ints)]
    if (length(ints) > 0L && is.finite(suppressWarnings(as.integer(new_event_id)))) {
      base <- as.integer(new_event_id)
      new_perps <- as.character(seq.int(from = base + 1L, length.out = K))
    } else if (length(ints) > 0L) {
      base <- max(ints)
      new_perps <- as.character(seq.int(from = base + 1L, length.out = K))
    } else {
      # Non-integer IDs: ensure uniqueness with a simple counter
      n0 <- length(existing)
      new_perps <- paste0("n", seq.int(from = n0 + 2L, length.out = K))  # +1 would be event
      if (new_event_id %in% new_perps) {
        # extremely unlikely, but be safe
        new_perps <- paste0(new_perps, "_p")
      }
    }
  }

  arrivals <- data.frame(
    id = c(new_event_id, new_perps),
    role = c("event", rep("perp", length(new_perps))),
    stringsAsFactors = FALSE
  )

  # --- deterministic edges to NEW perps ---
  det_edges <- if (length(new_perps) > 0L) {
    data.frame(i = rep.int(new_event_id, length(new_perps)), j = new_perps, stringsAsFactors = FALSE)
  } else {
    data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
  }

  # --- Bernoulli edges to OLD perps only ---
  p_old <- edge_probs_ba_bip(net, t_k = t_k, beta_edges = beta_edges, delta = delta)
  old_perps <- names(p_old)

  if (length(old_perps) == 0L) {
    return(list(arrivals = arrivals, edges = det_edges))
  }

  add <- stats::runif(length(p_old)) < p_old
  bern_edges <- if (any(add)) {
    data.frame(i = rep.int(new_event_id, sum(add)), j = old_perps[add], stringsAsFactors = FALSE)
  } else {
    data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
  }

  edges <- rbind(det_edges, bern_edges)
  list(arrivals = arrivals, edges = edges)
}
