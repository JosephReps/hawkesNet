# sim_ba_bip.R

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
#' Node roles are NOT stored in the return object (to match sim_mark_ba()).
#' In your Hawkes simulator, you should pass the roles into state_step() via role_k.
#'
#' @param state Pre-event state (must include named deg, born, role).
#' @param t_k Numeric scalar event time.
#' @param params List containing at least beta_edges and lambda_new.
#' @param delta Nonnegative scalar baseline for attractiveness.
#' @param new_event Optional character scalar to force the event node id (tests).
#'
#' @return list(arrivals = <chr vec>, edges = data.frame(i,j))
#' @export
sim_mark_ba_bip <- function(state, t_k, params, delta = 0.001, new_event = NULL) {
  stopifnot(is.list(state), all(c("deg", "born", "role") %in% names(state)))
  stopifnot(length(t_k) == 1L, is.numeric(t_k), is.finite(t_k))
  stopifnot(is.list(params))

  stopifnot("beta_edges" %in% names(params))
  stopifnot("lambda_new" %in% names(params))
  beta_edges <- params$beta_edges
  lambda_new <- params$lambda_new

  stopifnot(length(beta_edges) == 1L, is.numeric(beta_edges), is.finite(beta_edges), beta_edges >= 0)
  stopifnot(length(lambda_new) == 1L, is.numeric(lambda_new), is.finite(lambda_new), lambda_new >= 0)
  stopifnot(length(delta) == 1L, is.numeric(delta), is.finite(delta), delta >= 0)

  # --- choose/generate new EVENT node id ---
  if (is.null(new_event)) {
    new_event <- next_node_id(state)  # reuse existing helper
  }
  stopifnot(length(new_event) == 1L, nzchar(new_event))
  new_event <- as.character(new_event)

  # --- sample number of new PERP nodes ---
  K <- stats::rpois(1L, lambda_new)
  if (!is.finite(K) || K < 0) stop("Invalid K from rpois().", call. = FALSE)

  # --- choose/generate ids for new PERP nodes (sequentially after the event id) ---
  # Use the same id-space; ensure uniqueness by building a temporary state for id generation.
  new_perps <- character(0)
  if (K > 0L) {
    tmp <- state
    # pretend the new event exists so next_node_id increments beyond it
    tmp$deg[new_event] <- 0L
    tmp$born[new_event] <- t_k
    tmp$adj[new_event] <- list(character(0))
    tmp$role[new_event] <- "event"

    new_perps <- character(K)
    for (m in seq_len(K)) {
      nid <- next_node_id(tmp)
      new_perps[m] <- nid
      tmp$deg[nid] <- 0L
      tmp$born[nid] <- t_k
      tmp$adj[nid] <- list(character(0))
      tmp$role[nid] <- "perp"
    }
  }

  arrivals <- c(new_event, new_perps)

  # --- deterministic edges to new perps ---
  det_edges <- if (K > 0L) {
    data.frame(i = rep.int(new_event, K), j = new_perps, stringsAsFactors = FALSE)
  } else {
    data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
  }

  # --- Bernoulli edges to OLD perps only ---
  p_old <- edge_probs_ba_bip(state, t_k = t_k, beta_edges = beta_edges, delta = delta)
  old_perps <- names(p_old)

  if (length(old_perps) == 0L) {
    # No old perps to attach to: return only deterministic edges
    return(list(arrivals = arrivals, edges = det_edges))
  }

  add <- stats::runif(length(p_old)) < p_old
  bern_edges <- if (any(add)) {
    data.frame(i = rep.int(new_event, sum(add)), j = old_perps[add], stringsAsFactors = FALSE)
  } else {
    data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
  }

  edges <- rbind(det_edges, bern_edges)
  list(arrivals = arrivals, edges = edges)
}
