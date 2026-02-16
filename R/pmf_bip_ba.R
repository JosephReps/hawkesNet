# R/pmf_bip_ba.R

# Internal: undirected edge key helper (order-invariant)
.edge_key <- function(a, b) {
  a <- as.character(a); b <- as.character(b)
  lo <- ifelse(a <= b, a, b)
  hi <- ifelse(a <= b, b, a)
  paste0(lo, "|", hi)
}

.edge_keys_from_df <- function(edges_df) {
  if (is.null(edges_df) || nrow(edges_df) == 0L) return(character(0))
  stopifnot(is.data.frame(edges_df), all(c("i", "j") %in% names(edges_df)))
  .edge_key(edges_df$i, edges_df$j)
}

# Internal: compute BA-bip attachment probabilities over OLD perp nodes only.
# Returns a named numeric vector p_old_perp with names = perp node IDs.
edge_probs_ba_bip <- function(net, t_k, beta_edges = 0.5, delta = 0.001) {
  if (net %n% "n" == 0) return(numeric(0))

  roles <- network::get.vertex.attribute(net, "role")
  ids <- network::network.vertex.names(net)
  if (is.null(ids)) ids <- as.character(seq_len(net %n% "n"))

  old_perp_idx <- which(!is.na(roles) & roles == "perp")
  if (length(old_perp_idx) == 0L) return(numeric(0))

  born <- network::get.vertex.attribute(net, "time")
  if (is.null(born)) {
    stop("BA-bip requires vertex attribute `time` (birth time).", call. = FALSE)
  }

  # Age and degree for OLD perps
  age <- t_k - born[old_perp_idx]
  deg <- sna::degree(net)[old_perp_idx]

  w <- (deg + delta) * exp(-beta_edges * age)
  W <- sum(w)

  if (!is.finite(W) || W <= 0) {
    stop("Invalid BA-bip weights over old perps: sum(weights) must be > 0.", call. = FALSE)
  }

  p <- w / W
  names(p) <- ids[old_perp_idx]

  # Numerical safety
  p <- pmin(pmax(p, 0), 1)
  p
}

#' Log pmf for BA-bip mark (network-based)
#'
#' Likelihood contribution for one accepted event:
#' - Poisson(K_new | lambda_new) for the number of *new perp* nodes in this event.
#' - Deterministic edges between the new event node and each new perp (probability 1).
#' - Independent Bernoulli edges between the new event node and each *old perp* in the
#'   current network, with probabilities from edge_probs_ba_bip().
#'
#' Notes:
#' - We intentionally ignore any other edges in `new_edges` (e.g., event-event or perp-perp),
#'   but we DO enforce the deterministic event--new_perp edges. If they are missing, we return -Inf.
#' - This mirrors the "thinned" BA kernel style you used: minimal validation, but enough to keep
#'   the model identifiable/correct.
#'
#' @param net Pre-event network state.
#' @param new_nodes data.frame of nodes arriving at t_k. Must include `id` and (for BA-bip) `role`.
#' @param new_edges data.frame of edges arriving at t_k. Must include `i`,`j`.
#' @param t_k Numeric scalar event time.
#' @param params List containing at least beta_edges and lambda_new.
#' @param role_k Optional character vector of roles aligned to rows of new_nodes (legacy hook).
#' @param delta Nonnegative scalar baseline for attractiveness.
#'
#' @return numeric scalar log pmf.
#' @export
log_pmf_ba_bip <- function(net,
                           new_nodes,
                           new_edges,
                           t_k,
                           params,
                           delta = 0.001) {
  stopifnot(is.list(params))
  stopifnot(all(c("beta_edges", "lambda_new") %in% names(params)))

  beta_edges <- params$beta_edges
  lambda_new <- params$lambda_new

  stopifnot(length(beta_edges) == 1L, is.numeric(beta_edges), is.finite(beta_edges), beta_edges >= 0)
  stopifnot(length(lambda_new) == 1L, is.numeric(lambda_new), is.finite(lambda_new), lambda_new >= 0)
  stopifnot(length(t_k) == 1L, is.numeric(t_k), is.finite(t_k))
  stopifnot(length(delta) == 1L, is.numeric(delta), is.finite(delta), delta >= 0)

  stopifnot(is.data.frame(new_nodes), "id" %in% names(new_nodes))
  stopifnot(is.data.frame(new_edges), all(c("i", "j") %in% names(new_edges)))

  # Roles: prefer column in new_nodes, otherwise use role_k (legacy)
  if (!("role" %in% names(new_nodes))) {
    stop("BA-bip: expected `role` column in new_nodes.", call. = FALSE)
  } else {
    # Make sure its a character vector (is this really necessary??)
    new_nodes$role <- as.character(new_nodes$role)
  }

  # Identify new event node + new perps
  is_event <- !is.na(new_nodes$role) & new_nodes$role == "event"
  is_perp  <- !is.na(new_nodes$role) & new_nodes$role == "perp"

  if (sum(is_event) != 1L) stop("BA-bip strict: expected exactly one new EVENT node.", call. = FALSE)

  new_event <- as.character(new_nodes$id[which(is_event)[1]])
  new_perps <- as.character(new_nodes$id[which(is_perp)])

  # Poisson term for number of new perps
  K_new <- length(new_perps)
  lp <- stats::dpois(K_new, lambda_new, log = TRUE)

  # Deterministic edges: new_event - each new_perp (probability 1)
  if (K_new > 0L) {
    edge_keys <- .edge_keys_from_df(new_edges)
    need <- .edge_key(rep.int(new_event, K_new), new_perps)
    ok <- need %in% edge_keys
    if (!all(ok)) return(list(logp = -Inf, edge_probs = p_old))
  }

  # Bernoulli-product over OLD perps only
  p_old <- edge_probs_ba_bip(net, t_k = t_k, beta_edges = beta_edges, delta = delta)
  old_perps <- names(p_old)
  if (length(old_perps) == 0L) return(list(logp = as.numeric(lp), edge_probs = p_old))

  # Which old perps are connected to the new event in this mark?
  edge_keys <- .edge_keys_from_df(new_edges)
  present <- old_perps[.edge_key(rep.int(new_event, length(old_perps)), old_perps) %in% edge_keys]
  absent <- setdiff(old_perps, present)

  if (length(present) > 0L) lp <- lp + sum(log(p_old[present]))
  if (length(absent)  > 0L) lp <- lp + sum(log1p(-p_old[absent]))


  return(list(logp = as.numeric(lp), edge_probs = p_old))
}
