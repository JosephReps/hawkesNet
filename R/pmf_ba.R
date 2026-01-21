# R/pmf_ba.R

# Internal: compute BA attachment probabilities over existing (old) nodes.
# Returns a named numeric vector p_old with names = old node IDs.
edge_probs_ba <- function(state, t_k, beta_edges = 0, delta = 0.001) {
  old_nodes <- names(state$deg)
  if (is.null(old_nodes)) old_nodes <- character(0)

  if (length(old_nodes) == 0L) {
    return(structure(numeric(0), names = character(0)))
  }

  # Birth-time age of old nodes
  age <- t_k - state$born[old_nodes]
  if (any(!is.finite(age))) stop("Non-finite ages encountered in born times.", call. = FALSE)

  # Degrees of old nodes
  deg <- as.numeric(state$deg[old_nodes])
  if (any(!is.finite(deg))) stop("Non-finite degrees encountered.", call. = FALSE)

  # Attachment attractiveness score (must be >= 0)
  #   delta=0 results in pure degree-based preference
  w <- (deg + delta) * exp(-beta_edges * age)
  if (any(!is.finite(w))) stop("Non-finite BA weights encountered.", call. = FALSE)
  if (any(w < 0)) stop("Invalid BA weights: must be nonnegative.", call. = FALSE)

  W <- sum(w)
  if (!is.finite(W) || W <= 0) {
    stop("Invalid BA weights: sum(weights) must be > 0.", call. = FALSE)
  }

  # Normalise probabilities
  #   IMPORTANT: if only 1 old node, this implies probability of attachment = 1)
  p_old <- w / W
  names(p_old) <- old_nodes

  # Numerical safety
  p_old <- pmin(pmax(p_old, 0), 1)

  p_old
}

#' Log PMF of the BA mark kernel (Bernoulli-product form) for one event
#'
#' Computes the log probability of the observed edge set `edges_k` under the
#' Barabási–Albert (BA) **Bernoulli-product** mark model, for a single event at
#' time `t_k`. See `sim_mark_ba()` for the corresponding simulator.
#'
#' @param state Pre-event state (see `state_init()`). Must include named `deg` and `born`.
#' @param arrival_k Node ID of the node arriving at time `t_k`.
#' @param edges_k A data.frame with columns `i` and `j` giving observed undirected
#'   edges at time `t_k`. May be a 0-row data.frame.
#' @param t_k Numeric scalar event time.
#' @param beta_edges Nonnegative scalar controlling exponential ageing in the
#'   attachment weights.
#'
#' @return A numeric scalar: log probability of `edges_k` under the BA mark model.
#' @export
log_pmf_ba <- function(state,
                       arrival_k,
                       edges_k,
                       t_k,
                       params,
                       delta = 0.001) {
  stopifnot("beta_edges" %in% names(params))
  beta_edges <- params[["beta_edges"]]

  # if (beta_edges < 0) browser()

  stopifnot(is.list(state), all(c("deg", "born") %in% names(state)))
  stopifnot(is.data.frame(edges_k), all(c("i", "j") %in% names(edges_k)))
  stopifnot(length(t_k) == 1L, is.numeric(t_k), is.finite(t_k))
  stopifnot(length(beta_edges) == 1L, is.numeric(beta_edges),
            is.finite(beta_edges), beta_edges >= 0)
  stopifnot(length(delta) == 1L, is.numeric(delta), is.finite(delta), delta >= 0)

  # BA kernel requires exactly 1 arriving node
  stopifnot(length(arrival_k) == 1L, nzchar(arrival_k))
  new_node <- arrival_k[[1]]

  # Pre-event old nodes
  old_nodes <- names(state$deg)
  if (is.null(old_nodes)) old_nodes <- character(0)

  # No old nodes: no eligible attachment targets
  if (length(old_nodes) == 0L) {
    if (nrow(edges_k) == 0L) return(0)
    stop("BA event has no eligible existing nodes, but edges were observed.", call. = FALSE)
  }

  # Old-node attachment probabilities
  p_old <- edge_probs_ba(state, t_k = t_k, beta_edges = beta_edges, delta = delta)
  old_nodes <- names(p_old)

  # No edges present in the event
  if (nrow(edges_k) == 0L) {
    return(as.numeric(sum(log(1 - p_old))))
  }

  # Otherwise, validate & orient observed edges as (new_node -> old_node)
  ends <- validate_edges_ba(edges_k, arrival_k = new_node, old_nodes = old_nodes)
  old_end <- ends$old_end

  present <- old_end
  absent  <- setdiff(old_nodes, present)

  lp <- 0.0
  if (length(present) > 0L) lp <- lp + sum(log(p_old[present]))
  if (length(absent) > 0L) lp <- lp + sum(log(1 - p_old[absent]))

  as.numeric(lp)
}
