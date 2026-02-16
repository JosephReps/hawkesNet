# R/pmf_ba.R

# Internal: compute BA attachment probabilities over existing (old) nodes.
# Returns a named numeric vector p_old with names = old node IDs.
edge_probs_ba <- function(net, t_k, beta_edges = 0.5, delta = 0.001) {
  if (net %n% 'n' == 0) {
    return(numeric(0))
  }

  # Birth-time age of old nodes
  age <- t_k - network::get.vertex.attribute(net, "time")
  # Degrees of old nodes
  deg <- sna::degree(net)

  # Attachment attractiveness score (must be >= 0)
  #   delta=0 results in pure degree-based preference
  w <- (deg + delta) * exp(-beta_edges * age)

  # Normalizing term
  W <- sum(w)
  if (!is.finite(W) || W <= 0) {e
    stop("Invalid BA weights: sum(weights) must be > 0.", call. = FALSE)
  }

  # Normalise probabilities
  #   IMPORTANT: if only 1 old node, this implies probability of attachment = 1)
  p_old <- w / W
  names(p_old) <- network::network.vertex.names(net)

  # Numerical safety
  p_old <- pmin(pmax(p_old, 0), 1)

  p_old
}

# Internal: validate that edges_k are BA-valid and return oriented endpoints.
#
# Paper BA (strict):
# - Exactly one arriving node `arrival_k` (character scalar).
# - Every edge must be between that new node and one old node in `old_nodes`.
# - No endpoints outside {arrival_k} ∪ old_nodes.
# - No new-new or old-old edges.
# - No duplicate (new, old) pairs within the event.
#
# Returns:
#   list(new_end = <chr>, old_end = <chr>)
# one element per row of edges_k (same order).
validate_edges_ba <- function(edges_k, arrival_k, old_nodes) {
  stopifnot(is.data.frame(edges_k), all(c("i", "j") %in% names(edges_k)))

  if (nrow(edges_k) == 0L) {
    return(list(new_end = character(0), old_end = character(0)))
  }

  # Exactly one new node for BA
  stopifnot(length(arrival_k) == 1L, nzchar(arrival_k))
  new_node <- as.character(arrival_k[[1]])

  old_nodes <- as.character(old_nodes)
  # if old_nodes is empty, any edge is impossible
  if (length(old_nodes) == 0L) {
    stop("BA event has edges but there are no eligible old nodes.", call. = FALSE)
  }

  i <- as.character(edges_k$i)
  j <- as.character(edges_k$j)

  i_is_new <- (i == new_node)
  j_is_new <- (j == new_node)
  i_is_old <- (i %in% old_nodes)
  j_is_old <- (j %in% old_nodes)

  # Any endpoint must be either the new node or an eligible old node
  endpoint_ok <- (i_is_new | i_is_old) & (j_is_new | j_is_old)
  if (any(!endpoint_ok)) {
    stop("BA event has edge endpoint not in {arrival_k} ∪ old_nodes.", call. = FALSE)
  }

  # Edges must be new-old: exactly one endpoint is the new node, the other is old
  is_new_old <- (i_is_new & j_is_old) | (j_is_new & i_is_old)
  if (any(!is_new_old)) {
    stop("BA event has non-BA edges (must be between the new node and an old node).", call. = FALSE)
  }

  # Orient edges: new_end is always the new node; old_end is the old node
  new_end <- ifelse(i_is_new, i, j)
  old_end <- ifelse(i_is_old, i, j)

  # Sanity: all new_end should equal new_node (given strict checks)
  # (keep as stopifnot for internal consistency)
  stopifnot(all(new_end == new_node))

  # No duplicates: Bernoulli model can't contain repeated (new, old) within one event
  if (anyDuplicated(old_end) > 0L) {
    stop("Duplicate BA edges detected within event (same old node repeated).", call. = FALSE)
  }

  list(new_end = new_end, old_end = old_end)
}

#' Log PMF of the BA mark kernel for one event
#'
#' Computes the log probability of the observed new–old edge set (`new_edges`)
#' under the Barabási–Albert (BA) Bernoulli-product mark model at event time
#' `t_k`. Returns both the per-old-node attachment probabilities and the log PMF.
#'
#' @inheritParams new-event-params
#' @param params Named list of model parameters. Must contain `beta_edges`.
#' @param delta Nonnegative scalar added for numerical stability in attachment
#'   probabilities (passed to `edge_probs_ba()`).
#'
#' @return A list with:
#' \describe{
#'   \item{edge_probs}{Named numeric vector of attachment probabilities for each
#'   eligible existing (old) node.}
#'   \item{logp}{Scalar log probability of the observed edges for this event.}
#' }
#'
#' @details
#' The BA mark model assumes exactly one arriving node in `new_nodes`. Edges in
#' `new_edges` are treated as new–old attachments; if there are no existing nodes
#' in `net`, the only valid observation is `new_edges` empty.
#'
#' @export
log_pmf_ba <- function(net,
                       new_nodes,
                       new_edges,
                       t_k,
                       params,
                       delta = 0.001) {
  # Validate input
  stopifnot("beta_edges" %in% names(params))
  beta_edges <- params[["beta_edges"]]
  stopifnot(is.data.frame(new_edges), all(c("i", "j") %in% names(new_edges)))
  stopifnot(length(t_k) == 1L, is.numeric(t_k), is.finite(t_k))
  stopifnot(length(beta_edges) == 1L, is.numeric(beta_edges),
            is.finite(beta_edges), beta_edges >= 0)
  stopifnot(length(delta) == 1L, is.numeric(delta), is.finite(delta), delta >= 0)

  # BA kernel requires exactly 1 arriving node
  stopifnot(nrow(new_nodes) == 1)
  # new_node <- arrival_k[[1]]

  # Pre-event old nodes
  existing_nodes <- network::network.vertex.names(net)
  if (is.null(existing_nodes)) existing_nodes <- character(0)

  # No existing nodes: no eligible attachment targets
  if (length(existing_nodes) == 0L) {
    if (nrow(new_edges) == 0L) {
      # Have to make sure it is named appropriately, R is pretty clapped sometimes
      return(list(edge_probs = setNames(0, new_nodes$id), logp = 0))
    }
    stop("BA event has no eligible existing nodes, but edges were observed.", call. = FALSE)
  }

  # Old-node attachment probabilities
  p_old <- edge_probs_ba(net, t_k = t_k, beta_edges = beta_edges, delta = delta)

  # No edges present in the event
  if (nrow(new_edges) == 0L) {
    return(list(edge_probs = p_old, logp = as.numeric(sum(log(1 - p_old)))))
  }

  # Otherwise, validate & orient observed edges as (new_node -> old_node)
  # ends <- validate_edges_ba(edges_k, arrival_k = new_node, old_nodes = old_nodes)
  # old_end <- ends$old_end

  # Only 1 new node for BA kernel
  new_node_id <- new_nodes$id
  # Janky, will clean up when I've had proper sleep
  present <- intersect(existing_nodes, c(new_edges$i, new_edges$j))
  absent <- setdiff(existing_nodes, present)

  lp <- 0.0
  if (length(present) > 0L) lp <- lp + sum(log(p_old[present]))
  if (length(absent) > 0L) lp <- lp + sum(log(1 - p_old[absent]))

  return(list(edge_probs = p_old, logp = as.numeric(lp)))
}
