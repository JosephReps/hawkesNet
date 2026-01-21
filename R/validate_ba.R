# R/validate_ba.R

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

