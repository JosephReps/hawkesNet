#' Update deterministic state after an event's edge batch
#'
#' Also updates the cached `network::network` object stored in state$net.
#'
#' @param state list with elements deg, born, adj, role, net, idx
#' @param edges_k data.frame with columns i, j (may be 0-row)
#' @param t_k numeric scalar event time
#' @param update_adj logical; if TRUE, update adjacency lists
#' @return updated state
#' @export
state_add_edges <- function(state, edges_k, t_k, update_adj = TRUE) {
  stopifnot(
    is.list(state),
    all(c("deg", "born", "adj", "role", "net", "idx") %in% names(state))
  )
  stopifnot(is.data.frame(edges_k), all(c("i", "j") %in% names(edges_k)))
  stopifnot(length(t_k) == 1L, is.numeric(t_k), is.finite(t_k))
  if (!inherits(state$net, "network")) stop("state_add_edges(): state$net must be a 'network' object.", call. = FALSE)

  if (nrow(edges_k) == 0L) return(state)

  # Ensure node ids are character keys for naming
  i_chr <- as.character(edges_k$i)
  j_chr <- as.character(edges_k$j)

  if (anyNA(i_chr) || anyNA(j_chr) || any(!nzchar(i_chr)) || any(!nzchar(j_chr))) {
    stop("state_add_edges(): edges_k contains NA/empty node IDs.", call. = FALSE)
  }

  # No self loops
  if (any(i_chr == j_chr)) {
    stop("state_add_edges(): self-loops are not allowed.", call. = FALSE)
  }

  # No duplicate edges within event (undirected)
  key_evt <- paste(pmin(i_chr, j_chr), pmax(i_chr, j_chr), sep = "|")
  if (anyDuplicated(key_evt)) {
    stop("state_add_edges(): duplicate edges within event.", call. = FALSE)
  }

  # Make sure all edge endpoints already exist
  nodes_evt <- unique(c(i_chr, j_chr))
  existing <- names(state$deg)
  if (is.null(existing)) existing <- character(0)

  missing <- setdiff(nodes_evt, existing)
  if (length(missing) > 0L) {
    stop(
      "state_add_edges(): edge endpoints not in state. Call state_add_nodes() first. Missing nodes: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  # Also ensure idx mapping exists for all endpoints
  if (any(is.na(state$idx[nodes_evt])) || any(state$idx[nodes_evt] <= 0L)) {
    bad <- nodes_evt[is.na(state$idx[nodes_evt]) | state$idx[nodes_evt] <= 0L]
    stop("state_add_edges(): missing/invalid state$idx entries for: ", paste(bad, collapse = ", "), call. = FALSE)
  }

  # No duplicate edges against existing state (simple graph assumption)
  # Use adjacency lists (fast) if update_adj=TRUE or adj exists.
  # If adj is empty/missing, fall back to network edge check.
  for (r in seq_len(nrow(edges_k))) {
    a <- i_chr[r]
    b <- j_chr[r]
    if (!is.null(state$adj[[a]]) && b %in% state$adj[[a]]) {
      stop("state_add_edges(): edge already exists in state: ", a, " - ", b, call. = FALSE)
    }
    if (!is.null(state$adj[[b]]) && a %in% state$adj[[b]]) {
      stop("state_add_edges(): edge already exists in state: ", a, " - ", b, call. = FALSE)
    }
  }

  # Update degrees
  inc <- table(c(i_chr, j_chr))
  state$deg[names(inc)] <- state$deg[names(inc)] + as.integer(inc)

  # Update adjacency lists
  if (update_adj) {
    m <- nrow(edges_k)
    for (r in seq_len(m)) {
      a <- i_chr[r]
      b <- j_chr[r]
      state$adj[[a]] <- unique(c(as.character(state$adj[[a]]), b))
      state$adj[[b]] <- unique(c(as.character(state$adj[[b]]), a))
    }
  }

  # Update cached network object (add edges by vertex index)
  tail_idx <- as.integer(state$idx[i_chr])
  head_idx <- as.integer(state$idx[j_chr])

  # network::add.edges expects tails/heads vectors, adds each pair
  network::add.edges(state$net, tail = tail_idx, head = head_idx)

  state
}
