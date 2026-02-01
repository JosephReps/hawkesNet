# tests/testthat/helper-ba.R

# Helper: build a minimal state object for BA PMF tests.
# Keeps state construction consistent across tests.
make_state_test_helper <- function(deg, born) {
  stopifnot(is.numeric(born) || is.integer(born))
  stopifnot(is.numeric(deg)  || is.integer(deg))
  stopifnot(!is.null(names(deg)), !is.null(names(born)))
  stopifnot(identical(sort(names(deg)), sort(names(born))))

  st <- state_init()
  st$deg  <- deg
  st$born <- born

  # adj isn't used by log_pmf_ba(), but keep it present for consistency
  st$adj <- setNames(vector("list", length(deg)), names(deg))
  for (nm in names(st$adj)) st$adj[[nm]] <- character(0)

  st
}

# Helper: canonical empty edges data.frame (0-row) with required columns.
empty_edges <- function() {
  data.frame(i = character(0), j = character(0))
}

# Helper: enumerate all subsets of BA candidate edges for 1 arrival × old_nodes.
# Returns a list of edges_k data.frames.
enumerate_candidate_edge_subsets <- function(arrival_k, old_nodes) {
  stopifnot(length(arrival_k) == 1L, nzchar(arrival_k))
  arrival_k <- as.character(arrival_k[[1]])
  old_nodes <- as.character(old_nodes)

  if (length(old_nodes) == 0L) {
    return(list(empty_edges()))
  }

  m <- length(old_nodes)  # number of candidate edges (new -> each old)
  out <- vector("list", 2^m)

  for (mask in 0:(2^m - 1)) {
    sel <- which(as.logical(intToBits(mask)[seq_len(m)]))
    if (length(sel) == 0L) {
      out[[mask + 1L]] <- empty_edges()
    } else {
      out[[mask + 1L]] <- data.frame(i = rep(arrival_k, length(sel)),
                                     j = old_nodes[sel])
    }
  }

  out
}
