# R/sim_cs.R
#
# Simulator for CS mark kernel (arrivals + edges) for one event.
# Requires ernm + network and uses cached state$net + state$idx.

# Internal: generate the next m node ids from the given state
next_node_ids <- function(state, m) {
  stopifnot(length(m) == 1L, is.numeric(m), is.finite(m), m >= 0)
  m <- as.integer(m)
  if (m == 0L) return(character(0))

  existing <- names(state$deg)
  if (is.null(existing) || length(existing) == 0L) {
    return(as.character(seq_len(m)))
  }

  suppressWarnings(ints <- as.integer(existing))
  ints <- ints[is.finite(ints)]

  if (length(ints) == 0L) {
    # fallback for non-integer IDs
    start <- length(existing) + 1L
    return(paste0("n", seq.int(start, start + m - 1L)))
  }

  start <- max(ints) + 1L
  as.character(seq.int(start, start + m - 1L))
}

#' Simulate CS mark (arrivals + edges) for one event
#'
#' For one accepted event time t_k:
#' 1) Sample number of arrivals M ~ Pois(node_lambda) if t_k <= max_node_time, else M=0.
#' 2) Add M new nodes.
#' 3) For each new node u, consider candidate edges to old nodes v (deterministic truncation).
#' 4) Sample each candidate edge independently with probability p_uv from edge_probs_cs().
#'
#' @param state Pre-event state (must include deg, born, net, idx).
#' @param t_k Numeric scalar event time.
#' @param params List with beta_edges, node_lambda, CS_params.
#' @param formula_rhs RHS string for ERNM change stats, e.g. "edges + triangles".
#' @param truncation integer >= 1 or Inf; deterministic candidate set size per new node.
#' @param max_node_time numeric scalar; if t_k > max_node_time, set arrivals=0.
#' @param model_cache optional cache (env/list) holding ERNM model object.
#' @param new_nodes Optional character vector to force arriving node ids (useful for tests).
#'   If provided, overrides the Poisson draw and uses length(new_nodes) as M.
#'
#' @return list(arrivals = <chr vec>, edges = data.frame(i,j))
#' @export
sim_mark_cs <- function(state,
                        t_k,
                        params,
                        formula_rhs,
                        truncation = Inf,
                        max_node_time = Inf,
                        model_cache = NULL,
                        new_nodes = NULL,
                        eps = 1e-12) {
  # hard requirements
  if (!requireNamespace("network", quietly = TRUE)) {
    stop("CS kernel requires package 'network'. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("ernm", quietly = TRUE)) {
    stop("CS kernel requires package 'ernm'. Please install it.", call. = FALSE)
  }

  stopifnot(is.list(state), all(c("deg", "born", "net", "idx") %in% names(state)))
  stopifnot(inherits(state$net, "network"))
  stopifnot(length(t_k) == 1L, is.numeric(t_k), is.finite(t_k))
  stopifnot(is.list(params))

  for (nm in c("beta_edges", "node_lambda", "CS_params")) {
    if (!(nm %in% names(params))) stop("sim_mark_cs(): missing params$", nm, call. = FALSE)
  }
  beta_edges <- params$beta_edges
  node_lambda <- params$node_lambda
  CS_params <- params$CS_params

  stopifnot(is.numeric(beta_edges), length(beta_edges) == 1L, is.finite(beta_edges), beta_edges >= 0)
  stopifnot(is.numeric(node_lambda), length(node_lambda) == 1L, is.finite(node_lambda), node_lambda > 0)
  stopifnot(is.numeric(CS_params), all(is.finite(CS_params)))
  stopifnot(is.character(formula_rhs), length(formula_rhs) == 1L, nzchar(formula_rhs))
  stopifnot(is.numeric(max_node_time), length(max_node_time) == 1L, is.finite(max_node_time))

  # --- choose arrivals ---
  if (!is.null(new_nodes)) {
    new_nodes <- unique(as.character(new_nodes))
    new_nodes <- new_nodes[!is.na(new_nodes) & nzchar(new_nodes)]
    M <- length(new_nodes)
  } else {
    if (t_k > max_node_time) {
      M <- 0L
      new_nodes <- character(0)
    } else {
      M <- stats::rpois(1L, lambda = node_lambda)
      if (M == 0L) {
        new_nodes <- character(0)
      } else {
        new_nodes <- next_node_ids(state, M)
      }
    }
  }

  # No arrivals => no edges
  if (length(new_nodes) == 0L) {
    return(list(
      arrivals = character(0),
      edges = data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
    ))
  }

  # Pre-event old nodes (before adding)
  old_nodes <- names(state$deg)
  if (is.null(old_nodes)) old_nodes <- character(0)

  # Add new nodes to the state so idx/net include them for change stats
  # (this updates the real state in the caller typically via state_step; here we only return marks,
  #  so we do NOT mutate the passed-in state; we compute probs using an augmented copy)
  # We'll build a temporary augmented state just for probability calculations.
  net2 <- if (exists("network.copy", where = asNamespace("network"), inherits = FALSE)) {
    network::network.copy(state$net)
  } else {
    unserialize(serialize(state$net, NULL))
  }
  idx2 <- state$idx
  born2 <- state$born

  nv_old <- network::network.size(net2)
  network::add.vertices(net2, nv = length(new_nodes))
  nv_new <- network::network.size(net2)
  new_vids <- seq.int(nv_old + 1L, nv_new)

  network::set.vertex.attribute(net2, "name", value = new_nodes, v = new_vids)
  network::set.vertex.attribute(net2, "born", value = rep(t_k, length(new_nodes)), v = new_vids)

  idx2[new_nodes] <- as.integer(new_vids)
  born2[new_nodes] <- t_k

  state2 <- state
  state2$net <- net2
  state2$idx <- idx2
  state2$born <- born2

  # If no old nodes, then no edges possible
  if (length(old_nodes) == 0L) {
    return(list(
      arrivals = new_nodes,
      edges = data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
    ))
  }

  # Candidate probabilities
  p <- edge_probs_cs(
    state = state2,
    new_nodes = new_nodes,
    old_nodes = old_nodes,
    t_k = t_k,
    CS_params = CS_params,
    beta_edges = beta_edges,
    formula_rhs = formula_rhs,
    truncation = truncation,
    model_cache = model_cache,
    eps = eps
  )

  if (length(p) == 0L) {
    return(list(
      arrivals = new_nodes,
      edges = data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
    ))
  }

  add <- stats::runif(length(p)) < p
  if (!any(add)) {
    return(list(
      arrivals = new_nodes,
      edges = data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
    ))
  }

  keys <- names(p)[add]
  # split "u|v" back into columns
  parts <- strsplit(keys, "\\|", fixed = FALSE)
  ii <- vapply(parts, `[[`, character(1), 1L)
  jj <- vapply(parts, `[[`, character(1), 2L)

  edges <- data.frame(i = ii, j = jj, stringsAsFactors = FALSE)

  list(arrivals = new_nodes, edges = edges)
}
