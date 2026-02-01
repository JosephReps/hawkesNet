# R/sim_cs.R
#
# Simulator for CS mark kernel (arrivals + edges) for one event (NETWORK-BASED STATE).
#
# Parameter convention (no nested vectors):
#   - beta_edges   : numeric scalar >= 0
#   - node_lambda  : numeric scalar > 0
#   - CS_*         : numeric scalars, one per ERNM change-stat column
#                   (see ?log_pmf_cs / ?edge_probs_cs)
#
# Requires packages: network, ernm, sna
#
# NOTE:
#   This simulator returns:
#     list(arrivals = data.frame(id=...), edges = data.frame(i=..., j=...))
#   matching the BA / BA-bip kernels' "net_add_event()" input structure.

# Internal: generate the next m node ids from the given net
.next_node_ids_net <- function(net, m) {
  stopifnot(length(m) == 1L, is.numeric(m), is.finite(m), m >= 0)
  m <- as.integer(m)
  if (m == 0L) return(character(0))

  existing <- network::network.vertex.names(net)
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
#' 2) Add M new nodes (conceptually; we augment a copy for change-stats).
#' 3) For each new node u, consider candidate edges to old nodes v (deterministic truncation).
#' 4) Sample each candidate edge independently with probability p_uv from edge_probs_cs().
#'
#' @param net Pre-event network state (network::network).
#' @param t_k Numeric scalar event time.
#' @param params List with beta_edges, node_lambda, and CS_* parameters.
#' @param formula_rhs RHS string for ERNM change stats, e.g. "edges + triangles".
#' @param truncation integer >= 1 or Inf; deterministic candidate set size per new node.
#' @param max_node_time numeric scalar; if t_k > max_node_time, set arrivals=0.
#' @param model_cache optional cache (list) holding ERNM model; passed through to edge_probs_cs().
#' @param new_nodes optional character vector to force the arriving node ids (tests).
#' @param eps numeric in (0,0.5) for probability clamping.
#'
#' @return list(arrivals = data.frame(id=...), edges = data.frame(i,j))
#' @export
sim_mark_cs <- function(net,
                        t_k,
                        params,
                        formula_rhs,
                        truncation = Inf,
                        max_node_time = 1e10,
                        model_cache = NULL,
                        new_nodes = NULL,
                        eps = 1e-12) {
  stopifnot(inherits(net, "network"))
  stopifnot(length(t_k) == 1L, is.numeric(t_k), is.finite(t_k))
  stopifnot(is.list(params))
  stopifnot(is.character(formula_rhs), length(formula_rhs) == 1L, nzchar(formula_rhs))
  stopifnot(is.numeric(max_node_time), length(max_node_time) == 1L, is.finite(max_node_time))

  for (nm in c("beta_edges", "node_lambda")) {
    if (!(nm %in% names(params))) stop("sim_mark_cs(): missing params$", nm, call. = FALSE)
  }
  beta_edges <- params$beta_edges
  node_lambda <- params$node_lambda

  stopifnot(is.numeric(beta_edges), length(beta_edges) == 1L, is.finite(beta_edges), beta_edges >= 0)
  stopifnot(is.numeric(node_lambda), length(node_lambda) == 1L, is.finite(node_lambda), node_lambda > 0)

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
        new_nodes <- .next_node_ids_net(net, M)
      }
    }
  }

  # Return early if no arrivals
  if (length(new_nodes) == 0L) {
    return(list(
      arrivals = data.frame(id = character(0), stringsAsFactors = FALSE),
      edges = data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE),
      model_cache = model_cache
    ))
  }

  # Pre-event old nodes
  old_nodes <- network::network.vertex.names(net)
  if (is.null(old_nodes)) old_nodes <- character(0)

  # If no old nodes, then no edges possible
  if (length(old_nodes) == 0L) {
    return(list(
      arrivals = data.frame(id = new_nodes, stringsAsFactors = FALSE),
      edges = data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE),
      model_cache = model_cache
    ))
  }

  # Augment a copy of net so change-stats include new nodes
  if (!exists(".cs_augment_net_with_new_nodes", mode = "function")) {
    stop("sim_mark_cs(): .cs_augment_net_with_new_nodes() not found. Did you source pmf_cs.R first?", call. = FALSE)
  }
  net2 <- .cs_augment_net_with_new_nodes(net, new_nodes, t_k = t_k)

  # Candidate probabilities (CS_* params are read inside edge_probs_cs via colnames(change_stats))
  res <- edge_probs_cs(
    net = net2,
    new_nodes = new_nodes,
    old_nodes = old_nodes,
    t_k = t_k,
    params = params,
    formula_rhs = formula_rhs,
    truncation = truncation,
    model_cache = model_cache,
    beta_edges = beta_edges,
    eps = eps
  )

  # Backwards/forwards compatibility: edge_probs_cs() may return either a numeric vector
  # or list(p=..., model_cache=...)
  if (is.list(res) && !is.null(res$p)) {
    p <- res$p
    model_cache <- res$model_cache
  } else {
    p <- res
  }

  if (length(p) == 0L) {
    return(list(
      arrivals = data.frame(id = new_nodes, stringsAsFactors = FALSE),
      edges = data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE),
      model_cache = model_cache
    ))
  }

  add <- stats::runif(length(p)) < p
  if (!any(add)) {
    return(list(
      arrivals = data.frame(id = new_nodes, stringsAsFactors = FALSE),
      edges = data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE),
      model_cache = model_cache
    ))
  }

  keys <- names(p)[add]
  parts <- strsplit(keys, "\\|", fixed = FALSE)
  ii <- vapply(parts, `[[`, character(1), 1L)
  jj <- vapply(parts, `[[`, character(1), 2L)

  edges <- data.frame(i = ii, j = jj, stringsAsFactors = FALSE)

  list(
    arrivals = data.frame(id = new_nodes, stringsAsFactors = FALSE),
    edges = edges
  )
}
