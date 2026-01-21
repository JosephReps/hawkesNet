# R/pmf_cs.R
#
# CS (change-statistics) Bernoulli-product mark kernel.
# Requires:
#   state$net : network::network (undirected)
#   state$idx : named integer, mapping node-id -> vertex index in state$net
#   state$born: named numeric, node entrance times
#
# Key fix vs previous version:
#   DO NOT construct "UndirectedNet". Use createCppModel(net ~ rhs) and
#   model$setNetwork(ernm::as.BinaryNet(net)) like the original package.

.cs_require_pkgs <- function() {
  if (!requireNamespace("network", quietly = TRUE)) {
    stop("CS kernel requires package 'network'. Please install it.", call. = FALSE)
  }
  if (!requireNamespace("ernm", quietly = TRUE)) {
    stop("CS kernel requires package 'ernm'. Please install it.", call. = FALSE)
  }
}

.cs_net_copy <- function(net) {
  if (exists("network.copy", where = asNamespace("network"), inherits = FALSE)) {
    return(network::network.copy(net))
  }
  unserialize(serialize(net, NULL))
}

# deterministic candidate set: each new node considers the first `truncation`
# old nodes in sorted order (no randomness in the candidate set).
.cs_candidates_new_old <- function(new_nodes, old_nodes, truncation = Inf) {
  stopifnot(is.character(new_nodes), is.character(old_nodes))
  if (length(new_nodes) == 0L || length(old_nodes) == 0L) {
    return(list(tail = character(0), head = character(0)))
  }

  old_sorted <- sort(old_nodes)

  if (is.finite(truncation)) {
    truncation <- as.integer(truncation)
    if (!(length(truncation) == 1L && is.finite(truncation) && truncation >= 1L)) {
      stop("CS truncation must be an integer >= 1 (or Inf).", call. = FALSE)
    }
  }

  tail_out <- character(0)
  head_out <- character(0)

  for (u in new_nodes) {
    heads_u <- if (is.finite(truncation)) old_sorted[seq_len(min(truncation, length(old_sorted)))] else old_sorted
    tail_out <- c(tail_out, rep.int(u, length(heads_u)))
    head_out <- c(head_out, heads_u)
  }

  list(tail = tail_out, head = head_out)
}

# Validate observed edges: must be new<->old only, no duplicates, no loops.
validate_edges_cs <- function(edges_k, new_nodes, old_nodes) {
  stopifnot(is.data.frame(edges_k), all(c("i", "j") %in% names(edges_k)))

  new_nodes <- as.character(new_nodes)
  old_nodes <- as.character(old_nodes)

  if (nrow(edges_k) == 0L) return(list(new_end = character(0), old_end = character(0)))

  i <- as.character(edges_k$i)
  j <- as.character(edges_k$j)

  if (anyNA(i) || anyNA(j) || any(!nzchar(i)) || any(!nzchar(j))) {
    stop("CS edges_k contains NA/empty node IDs.", call. = FALSE)
  }
  if (any(i == j)) stop("CS kernel does not allow self-loops.", call. = FALSE)

  key <- paste(pmin(i, j), pmax(i, j), sep = "|")
  if (anyDuplicated(key)) stop("CS kernel observed duplicate edges in one event.", call. = FALSE)

  is_i_new <- i %in% new_nodes
  is_j_new <- j %in% new_nodes
  is_i_old <- i %in% old_nodes
  is_j_old <- j %in% old_nodes

  ok <- (is_i_new & is_j_old) | (is_j_new & is_i_old)
  if (!all(ok)) {
    bad <- which(!ok)[1]
    stop(
      "CS kernel requires each observed edge to be between a new node and an old node.\n",
      "Bad edge: (", i[[bad]], ", ", j[[bad]], ").",
      call. = FALSE
    )
  }

  new_end <- ifelse(is_i_new, i, j)
  old_end <- ifelse(is_i_new, j, i)

  list(new_end = as.character(new_end), old_end = as.character(old_end))
}

# Old-code-compatible: ensure at least 4 vertices + remove 'na' attribute before C++
.cs_prepare_net_for_cpp <- function(net) {
  .cs_require_pkgs()
  stopifnot(inherits(net, "network"))

  net2 <- net

  nv <- network::network.size(net2)
  if (nv < 4L) {
    network::add.vertices(net2, nv = 4L - nv)
  }

  # match old code: delete "na" attribute if present
  if ("na" %in% network::list.vertex.attributes(net2)) {
    network::delete.vertex.attribute(net2, "na")
  }

  net2
}

# Get/create (and cache) the ERNM cpp model, but ALWAYS setNetwork(as.BinaryNet(net))
.cs_get_model <- function(net, formula_rhs, model_cache = NULL) {
  .cs_require_pkgs()
  stopifnot(inherits(net, "network"))
  stopifnot(is.character(formula_rhs), length(formula_rhs) == 1L, nzchar(formula_rhs))

  # cache hit?
  if (!is.null(model_cache) && !is.null(model_cache$model) &&
      is.character(model_cache$rhs) && identical(model_cache$rhs, formula_rhs)) {
    mdl <- model_cache$model
  } else {
    # bind a symbol into the formula env like the old package expects
    f_env <- new.env(parent = emptyenv())
    f_env$new_net <- net
    f <- stats::as.formula(paste0("new_net ~ ", formula_rhs), env = f_env)
    mdl <- ernm::createCppModel(f)

    if (!is.null(model_cache)) {
      model_cache$rhs <- formula_rhs
      model_cache$model <- mdl
    }
  }

  # ALWAYS update network before use (like old code does via setNetwork)
  net_cpp <- .cs_prepare_net_for_cpp(net)
  mdl$setNetwork(ernm::as.BinaryNet(net_cpp))
  mdl$calculate()

  mdl
}

# Compute Bernoulli probs for candidate edges new->old
edge_probs_cs <- function(state,
                          new_nodes,
                          old_nodes,
                          t_k,
                          CS_params,
                          beta_edges = 0,
                          formula_rhs,
                          truncation = Inf,
                          model_cache = NULL,
                          eps = 1e-12) {
  .cs_require_pkgs()

  stopifnot(is.list(state), all(c("net", "idx", "born") %in% names(state)))
  stopifnot(inherits(state$net, "network"))
  stopifnot(is.numeric(t_k), length(t_k) == 1L, is.finite(t_k))
  stopifnot(is.numeric(beta_edges), length(beta_edges) == 1L, is.finite(beta_edges), beta_edges >= 0)
  stopifnot(is.numeric(CS_params), all(is.finite(CS_params)))
  stopifnot(is.character(formula_rhs), length(formula_rhs) == 1L, nzchar(formula_rhs))
  stopifnot(is.numeric(eps), length(eps) == 1L, is.finite(eps), eps > 0, eps < 0.5)

  new_nodes <- unique(as.character(new_nodes))
  old_nodes <- unique(as.character(old_nodes))

  if (length(new_nodes) == 0L || length(old_nodes) == 0L) {
    return(structure(numeric(0), names = character(0)))
  }

  cand <- .cs_candidates_new_old(new_nodes, old_nodes, truncation = truncation)
  tails <- cand$tail
  heads <- cand$head
  if (length(tails) == 0L) return(structure(numeric(0), names = character(0)))

  tail_idx <- state$idx[tails]
  head_idx <- state$idx[heads]
  if (any(is.na(tail_idx)) || any(is.na(head_idx))) {
    stop("edge_probs_cs(): missing idx mapping for some candidate nodes. Did you add arrivals to state first?",
         call. = FALSE)
  }

  mdl <- .cs_get_model(state$net, formula_rhs = formula_rhs, model_cache = model_cache)

  C <- mdl$computeChangeStats(as.integer(tail_idx), as.integer(head_idx))
  if (!is.matrix(C)) C <- as.matrix(C)

  if (ncol(C) != length(CS_params)) {
    stop(
      "edge_probs_cs(): length(CS_params) does not match number of change-stat columns.\n",
      "ncol(change_stats) = ", ncol(C), ", length(CS_params) = ", length(CS_params),
      call. = FALSE
    )
  }

  eta <- as.vector(C %*% CS_params)
  if (any(!is.finite(eta))) stop("edge_probs_cs(): non-finite linear predictors.", call. = FALSE)

  p0 <- stats::plogis(eta)

  # decay uses OLD node entrance times (node arrival times)
  age <- t_k - state$born[heads]
  if (any(!is.finite(age))) stop("edge_probs_cs(): non-finite ages from born times.", call. = FALSE)
  if (any(age < 0)) stop("edge_probs_cs(): negative ages encountered (born time after t_k).", call. = FALSE)

  p <- p0 * exp(-beta_edges * age)
  p <- pmin(pmax(p, eps), 1 - eps)

  names(p) <- paste(tails, heads, sep = "|")
  p
}

#' Log PMF for one CS event (Poisson arrivals + Bernoulli product over candidates)
#'
#' @export
log_pmf_cs <- function(state,
                       arrivals_k,
                       edges_k,
                       t_k,
                       params,
                       formula_rhs,
                       truncation = Inf,
                       max_node_time = Inf,
                       model_cache = NULL,
                       eps = 1e-12) {
  .cs_require_pkgs()

  stopifnot(is.list(state), all(c("net", "idx", "born") %in% names(state)))
  stopifnot(inherits(state$net, "network"))
  stopifnot(is.data.frame(edges_k), all(c("i", "j") %in% names(edges_k)))
  stopifnot(is.numeric(t_k), length(t_k) == 1L, is.finite(t_k))
  stopifnot(is.list(params))

  for (nm in c("beta_edges", "node_lambda", "CS_params")) {
    if (!(nm %in% names(params))) stop("log_pmf_cs(): missing params$", nm, call. = FALSE)
  }

  beta_edges  <- params$beta_edges
  node_lambda <- params$node_lambda
  CS_params   <- params$CS_params

  stopifnot(is.numeric(beta_edges), length(beta_edges) == 1L, is.finite(beta_edges), beta_edges >= 0)
  stopifnot(is.numeric(node_lambda), length(node_lambda) == 1L, is.finite(node_lambda), node_lambda > 0)
  stopifnot(is.numeric(CS_params), all(is.finite(CS_params)))
  stopifnot(is.character(formula_rhs), length(formula_rhs) == 1L, nzchar(formula_rhs))
  stopifnot(is.numeric(max_node_time), length(max_node_time) == 1L, is.finite(max_node_time))

  arrivals_k <- unique(as.character(arrivals_k))
  arrivals_k <- arrivals_k[!is.na(arrivals_k) & nzchar(arrivals_k)]

  # old nodes = current state node IDs (use idx names)
  old_nodes <- names(state$idx)
  if (is.null(old_nodes)) old_nodes <- character(0)
  old_nodes <- old_nodes[!is.na(old_nodes) & nzchar(old_nodes)]
  old_nodes <- unique(old_nodes)

  if (length(intersect(arrivals_k, old_nodes)) > 0L) {
    bad <- intersect(arrivals_k, old_nodes)
    stop("log_pmf_cs(): arrivals_k contains nodes already present in state: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }

  # Poisson arrivals term (like old code’s dpois(new_nodes-old_nodes, node_lambda)) :contentReference[oaicite:6]{index=6}
  if (t_k > max_node_time) {
    if (length(arrivals_k) > 0L) stop("log_pmf_cs(): t_k > max_node_time but arrivals_k is non-empty.", call. = FALSE)
    ll_nodes <- 0.0
  } else {
    ll_nodes <- stats::dpois(length(arrivals_k), lambda = node_lambda, log = TRUE)
  }

  # no old nodes => no edges possible
  if (length(old_nodes) == 0L) {
    if (nrow(edges_k) != 0L) stop("log_pmf_cs(): no old nodes exist, but edges were observed.", call. = FALSE)
    return(as.numeric(ll_nodes))
  }

  # no arrivals => no candidate edges; edges must be empty
  if (length(arrivals_k) == 0L) {
    if (nrow(edges_k) != 0L) stop("log_pmf_cs(): no arrivals but edges were observed.", call. = FALSE)
    return(as.numeric(ll_nodes))
  }

  # augment state with the new isolated vertices *before* computing change-stats
  net2  <- .cs_net_copy(state$net)
  idx2  <- state$idx
  born2 <- state$born

  nv_old <- network::network.size(net2)
  network::add.vertices(net2, nv = length(arrivals_k))
  new_vids <- seq.int(nv_old + 1L, nv_old + length(arrivals_k))

  # set vertex "name" if you want; mapping is what matters
  network::set.vertex.attribute(net2, "name", value = arrivals_k, v = new_vids)
  idx2[arrivals_k]  <- as.integer(new_vids)
  born2[arrivals_k] <- t_k

  state2 <- state
  state2$net  <- net2
  state2$idx  <- idx2
  state2$born <- born2

  # validate observed edges are new<->old only
  ends <- validate_edges_cs(edges_k, new_nodes = arrivals_k, old_nodes = old_nodes)

  p <- edge_probs_cs(
    state       = state2,
    new_nodes   = arrivals_k,
    old_nodes   = old_nodes,
    t_k         = t_k,
    CS_params   = CS_params,
    beta_edges  = beta_edges,
    formula_rhs = formula_rhs,
    truncation  = truncation,
    model_cache = model_cache,
    eps         = eps
  )

  if (length(p) == 0L) {
    if (nrow(edges_k) == 0L) return(as.numeric(ll_nodes))
    stop("log_pmf_cs(): no candidate edges but edges were observed.", call. = FALSE)
  }

  present_keys <- paste(ends$new_end, ends$old_end, sep = "|")

  missing <- setdiff(present_keys, names(p))
  if (length(missing) > 0L) {
    stop(
      "log_pmf_cs(): observed edges include pairs not in the candidate set (possibly due to truncation).\n",
      "Example missing candidate: ", missing[[1]],
      call. = FALSE
    )
  }

  in_mark <- names(p) %in% present_keys
  ll_edges <- sum(if (any(in_mark)) log(p[in_mark]) else 0) +
    sum(if (any(!in_mark)) log1p(-p[!in_mark]) else 0)

  as.numeric(ll_nodes + ll_edges)
}
