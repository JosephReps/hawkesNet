# R/sim_cs_bip.R

# Internal: generate the next m node ids from the given net (same helper as sim_cs.R)
.next_node_ids_net <- function(net, m) {
  stopifnot(length(m) == 1L, is.numeric(m), is.finite(m), m >= 0)
  m <- as.integer(m)
  if (m == 0L) return(character(0))
  ids <- network::network.vertex.names(net)
  if (is.null(ids) || length(ids) == 0L) {
    start <- 1L
  } else {
    suppressWarnings({
      nums <- as.integer(ids)
    })
    if (all(!is.na(nums))) start <- max(nums) + 1L else start <- length(ids) + 1L
  }
  as.character(seq.int(start, length.out = m))
}

#' Simulate CS-bip mark (arrivals + edges) for one event
#'
#' Rules (BA-bip style):
#' - exactly 1 new "event" node
#' - K_new new "perp" nodes, K_new ~ Poisson(lambda_new)
#' - deterministic edges: event -- each new perp
#' - optional edges: event -- each OLD perp with CS probability
#'
#' @param net pre-event network::network
#' @param t_k event time
#' @param params list with beta_edges, lambda_new, and CS_* parameters
#' @param formula_rhs ERNM RHS string, e.g. "edges + gwesp(0.5, fixed=TRUE)"
#' @param model_cache optional cache list(rhs=..., model=...) (passed to edge_probs_cs)
#' @param eps clamp
#' @return list(arrivals=data.frame(id,role), edges=data.frame(i,j), model_cache=...)
#' @export
sim_mark_cs_bip <- function(
    net, t_k, params,
    formula_rhs,
    model_cache = NULL,
    eps = 1e-12
) {
  stopifnot(is.numeric(t_k), length(t_k) == 1L, is.finite(t_k))
  stopifnot(is.list(params))
  if (!exists("edge_probs_cs", mode = "function")) {
    stop("sim_mark_cs_bip(): edge_probs_cs() not found. Source/load pmf_cs.R.", call. = FALSE)
  }
  if (!exists(".cs_augment_net_with_new_nodes", mode = "function")) {
    stop("sim_mark_cs_bip(): .cs_augment_net_with_new_nodes() not found. Source/load pmf_cs.R.", call. = FALSE)
  }

  beta_edges <- params[["beta_edges"]]
  lambda_new <- params[["lambda_new"]]

  stopifnot(is.numeric(beta_edges), length(beta_edges) == 1L, is.finite(beta_edges), beta_edges >= 0)
  stopifnot(is.numeric(lambda_new), length(lambda_new) == 1L, is.finite(lambda_new), lambda_new >= 0)

  # Sample arrivals
  K_new <- as.integer(stats::rpois(1L, lambda = lambda_new))

  # New IDs: 1 event + K_new perps (must allocate in ONE call to avoid duplicates)
  ids_all  <- .next_node_ids_net(net, 1L + K_new)
  new_event <- ids_all[1]
  new_perps <- if (K_new > 0L) ids_all[-1] else character(0)

  arrivals <- data.frame(
    id   = c(new_event, new_perps),
    role = c("event", rep("perp", length(new_perps))),
    stringsAsFactors = FALSE
  )

  # Deterministic edges: event -- new perps
  edges_det <- if (length(new_perps) == 0L) {
    data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
  } else {
    data.frame(i = rep(new_event, length(new_perps)), j = new_perps, stringsAsFactors = FALSE)
  }

  # OLD perps (role == "perp")
  roles <- network::get.vertex.attribute(net, "role")
  ids <- network::network.vertex.names(net)
  old_perps <- character(0)
  if (!is.null(ids) && length(ids) > 0L && !is.null(roles)) {
    roles <- as.character(roles)
    old_perps <- as.character(ids[which(!is.na(roles) & roles == "perp")])
  }

  edges_opt <- data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)

  if (length(old_perps) > 0L) {
    # Augment net with new nodes for change-stats calculations
    net2 <- .cs_augment_net_with_new_nodes(net, c(new_event, new_perps), t_k = t_k)

    # Candidate probs only for: new_event -> each old perp
    res <- edge_probs_cs(
      net = net2,
      new_nodes = new_event,
      old_nodes = old_perps,
      t_k = t_k,
      params = params,
      formula_rhs = formula_rhs,
      truncation = Inf,
      model_cache = model_cache,
      beta_edges = beta_edges,
      eps = eps
    )

    if (is.list(res) && !is.null(res$p)) {
      p <- res$p
      model_cache <- res$model_cache
    } else {
      p <- res
    }

    if (length(p) != length(old_perps)) {
      stop("sim_mark_cs_bip(): edge_probs_cs returned wrong length.", call. = FALSE)
    }

    keep <- stats::runif(length(p)) < p
    if (any(keep)) {
      edges_opt <- data.frame(i = rep(new_event, sum(keep)), j = old_perps[keep], stringsAsFactors = FALSE)
    }
  }

  edges <- rbind(edges_det, edges_opt)

  list(arrivals = arrivals, edges = edges, model_cache = model_cache)
}
