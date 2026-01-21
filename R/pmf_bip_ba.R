# pmf_ba_bip.R

# Internal: canonicalize undirected edge endpoints and build keys
.edge_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "|")

# Internal: ensure edges_k has i/j, coerce to character, and return canonical keys
.edge_keys_from_df <- function(edges_k) {
  edges_k <- as.data.frame(edges_k, stringsAsFactors = FALSE)
  if (nrow(edges_k) == 0L) return(character(0))

  if (!all(c("i", "j") %in% names(edges_k))) {
    stop("edges_k must be a data.frame with columns i and j.", call. = FALSE)
  }

  i <- as.character(edges_k$i)
  j <- as.character(edges_k$j)

  if (anyNA(i) || anyNA(j) || any(!nzchar(i)) || any(!nzchar(j))) {
    stop("edges_k contains NA/empty endpoint ids.", call. = FALSE)
  }
  if (any(i == j)) stop("Self-loops are not allowed in edges_k.", call. = FALSE)

  key <- .edge_key(i, j)
  if (anyDuplicated(key)) stop("Duplicate undirected edges found in edges_k.", call. = FALSE)
  key
}

# Internal: compute BA-bip attachment probabilities over existing OLD PERP nodes.
# Returns a named numeric vector p_perp with names = old perp node IDs.
edge_probs_ba_bip <- function(state, t_k, beta_edges = 0, delta = 0.001) {
  old_nodes <- names(state$deg)
  if (is.null(old_nodes)) old_nodes <- character(0)
  if (length(old_nodes) == 0L) {
    return(structure(numeric(0), names = character(0)))
  }

  if (!("role" %in% names(state))) stop("state must include `role` for BA-bip.", call. = FALSE)

  old_perps <- old_nodes[
    old_nodes %in% names(state$role) &
      !is.na(state$role[old_nodes]) &
      state$role[old_nodes] == "perp"
  ]

  if (length(old_perps) == 0L) {
    return(structure(numeric(0), names = character(0)))
  }

  age <- t_k - state$born[old_perps]
  if (any(!is.finite(age))) stop("Non-finite ages encountered in born times.", call. = FALSE)

  deg <- as.numeric(state$deg[old_perps])
  if (any(!is.finite(deg))) stop("Non-finite degrees encountered.", call. = FALSE)

  w <- (deg + delta) * exp(-beta_edges * age)
  if (any(!is.finite(w))) stop("Non-finite weights encountered.", call. = FALSE)
  if (any(w < 0)) stop("Negative weights encountered.", call. = FALSE)

  W <- sum(w)
  if (!is.finite(W) || W <= 0) stop("Non-positive total weight encountered.", call. = FALSE)

  p <- w / W
  names(p) <- old_perps

  # keep numeric stability
  p <- pmin(pmax(p, 0), 1)
  p
}

#' Validate BA-bip mark structure at a single event time
#'
#' STRICT rules:
#' - role_k must be provided for arrivals_k and contain exactly one "event" and
#'   all remaining arrivals must be "perp" (no NA roles allowed).
#' - edges_k must contain:
#'   * all deterministic edges (new_event -- each new_perp)
#'   * optionally, Bernoulli edges (new_event -- old_perp)
#' - every edge in edges_k must be incident to the new_event
#' - every non-event endpoint must be a perp (old or new)
#'
#' @return list(new_event=..., new_perps=..., old_perps=..., edge_keys=...)
validate_edges_ba_bip <- function(state, arrivals_k, edges_k, t_k, role_k) {
  stopifnot(is.list(state), all(c("deg", "born", "role") %in% names(state)))
  stopifnot(length(t_k) == 1L, is.numeric(t_k), is.finite(t_k))

  arrivals_k <- as.character(arrivals_k)
  arrivals_k <- arrivals_k[!is.na(arrivals_k) & nzchar(arrivals_k)]
  if (length(arrivals_k) == 0L) stop("BA-bip requires arrivals_k to contain a new event node.", call. = FALSE)

  # role_k: strict, must be provided and aligned
  if (is.null(role_k)) stop("BA-bip requires `role_k` for arrivals_k (strict).", call. = FALSE)

  role_k <- as.character(role_k)
  if (length(role_k) == 1L) {
    role_k <- rep(role_k, length(arrivals_k))
  } else if (length(role_k) != length(arrivals_k)) {
    stop("role_k must be length 1 or length(arrivals_k).", call. = FALSE)
  }

  if (anyNA(role_k)) stop("BA-bip strict: role_k cannot contain NA.", call. = FALSE)
  ok <- role_k %in% c("event", "perp")
  if (!all(ok)) stop("BA-bip strict: role_k values must be 'event' or 'perp'.", call. = FALSE)

  if (sum(role_k == "event") != 1L) {
    stop("BA-bip strict: arrivals_k must contain exactly one 'event' node.", call. = FALSE)
  }

  new_event <- arrivals_k[role_k == "event"]
  new_perps <- arrivals_k[role_k == "perp"]

  # old perps are those in state (role=="perp") excluding any arrivals (new nodes)
  old_nodes <- names(state$deg)
  if (is.null(old_nodes)) old_nodes <- character(0)
  old_perps <- old_nodes[
    old_nodes %in% names(state$role) &
      !is.na(state$role[old_nodes]) &
      state$role[old_nodes] == "perp"
  ]
  old_perps <- setdiff(old_perps, arrivals_k)

  # edges_k validation + keys
  edge_keys <- .edge_keys_from_df(edges_k)

  # every edge must be incident to new_event
  if (length(edge_keys) > 0L) {
    edges_k <- as.data.frame(edges_k, stringsAsFactors = FALSE)
    i <- as.character(edges_k$i)
    j <- as.character(edges_k$j)
    if (!all(i == new_event | j == new_event)) {
      stop("BA-bip strict: every edge in edges_k must be incident to the new event node.", call. = FALSE)
    }

    other <- ifelse(i == new_event, j, i)

    # every other endpoint must be a perp (old or new)
    perps_all <- c(old_perps, new_perps)
    if (!all(other %in% perps_all)) {
      bad <- setdiff(unique(other), perps_all)
      stop("BA-bip strict: edge endpoints other than new_event must be perps. Bad: ",
           paste(bad, collapse = ", "), call. = FALSE)
    }
  }

  # deterministic edges (new_event -- each new_perp) must all be present
  if (length(new_perps) > 0L) {
    det_keys <- .edge_key(rep.int(new_event, length(new_perps)), new_perps)
    missing <- setdiff(det_keys, edge_keys)
    if (length(missing) > 0L) {
      stop("BA-bip strict: missing deterministic edges from new event to new perps.", call. = FALSE)
    }
  }

  list(
    new_event = new_event,
    new_perps = new_perps,
    old_perps = old_perps,
    edge_keys = edge_keys
  )
}

#' Log PMF for the BA-bip mark at one event time
#'
#' log p(mark | history) =
#'   log Pois(K_new | lambda_new)
#' + sum_{old perps i} [ z_i log p_i + (1 - z_i) log(1 - p_i) ]
#'
#' where z_i indicates whether the new event node connects to old perp i at this time.
#' Deterministic edges (new_event -- new_perp) contribute probability 1 and drop out,
#' but are enforced by strict validation.
#'
#' @param state pre-event state
#' @param arrivals_k character vector of arriving node ids
#' @param edges_k data.frame of new edges at this time (i,j)
#' @param t_k event time
#' @param params list containing at least beta_edges and lambda_new
#' @param role_k required (strict): roles for arrivals_k ("event"/"perp"), length 1 or length(arrivals_k)
#' @param delta nonnegative numeric scalar
#' @return numeric scalar log probability (may be -Inf)
#' @export
log_pmf_ba_bip <- function(state, arrivals_k, edges_k, t_k, params, role_k, delta = 0.001) {
  stopifnot(is.list(params))
  stopifnot("beta_edges" %in% names(params))
  stopifnot("lambda_new" %in% names(params))

  beta_edges <- params$beta_edges
  lambda_new <- params$lambda_new

  stopifnot(length(beta_edges) == 1L, is.numeric(beta_edges), is.finite(beta_edges), beta_edges >= 0)
  stopifnot(length(lambda_new) == 1L, is.numeric(lambda_new), is.finite(lambda_new), lambda_new >= 0)
  stopifnot(length(delta) == 1L, is.numeric(delta), is.finite(delta), delta >= 0)

  info <- validate_edges_ba_bip(state, arrivals_k, edges_k, t_k, role_k)

  new_event <- info$new_event
  new_perps <- info$new_perps
  old_perps <- info$old_perps
  edge_keys <- info$edge_keys

  # Poisson term for number of new perps
  K_new <- length(new_perps)
  lp <- stats::dpois(K_new, lambda_new, log = TRUE)

  # Bernoulli-product over OLD perps only
  if (length(old_perps) == 0L) return(lp)

  p <- edge_probs_ba_bip(state, t_k = t_k, beta_edges = beta_edges, delta = delta)

  # Ensure p is aligned to old_perps (should be by construction)
  p <- p[old_perps]
  if (any(is.na(p))) {
    # This can happen if state roles/deg are inconsistent; treat as error
    stop("Internal error: missing probabilities for some old perps.", call. = FALSE)
  }

  # z_i: edge exists between new_event and old_perp i
  z <- .edge_key(rep.int(new_event, length(old_perps)), old_perps) %in% edge_keys

  # Sum log terms (allow -Inf naturally when p==1 but z==FALSE etc.)
  lp <- lp + sum(ifelse(z, log(p), log1p(-p)))
  lp
}
