#' Add node arrivals to the deterministic state
#'
#' Initializes degree/birth/adjacency entries for nodes not yet present.
#' Optionally stores a node "role" (e.g. "perp"/"event") for bipartite kernels.
#' Also updates a cached `network::network` object stored in state$net.
#'
#' @param state list with elements deg, born, adj, role, net, idx
#' @param node_ids character vector of node ids arriving at time `t_k`
#' @param t_k numeric scalar event time (birth time assigned to new nodes)
#' @param role optional character; either length 1 (recycled) or same length as node_ids.
#'   Values should be "perp" or "event" (or NA).
#' @return updated state
#' @export
state_add_nodes <- function(state, node_ids, t_k, role = NULL) {
  stopifnot(
    is.list(state),
    all(c("deg", "born", "adj", "role", "net", "idx") %in% names(state))
  )
  stopifnot(length(t_k) == 1L, is.numeric(t_k), is.finite(t_k))
  if (!inherits(state$net, "network")) stop("state_add_nodes(): state$net must be a 'network' object.", call. = FALSE)

  node_ids_raw <- as.character(node_ids)

  # keep only valid ids, but preserve alignment for role matching
  keep <- !is.na(node_ids_raw) & nzchar(node_ids_raw)
  node_ids <- node_ids_raw[keep]
  if (length(node_ids) == 0L) return(state)

  # normalize role vector for the kept ids
  if (!is.null(role)) {
    role_raw <- as.character(role)
    if (length(role_raw) == 1L) {
      role_kept <- rep(role_raw, sum(keep))
    } else {
      if (length(role_raw) != length(node_ids_raw)) {
        stop("state_add_nodes(): if `role` is not length 1, it must match length(node_ids).", call. = FALSE)
      }
      role_kept <- role_raw[keep]
    }

    ok <- is.na(role_kept) | role_kept %in% c("perp", "event")
    if (!all(ok)) {
      bad <- unique(role_kept[!ok])
      stop(
        "state_add_nodes(): invalid role values: ", paste(bad, collapse = ", "),
        ". Allowed: 'perp', 'event', or NA.",
        call. = FALSE
      )
    }
  } else {
    role_kept <- rep(NA_character_, length(node_ids))
  }

  # existing nodes
  existing <- names(state$deg)
  if (is.null(existing)) existing <- character(0)

  # only add genuinely new nodes
  node_ids_u <- unique(node_ids)
  new_nodes  <- setdiff(node_ids_u, existing)
  if (length(new_nodes) == 0L) return(state)

  # map roles to unique ids: take first occurrence
  role_map <- stats::setNames(role_kept, node_ids)
  role_new <- unname(role_map[new_nodes])
  if (length(role_new) != length(new_nodes)) role_new <- rep(NA_character_, length(new_nodes))

  # --- Update deterministic vectors/lists ---
  state$deg[new_nodes]  <- 0L
  state$born[new_nodes] <- t_k

  if (is.null(names(state$adj))) names(state$adj) <- character(0)
  state$adj[new_nodes] <- replicate(length(new_nodes), character(0), simplify = FALSE)

  if (is.null(names(state$role))) names(state$role) <- character(0)
  state$role[new_nodes] <- role_new

  # --- Update cached network and idx map ---
  nv_old <- network::network.size(state$net)
  network::add.vertices(state$net, nv = length(new_nodes))
  nv_new <- network::network.size(state$net)

  new_vids <- seq.int(nv_old + 1L, nv_new)

  # idx map: id -> vertex index
  state$idx[new_nodes] <- as.integer(new_vids)

  # Set vertex attributes in net to mirror state
  network::set.vertex.attribute(state$net, "name", value = as.character(new_nodes), v = new_vids)
  network::set.vertex.attribute(state$net, "born", value = rep(t_k, length(new_nodes)), v = new_vids)
  network::set.vertex.attribute(state$net, "role", value = role_new, v = new_vids)

  state
}
