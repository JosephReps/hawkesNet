#' Construct an initial network state (baseline graph) at observation start T0
#'
#' @export
make_state <- function(
    nodes0 = NULL,
    edges0 = NULL,
    T0 = 0,
    implicit_nodes = TRUE
) {
  if (!(is.numeric(T0) && length(T0) == 1L && is.finite(T0))) {
    stop("T0 must be a finite numeric scalar.")
  }

  node_ids   <- character(0)
  node_times <- NULL
  node_roles <- NULL

  if (!is.null(nodes0)) {
    if (is.data.frame(nodes0)) {
      if (!("id" %in% names(nodes0))) stop("nodes0 must have column 'id'.")
      node_ids <- as.character(nodes0$id)

      if (anyNA(node_ids) || any(!nzchar(node_ids))) stop("nodes0$id contains NA/empty IDs.")
      if (anyDuplicated(node_ids)) stop("nodes0$id contains duplicates. Provide unique baseline nodes.")

      if ("time" %in% names(nodes0)) {
        node_times <- as.numeric(nodes0$time)
        if (length(node_times) != length(node_ids)) stop("nodes0$time must match nodes0$id length.")
        if (any(!is.finite(node_times))) stop("nodes0$time must be finite for all nodes.")
        if (any(node_times > T0)) stop("nodes0$time must be <= T0 for all baseline nodes.")
      }

      if ("role" %in% names(nodes0)) {
        node_roles <- as.character(nodes0$role)
        if (length(node_roles) != length(node_ids)) stop("nodes0$role must match nodes0$id length.")
        ok <- is.na(node_roles) | node_roles %in% c("perp", "event")
        if (!all(ok)) {
          bad <- unique(node_roles[!ok])
          stop("nodes0$role contains invalid values: ", paste(bad, collapse = ", "),
               ". Allowed: 'perp', 'event', or NA.")
        }
      }
    } else {
      node_ids <- as.character(nodes0)
      if (anyNA(node_ids) || any(!nzchar(node_ids))) stop("nodes0 contains NA/empty IDs.")
      if (anyDuplicated(node_ids)) stop("nodes0 contains duplicates. Provide unique baseline nodes.")
    }
  }

  if (is.null(edges0)) {
    edges0 <- data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
  } else {
    if (!is.data.frame(edges0)) stop("edges0 must be a data.frame with columns i, j.")
    if (!all(c("i", "j") %in% names(edges0))) stop("edges0 must have columns 'i' and 'j'.")
    edges0$i <- as.character(edges0$i)
    edges0$j <- as.character(edges0$j)

    if (anyNA(edges0$i) || anyNA(edges0$j) || any(!nzchar(edges0$i)) || any(!nzchar(edges0$j))) {
      stop("edges0 contains NA/empty node IDs.")
    }
    if (any(edges0$i == edges0$j)) stop("edges0 contains self-loops (not allowed).")
  }

  edge_nodes <- unique(c(edges0$i, edges0$j))
  edge_nodes <- edge_nodes[nzchar(edge_nodes)]

  if (implicit_nodes) {
    all_nodes <- unique(c(node_ids, edge_nodes))
  } else {
    all_nodes <- node_ids
    missing <- setdiff(edge_nodes, all_nodes)
    if (length(missing) > 0L) {
      stop("edges0 references nodes not present in nodes0: ", paste(missing, collapse = ", "))
    }
  }

  if (length(all_nodes) == 0L) {
    return(state_init())
  }

  born <- rep(T0, length(all_nodes))
  names(born) <- all_nodes
  if (!is.null(node_times)) born[node_ids] <- node_times

  st <- state_init()

  # Add baseline nodes, with optional roles if present.
  if (!is.null(node_roles)) {
    role_named <- stats::setNames(node_roles, node_ids)
    st <- state_add_nodes(st, all_nodes, t_k = T0, role = role_named[all_nodes])
  } else {
    st <- state_add_nodes(st, all_nodes, t_k = T0, role = NA_character_)
  }

  # Override born times if nodes0 provided explicit times (<=T0)
  st$born[all_nodes] <- born

  # Keep net vertex attribute in sync as well
  if (length(all_nodes) > 0L) {
    vids <- unname(st$idx[all_nodes])
    if (any(is.na(vids))) stop("make_state(): internal error: missing idx for some nodes.", call. = FALSE)
    network::set.vertex.attribute(st$net, "born", value = as.numeric(born[all_nodes]), v = as.integer(vids))
  }

  if (nrow(edges0) > 0L) {
    st <- state_add_edges(st, edges0, t_k = T0)
  }

  st
}
