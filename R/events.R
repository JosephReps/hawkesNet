# events.R

#' Create a marked network-growth event stream
#'
#' Construct an event stream from timestamped node arrivals and edge additions.
#' Baseline nodes can be provided via `initial_network` for validation only (they
#' are not added as events).
#'
#' Extra columns in `nodes` and `edges` are preserved and passed to the internal
#' network object as vertex and edge attributes. The returned object uses
#' `event_id` (index into `times`) rather than raw `time` columns.
#'
#' ## Input contract (requirements)
#' The following requirements are enforced via errors (and covered by tests):
#'
#' **edges**
#'
#' - `edges` must be provided (may be empty), and be coercible to a data.frame.
#' - Must contain columns: `i`, `j`, `time`.
#' - `i` and `j` must be non-missing, non-empty (after coercion to character).
#' - No self-loops: `i != j`.
#' - `time` must be finite numeric.
#' - If `allow_multi_edges = FALSE` (default), each undirected edge `{i,j}` may
#'    appear at most once (duplicates after undirected canonicalization error).
#'
#' **nodes**
#' - If supplied, `nodes` must have columns: `id`, `time`.
#' - `id` must be non-missing, non-empty (after coercion to character) and unique.
#' - `time` must be finite numeric.
#'
#' **initial_network**
#' - `initial_network` must be a `network` object (or NULL).
#' - Baseline node ids (vertex names) must not also appear in `nodes$id`.
#' - If `initial_network` has vertex attribute `"time"`, it must be numeric
#'     and finite for all baseline nodes (used for edge-after-birth validation).
#'     If absent, baseline nodes are treated as born at `-Inf`.
#'
#' **cross-constraints**
#' - If `allow_implicit_birth = FALSE`, every edge endpoint must appear in
#'     `nodes$id` or `initial_network`.
#' - Every edge must occur after both endpoints are born (from `nodes` times
#'     and baseline `"born"` times, or `-Inf` if not provided).
#'
#' @param edges data.frame with columns `i`, `j`, `time`. Extra columns allowed.
#' @param nodes Optional data.frame with columns `id`, `time`. Extra columns allowed.
#' @param allow_implicit_birth Logical; if TRUE (default), endpoints missing from
#'   `nodes` and missing from `initial_network` are treated as born at their first
#'   incident edge time.
#' @param allow_multi_edges Logical; if FALSE (default), an undirected edge may be added at most once.
#' @param initial_network Optional baseline `network` object used for validation only.
#'
#' @return An object of class `events` with data.frames `times`, `nodes`, `edges`.
#' @export
make_events <- function(edges,
                        nodes = NULL,
                        allow_implicit_birth = TRUE,
                        allow_multi_edges = FALSE,
                        initial_network = net_init()) {
  # Make sure edges provided
  if (is.null(edges)) stop("`edges` must be provided (may be empty).")

  # Required columns check for edges
  edges <- as.data.frame(edges)
  if (!all(c("i", "j", "time") %in% names(edges))) {
    stop("`edges` must have columns: i, j, time")
  }

  # Required columns check for nodes, or create if NULL
  if (is.null(nodes)) {
    nodes <- data.frame(id = character(0), time = numeric(0))
  } else {
    nodes <- as.data.frame(nodes)
    if (!all(c("id", "time") %in% names(nodes))) {
      stop("`nodes` must have columns: id, time")
    }
  }

  # Coerce ID's to character
  nodes$id <- as.character(nodes$id)
  edges$i  <- as.character(edges$i)
  edges$j  <- as.character(edges$j)

  # Validate provided times
  nodes$time <- as.numeric(nodes$time)
  edges$time <- as.numeric(edges$time)
  if (length(nodes$time) > 0L && any(!is.finite(nodes$time))) stop("`nodes$time` must be finite numeric.")
  if (length(edges$time) > 0L && any(!is.finite(edges$time))) stop("`edges$time` must be finite numeric.")

  # Validate provided ID's
  if (anyNA(nodes$id) || any(nodes$id == "")) stop("`nodes$id` contains missing/empty ids.")
  if (anyNA(edges$i) || anyNA(edges$j)) stop("`edges` contains missing endpoints.")
  if (any(edges$i == "" | edges$j == "")) stop("`edges` contains empty endpoint ids.")

  # If an initial network has been provided, grab all the existing information
  known_nodes <- .known_nodes_from_network(initial_network)
  known_born  <- .known_born_from_network(initial_network, known_nodes)

  # Make sure baseline nodes DO NOT also appear in provided nodes
  if (length(known_nodes) > 0L && length(nodes$id) > 0L) {
    overlap <- intersect(unique(nodes$id), known_nodes)
    if (length(overlap) > 0L) {
      stop(
        "Nodes appear in both `nodes` and `initial_network`: ",
        paste(overlap, collapse = ", "),
        ". Provide baseline nodes only in `initial_network`, and in-window arrivals only in `nodes`."
      )
    }
  }

  # Calculate implicit births -
  # Fill nodes missing from nodes/initial_network but present in edges
  if (allow_implicit_birth && nrow(edges) > 0L) {
    # Grab the provided node ID's in edges and nodes
    in_nodes <- unique(nodes$id)
    in_edges <- unique(c(edges$i, edges$j))
    # Figure out which ID's are missing from nodes union initial network
    missing  <- setdiff(in_edges, union(in_nodes, known_nodes))

    # Grab the first time the missing nodes appear in edges
    if (length(missing) > 0L) {
      first_time <- vapply(
        missing,
        function(id) min(edges$time[edges$i == id | edges$j == id]),
        numeric(1)
      )

      # Make sure to warn the user that additional columns implicitly birthed
      # default have NA entries by default
      extra_cols <- setdiff(names(nodes), c("id", "time"))
      if (length(extra_cols) > 0L) {
        warning(
          "Created ", length(missing), " implicit node births. ",
          "Extra `nodes` columns filled with NA: ",
          paste(extra_cols, collapse = ", "),
          call. = FALSE
        )
      }

      n_missing <- length(missing)

      implicit <- nodes[rep(NA_integer_, n_missing), , drop = FALSE]

      implicit$id   <- missing
      implicit$time <- unname(first_time)

      nodes <- rbind(nodes, implicit)
    }
  }

  # Each node born once
  if (any(duplicated(nodes$id))) stop("Each node id may appear at most once in `nodes`.")

  # Default to the lower of the two ID's appearing first and check for
  # self-loops.
  # I think defaulting to this behaivour is a terrible idea, will probably
  # change at some point
  if (nrow(edges) > 0L) {
    if (any(edges$i == edges$j)) stop("Self-loops (i == j) are not allowed.")

    a <- edges$i
    b <- edges$j
    swap <- a > b
    edges$i[swap] <- b[swap]
    edges$j[swap] <- a[swap]
  }

  # Don't not allow duplicate edges if allow_multi_edges = FALSE
  if (!allow_multi_edges && nrow(edges) > 0L) {
    key <- paste0(edges$i, "|", edges$j)
    if (any(duplicated(key))) {
      stop("Duplicate undirected edges found. If you want to allow repeated additions, set allow_multi_edges=TRUE.")
    }
  }

  # Edges must reference known or born nodes, and occur after birth
  born_time <- stats::setNames(nodes$time, nodes$id)
  if (length(known_born) > 0L) {
    born_time <- c(born_time, known_born)
  }

  if (nrow(edges) > 0L) {
    present <- union(nodes$id, known_nodes)

    if (!allow_implicit_birth) {
      if (any(!(edges$i %in% present) | !(edges$j %in% present))) {
        stop("Some edge endpoints are not present in `nodes` or `initial_network` (and allow_implicit_birth=FALSE).")
      }
    } else {
      missing_endpoints <- setdiff(unique(c(edges$i, edges$j)), names(born_time))
      # If this triggers, something went wrong in implicit birth calculation
      if (length(missing_endpoints) > 0L) {
        stop("Internal error: missing birth times for endpoints: ", paste(missing_endpoints, collapse = ", "))
      }
    }

    # Validated these earlier
    ti <- born_time[edges$i]
    tj <- born_time[edges$j]

    bad <- which(edges$time < ti | edges$time < tj)
    if (length(bad) > 0L) {
      k <- bad[1]
      stop(
        "Some edges occur before one or both endpoints are born. ",
        "Example: (", edges$i[k], ", ", edges$j[k], ") at time ", edges$time[k],
        " but births are (", ti[k], ", ", tj[k], ")."
      )
    }
  }

  # Build event grid as union of node+edge times
  all_times <- sort(unique(c(nodes$time, edges$time)))
  if (length(all_times) == 0L) stop("No events: both `nodes` and `edges` are empty.")

  times <- data.frame(event_id = seq_along(all_times), t = all_times)

  # Map times to event_id using match()
  nodes$event_id <- match(nodes$time, all_times)
  edges$event_id <- if (nrow(edges) > 0L) match(edges$time, all_times) else integer(0)

  if (anyNA(nodes$event_id)) stop("Failed to map some node times onto the event grid.")
  if (anyNA(edges$event_id)) stop("Failed to map some edge times onto the event grid.")

  # Preserve extra columns: replace `time` with `event_id`
  nodes_keep <- setdiff(names(nodes), "time")
  nodes_keep <- c("event_id", "id", setdiff(nodes_keep, c("event_id", "id")))
  nodes_df <- nodes[, nodes_keep, drop = FALSE]

  edges_df <- if (nrow(edges) > 0L) {
    edges_keep <- setdiff(names(edges), "time")
    edges_keep <- c("event_id", "i", "j", setdiff(edges_keep, c("event_id", "i", "j")))
    edges[, edges_keep, drop = FALSE]
  } else {
    data.frame(event_id = integer(0), i = character(0), j = character(0))
  }

  events <- list(times = times, nodes = nodes_df, edges = edges_df)
  class(events) <- c("events", "list")

  events
}

# ---- internal helpers --------------------------------------------------------

.known_nodes_from_network <- function(initial_network) {
  if (is.null(initial_network)) return(character(0))

  if (!inherits(initial_network, "network")) {
    stop("`initial_network` must be a `network` object (e.g., from net_init()).")
  }

  ids <- network::network.vertex.names(initial_network)
  if (is.null(ids)) ids <- character(0)
  ids <- as.character(ids)

  if (anyNA(ids) || any(ids == "")) stop("`initial_network` contains missing/empty vertex names.")
  if (anyDuplicated(ids)) stop("`initial_network` contains duplicated vertex names.")

  ids
}

.known_born_from_network <- function(initial_network, known_nodes) {
  if (length(known_nodes) == 0L) return(stats::setNames(numeric(0), character(0)))
  if (is.null(initial_network)) {
    return(stats::setNames(rep(-Inf, length(known_nodes)), known_nodes))
  }

  # Optional vertex attribute "born"
  born <- network::get.vertex.attribute(initial_network, "born")

  if (is.null(born)) {
    return(stats::setNames(rep(-Inf, length(known_nodes)), known_nodes))
  }

  if (!is.numeric(born)) stop("`initial_network` vertex attribute 'born' must be numeric if provided.")

  # Ensure length matches vertices
  if (length(born) != length(known_nodes)) {
    stop("`initial_network` vertex attribute 'born' must have length equal to number of baseline vertices.")
  }

  # Require finite for all baseline nodes (contract); if you want NA -> -Inf, loosen this later.
  if (any(!is.finite(born))) stop("`initial_network` vertex attribute 'born' must be finite for all baseline nodes.")

  stats::setNames(as.numeric(born), known_nodes)
}
