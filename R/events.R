# events.R

#' Create a marked network-growth event stream
#'
#' This package models a growing network as a sequence of **event times**.
#' At each time, two things may happen:
#' 1) Node arrivals (possibly zero)
#' 2) Edge additions (possibly zero)
#'
#' The event stream is constructed from timestamped `nodes` and `edges`.
#' Optionally, you may provide an initial baseline network state `state0`
#' (built by `make_state()`), which defines nodes that already exist at the
#' observation start. Baseline nodes in `state0` are treated as *known endpoints*
#' for validation, but they are NOT added as events.
#'
#' Extra columns in `nodes` and `edges` are preserved and carried through into the
#' returned `hg_events` object. In the returned object, `time` is replaced by
#' `event_id` (an integer index into `times`).
#'
#' @param nodes Optional data.frame with columns `id`, `time`.
#'   Additional columns are allowed and will be preserved.
#' @param edges data.frame with columns `i`, `j`, `time`.
#'   Additional columns are allowed and will be preserved.
#'   (May be empty but must be provided.)
#' @param allow_implicit_birth Logical; if TRUE (default), endpoints missing from `nodes`
#'   and missing from `state0` are treated as born at their first incident edge time.
#' @param allow_multi_edges Logical; if FALSE (default), an undirected edge may be added at most once.
#' @param state0 Optional initial state (e.g., from `make_state()`), used only to
#'   define baseline nodes and their birth times for validation. Must be a state
#'   list with at least `born` (named numeric) or `deg` (named integer).
#'
#' @return An object of class `hg_events` with data.frames `times`, `nodes`, `edges`.
#' @export
make_events <- function(nodes = NULL,
                        edges,
                        allow_implicit_birth = TRUE,
                        allow_multi_edges = FALSE,
                        state0 = state_init()) {
  if (is.null(edges)) stop("`edges` must be provided (may be empty).")

  edges <- as.data.frame(edges)
  if (!all(c("i", "j", "time") %in% names(edges))) {
    stop("`edges` must have columns: i, j, time")
  }

  if (is.null(nodes)) {
    nodes <- data.frame(id = character(0), time = numeric(0))
  } else {
    nodes <- as.data.frame(nodes)
    if (!all(c("id", "time") %in% names(nodes))) {
      stop("`nodes` must have columns: id, time")
    }
  }

  # --- coerce IDs to character keys ---
  nodes$id <- as.character(nodes$id)
  edges$i  <- as.character(edges$i)
  edges$j  <- as.character(edges$j)

  # --- time sanity ---
  nodes$time <- as.numeric(nodes$time)
  edges$time <- as.numeric(edges$time)

  if (length(nodes$time) > 0L && any(!is.finite(nodes$time))) stop("`nodes$time` must be finite numeric.")
  if (length(edges$time) > 0L && any(!is.finite(edges$time))) stop("`edges$time` must be finite numeric.")

  # --- basic id sanity (before adding implicit births) ---
  if (anyNA(nodes$id) || any(nodes$id == "")) stop("`nodes$id` contains missing/empty ids.")
  if (anyNA(edges$i) || anyNA(edges$j)) stop("`edges` contains missing endpoints.")
  if (any(edges$i == "" | edges$j == "")) stop("`edges` contains empty endpoint ids.")

  # --- extract baseline nodes/births from state0 (for validation only) ---
  known_nodes <- .known_nodes_from_state(state0)
  known_born  <- .known_born_from_state(state0, known_nodes)

  # enforce "one way to input": baseline nodes must NOT also appear as in-window arrivals
  if (length(known_nodes) > 0L && length(nodes$id) > 0L) {
    overlap <- intersect(unique(nodes$id), known_nodes)
    if (length(overlap) > 0L) {
      stop("Nodes appear in both `nodes` and `state0`: ",
           paste(overlap, collapse = ", "),
           ". Provide baseline nodes only in `state0`, and in-window arrivals only in `nodes`.")
    }
  }

  # --- implicit births (optional): fill nodes missing from nodes/state0 but present in edges ---
  if (allow_implicit_birth && nrow(edges) > 0L) {
    in_nodes <- unique(nodes$id)
    in_edges <- unique(c(edges$i, edges$j))
    missing  <- setdiff(in_edges, union(in_nodes, known_nodes))

    if (length(missing) > 0L) {
      first_time <- vapply(
        missing,
        function(id) min(edges$time[edges$i == id | edges$j == id]),
        numeric(1)
      )

      # Create implicit rows with the same schema as `nodes`.
      # Extra columns (if any) will be NA for implicit nodes.
      implicit <- nodes[0, , drop = FALSE]
      implicit <- implicit[rep(1, length(missing)), , drop = FALSE]
      implicit$id   <- missing
      implicit$time <- unname(first_time)

      nodes <- rbind(nodes, implicit)
    }
  }

  # --- node uniqueness: each id born once (in-window nodes only) ---
  if (any(duplicated(nodes$id))) stop("Each node id may appear at most once in `nodes`.")

  # --- canonicalize undirected edges + edge sanity ---
  if (nrow(edges) > 0L) {
    if (any(edges$i == edges$j)) stop("Self-loops (i == j) are not allowed.")

    a <- edges$i
    b <- edges$j
    swap <- a > b
    edges$i[swap] <- b[swap]
    edges$j[swap] <- a[swap]
  }

  # --- edge uniqueness over time (growth only) ---
  if (!allow_multi_edges && nrow(edges) > 0L) {
    key <- paste0(edges$i, "|", edges$j)
    if (any(duplicated(key))) {
      stop("Duplicate undirected edges found. If you want to allow repeated additions, set allow_multi_edges=TRUE.")
    }
  }

  # --- edges must reference known or born nodes, and occur after birth ---
  born_time <- stats::setNames(nodes$time, nodes$id)
  if (length(known_born) > 0L) {
    # names are unique by construction
    born_time <- c(born_time, known_born)
  }

  if (nrow(edges) > 0L) {
    present <- union(nodes$id, known_nodes)

    # If allow_implicit_birth=FALSE, this catches missing endpoints not in nodes or state0
    if (!allow_implicit_birth) {
      if (any(!(edges$i %in% present) | !(edges$j %in% present))) {
        stop("Some edge endpoints are not present in `nodes` or `state0` (and allow_implicit_birth=FALSE).")
      }
    } else {
      # even with implicit births, born_time should cover everything by now
      missing_endpoints <- setdiff(unique(c(edges$i, edges$j)), names(born_time))
      if (length(missing_endpoints) > 0L) {
        stop("Internal error: missing birth times for endpoints: ", paste(missing_endpoints, collapse = ", "))
      }
    }

    ti <- born_time[edges$i]
    tj <- born_time[edges$j]

    if (anyNA(ti) || anyNA(tj)) {
      stop("Internal error: birth time lookup returned NA for some endpoints.")
    }

    if (any(edges$time < ti) || any(edges$time < tj)) {
      stop("Some edges occur before one or both endpoints are born.")
    }
  }

  # --- build event grid as union of node+edge times ---
  all_times <- sort(unique(c(nodes$time, edges$time)))
  if (length(all_times) == 0L) stop("No events: both `nodes` and `edges` are empty.")

  times <- data.frame(event_id = seq_along(all_times), t = all_times)

  # --- map times to event_id using match() ---
  nodes$event_id <- match(nodes$time, all_times)
  edges$event_id <- if (nrow(edges) > 0L) match(edges$time, all_times) else integer(0)

  if (anyNA(nodes$event_id)) stop("Failed to map some node times onto the event grid.")
  if (anyNA(edges$event_id)) stop("Failed to map some edge times onto the event grid.")

  # --- preserve extra columns: replace `time` with `event_id` ---
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

  validate_events(events)
  events
}

# ---- internal helpers --------------------------------------------------------

.known_nodes_from_state <- function(state0) {
  if (is.null(state0)) return(character(0))

  if (!is.list(state0) || is.null(names(state0))) {
    stop("`state0` must be a state list (e.g., from make_state()/state_init()).")
  }

  if (!is.null(state0$born)) {
    ids <- names(state0$born)
  } else if (!is.null(state0$deg)) {
    ids <- names(state0$deg)
  } else {
    stop("`state0` must contain either `$born` (named numeric) or `$deg` (named integer).")
  }

  if (is.null(ids)) ids <- character(0)
  ids <- as.character(ids)

  if (anyNA(ids) || any(ids == "")) stop("`state0` contains missing/empty node ids.")
  if (anyDuplicated(ids)) stop("`state0` contains duplicated node ids.")

  ids
}

.known_born_from_state <- function(state0, known_nodes) {
  if (length(known_nodes) == 0L) return(stats::setNames(numeric(0), character(0)))

  if (!is.null(state0$born)) {
    b <- state0$born
    if (!is.numeric(b) || is.null(names(b))) stop("`state0$born` must be a named numeric vector.")
    if (any(!is.finite(b))) stop("`state0$born` must be finite for all baseline nodes.")
    # Ensure it covers known_nodes exactly
    missing <- setdiff(known_nodes, names(b))
    extra   <- setdiff(names(b), known_nodes)
    if (length(missing) > 0L) stop("`state0$born` is missing baseline nodes: ", paste(missing, collapse = ", "))
    if (length(extra) > 0L)   stop("`state0$born` has extra names not in baseline: ", paste(extra, collapse = ", "))
    return(b[known_nodes])
  }

  # If no born vector exists, we cannot do "edge after birth" checks relative to baseline.
  # We set baseline births to -Inf so any finite edge time is after birth.
  stats::setNames(rep(-Inf, length(known_nodes)), known_nodes)
}

#' Validate an events object
#'
#' @param events hg_events
#' @param strict_cols if TRUE, disallow extra columns beyond the required ones
#' @return Invisibly TRUE; otherwise errors.
#' @export
validate_events <- function(events, strict_cols = FALSE) {
  if (!is.list(events) || is.null(names(events))) {
    stop("`events` must be a named list with elements `times`, `nodes`, `edges`.")
  }
  if (!all(c("times", "nodes", "edges") %in% names(events))) {
    stop("`events` must contain elements named `times`, `nodes`, and `edges`.")
  }
  if (!is.data.frame(events$times)) stop("`times` must be a data.frame.")
  if (!is.data.frame(events$nodes)) stop("`nodes` must be a data.frame.")
  if (!is.data.frame(events$edges)) stop("`edges` must be a data.frame.")

  times <- events$times
  nodes <- events$nodes
  edges <- events$edges

  req_times <- c("event_id", "t")
  req_nodes <- c("event_id", "id")
  req_edges <- c("event_id", "i", "j")

  miss <- setdiff(req_times, names(times))
  if (length(miss)) stop("`times` is missing columns: ", paste(miss, collapse = ", "))
  miss <- setdiff(req_nodes, names(nodes))
  if (length(miss)) stop("`nodes` is missing columns: ", paste(miss, collapse = ", "))
  miss <- setdiff(req_edges, names(edges))
  if (length(miss)) stop("`edges` is missing columns: ", paste(miss, collapse = ", "))

  if (strict_cols) {
    extra <- setdiff(names(times), req_times)
    if (length(extra)) stop("`times` has extra columns not allowed: ", paste(extra, collapse = ", "))
    extra <- setdiff(names(nodes), req_nodes)
    if (length(extra)) stop("`nodes` has extra columns not allowed: ", paste(extra, collapse = ", "))
    extra <- setdiff(names(edges), req_edges)
    if (length(extra)) stop("`edges` has extra columns not allowed: ", paste(extra, collapse = ", "))
  }

  if (anyNA(times$event_id) || anyNA(times$t)) stop("`times` contains NA in `event_id` or `t`.")
  if (anyNA(nodes$event_id) || anyNA(nodes$id)) stop("`nodes` contains NA in `event_id` or `id`.")
  if (anyNA(edges$event_id) || anyNA(edges$i) || anyNA(edges$j)) stop("`edges` contains NA in `event_id`, `i`, or `j`.")

  if (!is.numeric(times$t) || any(!is.finite(times$t))) stop("`times$t` must be finite numeric.")
  if (!is.numeric(times$event_id) || any(times$event_id %% 1 != 0)) stop("`times$event_id` must be integer-ish.")
  if (!is.numeric(nodes$event_id) || any(nodes$event_id %% 1 != 0)) stop("`nodes$event_id` must be integer-ish.")
  if (!is.numeric(edges$event_id) || any(edges$event_id %% 1 != 0)) stop("`edges$event_id` must be integer-ish.")

  n <- nrow(times)
  if (n == 0L) stop("`times` must have at least one row.")
  if (!identical(as.integer(times$event_id), seq_len(n))) {
    stop("`times$event_id` must be exactly 1:nrow(times) in order (no gaps, no reordering).")
  }
  if (is.unsorted(times$t, strictly = TRUE)) {
    stop("`times$t` must be strictly increasing (no ties).")
  }

  # event_id references
  if (nrow(nodes) > 0L) {
    bad <- setdiff(unique(as.integer(nodes$event_id)), seq_len(n))
    if (length(bad)) stop("`nodes$event_id` contains ids not in `times$event_id`: ", paste(bad, collapse = ", "))
  }
  if (nrow(edges) > 0L) {
    bad <- setdiff(unique(as.integer(edges$event_id)), seq_len(n))
    if (length(bad)) stop("`edges$event_id` contains ids not in `times$event_id`: ", paste(bad, collapse = ", "))
  }

  # nodes: unique birth
  nodes$id <- as.character(nodes$id)
  if (any(nodes$id == "")) stop("Empty node ids not allowed.")
  if (any(duplicated(nodes$id))) stop("Duplicate node ids in `nodes`.")

  # edges: no self loops, ids non-empty
  edges$i <- as.character(edges$i)
  edges$j <- as.character(edges$j)
  if (any(edges$i == "" | edges$j == "")) stop("Empty edge endpoint ids not allowed.")
  if (any(edges$i == edges$j)) stop("Self-loops (i == j) are not allowed.")
  if (any(duplicated(edges[c("event_id", "i", "j")]))) stop("Duplicate (event_id,i,j) rows found.")

  invisible(TRUE)
}
