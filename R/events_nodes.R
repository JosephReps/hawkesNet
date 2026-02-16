#' Split nodes into per-event batches
#'
#' @param events An events object created by `make_events()`.
#'
#' @return A list with one element per event. Element `[[k]]` is a data.frame of
#'   nodes for event `k` with column `id` (plus any additional node attributes),
#'   excluding `event_id`. Events with zero nodes return a 0-row data.frame with
#'   the same columns.
#'
#' @details
#' Node ids `id` are coerced to character, and rows with missing or empty ids are
#' dropped.
#'
#' @export
nodes_by_event <- function(events) {
  if (!inherits(events, "events")) {
    stop("`events` must be an object of class 'events', typically created by `make_events()`.")
  }
  if (!is.list(events) || is.null(events$times) || is.null(events$nodes)) {
    stop("`events` must contain `$times` and `$nodes`. Did you create this with `make_events()`?")
  }
  if (!is.data.frame(events$times) || !is.data.frame(events$nodes)) {
    stop("`events$times` and `events$nodes` must be data.frames. Did you create this with `make_events()`?")
  }

  n <- nrow(events$times)

  nd <- events$nodes
  if (n == 0L) {
    return(list())
  }

  if (!"event_id" %in% names(nd)) {
    stop("`events$nodes` must contain an `event_id` column. Did you create this with `make_events()`?")
  }
  if (!"id" %in% names(nd)) {
    stop("`events$nodes` must contain column `id`. Did you create this with `make_events()`?")
  }

  keep_cols <- setdiff(names(nd), "event_id")
  empty <- nd[0, keep_cols, drop = FALSE]
  out <- replicate(n, empty, simplify = FALSE)

  if (nrow(nd) == 0L) return(out)

  spl <- split(nd[, keep_cols, drop = FALSE], nd$event_id)

  spl <- lapply(spl, function(df) {
    df$id <- as.character(df$id)
    ok <- !is.na(df$id) & nzchar(df$id)
    df <- df[ok, , drop = FALSE]
    rownames(df) <- NULL
    df
  })

  idx <- suppressWarnings(as.integer(names(spl)))
  if (anyNA(idx)) {
    stop("`events$nodes$event_id` must be integer-valued. Did you create this with `make_events()`?")
  }
  if (any(idx < 1L | idx > n)) {
    stop("`events$nodes$event_id` must be in 1..nrow(events$times). Did you create this with `make_events()`?")
  }

  out[idx] <- spl
  out
}
