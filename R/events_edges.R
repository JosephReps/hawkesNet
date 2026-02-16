#' Split edges into per-event batches
#'
#' @param events An events object created by `make_events()`.
#'
#' @return A list with one element per event. Element `[[k]]` is a data.frame of
#'   edges for event `k` with columns `i` and `j` (plus any additional edge
#'   attributes), excluding `event_id`. Events with zero edges return a 0-row
#'   data.frame with the same columns.
#'
#' @details
#' Edge endpoints `i` and `j` are coerced to character, and rows with missing or
#' empty endpoints are dropped.
#'
#' @export
edges_by_event <- function(events) {
  if (!inherits(events, "events")) {
    stop("`events` must be an object of class 'events', typically created by `make_events()`.")
  }
  if (!is.list(events) || is.null(events$times) || is.null(events$edges)) {
    stop("`events` must contain `$times` and `$edges`. Did you create this with `make_events()`?")
  }
  if (!is.data.frame(events$times) || !is.data.frame(events$edges)) {
    stop("`events$times` and `events$edges` must be data.frames. Did you create this with `make_events()`?")
  }

  n <- nrow(events$times)

  e <- events$edges
  if (n == 0L) {
    # Degenerate case: no events. Return an empty list.
    return(list())
  }

  if (!"event_id" %in% names(e)) {
    stop("`events$edges` must contain an `event_id` column. Did you create this with `make_events()`?")
  }
  if (!all(c("i", "j") %in% names(e))) {
    stop("`events$edges` must contain columns `i` and `j`. Did you create this with `make_events()`?")
  }

  keep_cols <- setdiff(names(e), "event_id")
  empty <- e[0, keep_cols, drop = FALSE]
  out <- replicate(n, empty, simplify = FALSE)

  if (nrow(e) == 0L) return(out)

  # Split full edge rows by event_id (keep all cols except event_id)
  spl <- split(e[, keep_cols, drop = FALSE], e$event_id)

  # Clean each df (i/j as character; drop NA/empty endpoints)
  spl <- lapply(spl, function(df) {
    df$i <- as.character(df$i)
    df$j <- as.character(df$j)
    ok <- !is.na(df$i) & nzchar(df$i) & !is.na(df$j) & nzchar(df$j)
    df <- df[ok, , drop = FALSE]
    rownames(df) <- NULL
    df
  })

  # Map split results back into 1..n slots
  idx <- suppressWarnings(as.integer(names(spl)))
  if (anyNA(idx)) {
    stop("`events$edges$event_id` must be integer-valued. Did you create this with `make_events()`?")
  }
  if (any(idx < 1L | idx > n)) {
    stop("`events$edges$event_id` must be in 1..nrow(events$times). Did you create this with `make_events()`?")
  }

  out[idx] <- spl
  out
}
