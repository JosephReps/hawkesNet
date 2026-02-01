#' Split edges into per-event batches
#'
#' @param events An `hg_events` object (from `as_events()` or `make_events()`).
#' @return A list of length `nrow(events$times)`. Element `[[k]]` is a data.frame
#'   with columns `i, j` (and any additional edge attributes such as `etype`, `w`,
#'   etc.) for event `k`. Events with zero edges return a 0-row data.frame with
#'   the same columns (minus `event_id`).
#' @export
edges_by_event <- function(events) {
  validate_events(events)

  n <- nrow(events$times)

  e <- events$edges
  keep_cols <- setdiff(names(e), "event_id")

  # 0-row template with correct columns
  empty <- e[0, keep_cols, drop = FALSE]

  # Ensure events with zero edges still have an entry in the list
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

  idx <- as.integer(names(spl))
  out[idx] <- spl

  out
}

