#' Split node arrivals into per-event batches
#'
#' @param events An `hg_events` object (from `make_events()` or `as_events()`).
#' @return A list of length `nrow(events$times)`. Element `[[k]]` is a data.frame
#'   of node arrivals for event `k` (must contain column `id`; may contain `role`,
#'   `score`, etc). Events with zero arrivals return a 0-row data.frame with the
#'   same columns (minus `event_id`).
#' @export
nodes_by_event <- function(events) {
  validate_events(events)

  n <- nrow(events$times)

  nd <- events$nodes
  keep_cols <- setdiff(names(nd), "event_id")

  # 0-row template with correct columns
  empty <- nd[0, keep_cols, drop = FALSE]

  # Ensure events with zero nodes still have an entry in the list
  out <- replicate(n, empty, simplify = FALSE)

  if (nrow(nd) == 0L) return(out)

  # Split full node rows by event_id (keep all cols except event_id)
  spl <- split(nd[, keep_cols, drop = FALSE], nd$event_id)

  # Clean each df (id as character; drop NA/empty)
  spl <- lapply(spl, function(df) {
    df$id <- as.character(df$id)
    df <- df[!is.na(df$id) & nzchar(df$id), , drop = FALSE]
    rownames(df) <- NULL
    df
  })

  idx <- as.integer(names(spl))
  out[idx] <- spl

  out
}

