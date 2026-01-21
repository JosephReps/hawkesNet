#' Split edges into per-event batches
#'
#' @param events An `hg_events` object (from `as_events()`).
#' @return A list of length `nrow(events$times)`. Element `[[k]]` is a data.frame
#'   with columns `i, j` containing all edges for event `k`. Events with zero
#'   edges return a 0-row data.frame.
#' @export
edges_by_event <- function(events) {
  validate_events(events)

  n <- nrow(events$times)

  # Ensure events with zero edges still have an entry in the list
  out <- replicate(
    n,
    data.frame(i = character(), j = character()),
    simplify = FALSE
  )

  e <- events$edges
  if (nrow(e) == 0L) return(out)

  # Split the edges data frame by event_id, but only keep the i & j columns
  spl <- split(e[, c("i", "j")], e$event_id)

  # Make sure to clear the row names
  spl <- lapply(spl, function(df) {
    rownames(df) <- NULL
    df
  })

  # event_id are 1..n (forced constraint on data), and split() assigns the
  # event_id as the name of each list element so we can just assign directly
  idx <- as.integer(names(spl))
  out[idx] <- spl

  out
}
