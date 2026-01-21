#' Split node arrivals into per-event batches
#'
#' @param events An `hg_events` object (from `make_events()` or `as_events()`).
#' @return A list of length `nrow(events$times)`. Element `[[k]]` is a character
#'   vector of node ids arriving at event `k`. Events with zero arrivals return
#'   `character(0)`.
#' @export
nodes_by_event <- function(events) {
  validate_events(events)
  n <- nrow(events$times)
  out <- replicate(n, character(0), simplify = FALSE)

  nd <- events$nodes
  if (nrow(nd) == 0L) return(out)

  spl <- split(as.character(nd$id), nd$event_id)
  idx <- as.integer(names(spl))
  out[idx] <- lapply(spl, function(x) {
    x <- as.character(x)
    x[!is.na(x) & nzchar(x)]
  })

  out
}
