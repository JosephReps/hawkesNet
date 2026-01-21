#' Initialize deterministic state before any events
#'
#' @return A state list with empty degree, birth-time, adjacency, role, and a cached network object.
#' @export
state_init <- function() {
  if (!requireNamespace("network", quietly = TRUE)) {
    stop("Package 'network' is required (for CS kernel / cached network state). Please install it.",
         call. = FALSE)
  }

  net <- network::network.initialize(0, directed = FALSE)

  list(
    deg  = structure(integer(0), names = character(0)),
    born = structure(numeric(0),  names = character(0)),
    adj  = structure(list(),      names = character(0)),
    role = structure(character(0), names = character(0)),

    # Cached network representation for CS kernels:
    net  = net,

    # Map node-id (character) -> vertex index in `net`
    idx = structure(integer(0), names = character(0))
  )
}
