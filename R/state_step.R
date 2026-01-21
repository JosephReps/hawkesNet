#' Advance state by one event (arrivals then edges)
#'
#' @param state deterministic state (see `state_init`)
#' @param arrivals_k character vector of node ids arriving at this event
#' @param edges_k data.frame with columns `i`, `j` (may be 0-row)
#' @param t_k numeric scalar event time
#' @param implicit_birth Logical; if TRUE, any edge endpoint not already present
#'   in state is treated as born at `t_k`.
#' @param role_k optional character roles for arrivals_k (length 1 or length(arrivals_k)).
#' @return updated state
#' @export
state_step <- function(state, arrivals_k, edges_k, t_k, implicit_birth = TRUE, role_k = NULL) {
  # First add explicit arrivals (with optional roles)
  state <- state_add_nodes(state, arrivals_k, t_k, role = role_k)

  # Then infer any other new nodes based on the new edges (roles unknown)
  if (implicit_birth && nrow(edges_k) > 0L) {
    endpoints <- unique(c(as.character(edges_k$i), as.character(edges_k$j)))
    missing <- setdiff(endpoints, names(state$deg))
    state <- state_add_nodes(state, missing, t_k, role = NA_character_)
  }

  # Lastly, add edges
  state_add_edges(state, edges_k, t_k)
}
