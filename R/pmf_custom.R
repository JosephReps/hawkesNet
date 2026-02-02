#' Template: log PMF for a mark kernel (one event)
#'
#' This is a minimal skeleton for implementing a custom mark kernel. A mark
#' kernel assigns a probability to the observed edge set for a single event at
#' time `t_k`, conditional on the pre-event network `net` and the event mark
#' (`new_nodes`, `new_edges`).
#'
#' @param net network object, the current network at time t_k
#' @param new_nodes data.frame with at least column "id". Any additional
#'                  columns will be treated as vertex attributes.
#' @param new_edges data.frame with at least columns "i" and "j". Any additional
#'                  columns will be treated as edge attributes.
#' @param t_k numeric scalar, the event timestamp
#' @param params Named list of model parameters (kernel-specific, any others will be ignored).
#' @param ... Optional kernel-specific arguments (e.g. tuning constants, caches).
#'
#' @return A list with:
#' \describe{
#'   \item{logp}{Scalar log probability of observing `new_edges` (and any other
#'   mark components) at this event.}
#'   \item{edge_probs}{Named numeric vector (or matrix) of attachment
#'   probabilities for eligible targets. This is returned for debugging /
#'   diagnostics and may be `NULL` if not meaningful for the kernel.}
#'   \item{model_cache}{Optional updated cache object (useful for expensive mark
#'   models). Include only if you accept/use a cache via `...`.}
#' }
#'
#' @examples
#' # See log_pmf_ba() for a working implementation.
custom_pmf <- function(
    net,
    new_nodes,
    new_edges,
    t_k,
    params,
    ...
) {
  stopifnot(inherits(net, "network"))
  stopifnot(length(t_k) == 1L, is.numeric(t_k), is.finite(t_k))

  # Validate your input here, eg:
  #     - check parameter values are appropriate1
  #     - check new_nodes and new_edges are possible under the kernel specs
  # I think it would be hygiene to create a separate "validate_custom...()"
  # but you could just throw it all in here.

  stopifnot(is.data.frame(new_edges), all(c("i", "j") %in% names(new_edges))) # Alter this to include your expected columns
  stopifnot(is.data.frame(new_nodes), all(c("id") %in% names(new_nodes))) # Alter this to include your expected columns

  # Next I would cover your edge cases & compute edge-candidates, e.g
  #     - There are no existing nodes in the network
  #     - There are existing nodes, but no elgible edge-candidates

  # And then finally, compute the actual edge probabilities
  edge_probs <- calculate_custom_edge_probs(net, t_k, params)

  # Placeholder
  logp <- NA_real_

  return(
    list(
      logp = as.numeric(logp),
      edge_probs = edge_probs
      )
    )
}
