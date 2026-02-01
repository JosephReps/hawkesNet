#' Add one event to the network (nodes, then edges)
#'
#' @param net network object
#' @param new_nodes data.frame of nodes to add (must have id; may have role)
#' @param new_edges data.frame of edges to add (must have i, j; may have attrs)
#' @param t numeric scalar time
#'
#' @return updated network object
#' @export
net_add_event <- function(net, new_nodes, new_edges, t) {
  if (!is.null(new_nodes) && nrow(new_nodes) > 0L) {
    net <- net_add_nodes(net, new_nodes, t)
  }
  if (!is.null(new_edges) && nrow(new_edges) > 0L) {
    net <- net_add_edges(net, new_edges, t)
  }
  net
}
