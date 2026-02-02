#' Add nodes to the network object
#'
#' @inheritParams new-event-params
#'
#' @return Updated `network` object
#' @export
#'
#' @examples
#' # net <- network::network.initialize(0, directed = FALSE)
#' # net <- net_add_nodes(net, data.frame(id=c("a","b"), stringsAsFactors=FALSE), t_k=0)
#' # net <- net_add_edges(net, data.frame(i="a", j="b", w=1), t_k=0.5)
net_add_nodes <- function(net, new_nodes, t_k) {
  # Input validation
  if (!inherits(net, "network")) {
    stop("net_add_nodes(): net must be a 'network' object.", call. = FALSE)
  }
  if (!(is.numeric(t_k) && length(t_k) == 1L && is.finite(t_k))) {
    stop("net_add_nodes(): t_k must be a single finite numeric.", call. = FALSE)
  }
  if (!is.data.frame(new_nodes)) {
    stop("net_add_nodes(): new_nodes must be a data.frame", call. = FALSE)
  }
  if (!("id" %in% names(new_nodes))) {
    stop("net_add_nodes(): new_nodes must have column 'id'.", call. = FALSE)
  }
  if (anyNA(new_nodes$id)) {
    stop("net_add_nodes(): new_nodes$id contains NA.", call. = FALSE)
  }
  if (any(!nzchar(new_nodes$id))) {
    stop("net_add_nodes(): new_nodes$id contains empty strings.", call. = FALSE)
  }

  # If there are no new nodes, just return the network as-is
  if (length(new_nodes$id) == 0L) return(net)

  # Ensure uniqueness within the provided new_nodes df
  if (anyDuplicated(new_nodes$id)) {
    stop("net_add_nodes(): duplicate node ids provided in new_nodes.", call. = FALSE)
  }

  # Make sure new nodes do not already exist
  existing <- network::network.vertex.names(net)
  if (is.null(existing)) existing <- character(0)
  if (length(intersect(existing, new_nodes$id)) > 0) {
    stop("net_add_nodes(): node ids provided in new_nodes already exist in the network.", call. = FALSE)
  }

  # Construct vertex attributes
  new_nodes$time <- t_k
  names(new_nodes)[names(new_nodes) == "id"] <- "vertex.names"
  # This is some dumb shit tbh, I'm sure there is a good reason this shape is
  # required but fuck me it's annoying
  vat <- lapply(seq_len(nrow(new_nodes)), function(i) {
    as.list(new_nodes[i, , drop = FALSE])
  })

  # Add the new vertices
  network::add.vertices(net, nv = nrow(new_nodes), vattr = vat)

  return(net)
}

