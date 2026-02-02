#' Add edges to the network object
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
net_add_edges <- function(net, new_edges, t_k) {
  # Input validation
  if (!inherits(net, "network")) {
    stop("net_add_edges(): net must be a 'network' object.", call. = FALSE)
  }
  if (!(is.numeric(t_k) && length(t_k) == 1L && is.finite(t_k))) {
    stop("net_add_edges(): t_k must be a single finite numeric.", call. = FALSE)
  }
  if (!is.data.frame(new_edges)) {
    stop("net_add_edges(): new_edges must be a data.frame", call. = FALSE)
  }
  if (!(all(c("i", "j") %in% names(new_edges)))) {
    stop("net_add_nodes(): new_edges must have columns 'i' and 'j'.", call. = FALSE)
  }

  # If there are no new edges, just return the network as-is
  if (nrow(new_edges) == 0L) return(net)

  # Ensure node ids are character keys for naming
  i_chr <- as.character(new_edges$i)
  j_chr <- as.character(new_edges$j)

  if (anyNA(i_chr) || anyNA(j_chr) || any(!nzchar(i_chr)) || any(!nzchar(j_chr))) {
    stop("net_add_edges(): new_edges contains NA/empty node IDs.", call. = FALSE)
  }

  # No duplicate edges within event (undirected)
  key_evt <- paste(pmin(i_chr, j_chr), pmax(i_chr, j_chr), sep = "|")
  if (anyDuplicated(key_evt)) {
    stop("net_add_edges(): duplicate edges within event.", call. = FALSE)
  }

  # Make sure all edge endpoints already exist
  existing <- network::network.vertex.names(net)
  if (is.null(existing)) existing <- character(0)

  nodes_evt <- unique(c(i_chr, j_chr))
  missing_nodes <- setdiff(nodes_evt, existing)
  if (length(missing_nodes) > 0L) {
    stop(
      "net_add_edges(): edge endpoints not in network. Call net_add_nodes() first. Missing nodes: ",
      paste(missing_nodes, collapse = ", "),
      call. = FALSE
    )
  }

  # Map vertex names -> vertex indices in the network
  # (network.vertex.names order matches vertex indices 1..n)
  idx_map <- stats::setNames(seq_along(existing), existing)
  tail_idx <- as.integer(idx_map[i_chr])
  head_idx <- as.integer(idx_map[j_chr])

  if (anyNA(tail_idx) || anyNA(head_idx) || any(tail_idx <= 0L) || any(head_idx <= 0L)) {
    stop("net_add_edges(): internal error mapping node ids to vertex indices.", call. = FALSE)
  }

  # VERY IMPORTANT, we are assuming that the edges passed do not already exist
  # in the network. If you wanted to add a check for that, this is the place to
  # do it.

  # Add the edges to the network
  m <- length(tail_idx)
  e_old <- network::network.edgecount(net)
  network::add.edges(net, tail = tail_idx, head = head_idx)
  e_new <- network::network.edgecount(net)

  if (e_new != e_old + m) {
    stop(
      "net_add_edges(): unexpected edgecount after adding edges. ",
      "Expected +", m, " edges, got +", (e_new - e_old), ".",
      call. = FALSE
    )
  }

  # Edge IDs of newly added edges (network increments edge IDs sequentially)
  new_eids <- seq.int(e_old + 1L, e_new)

  # Set edge attributes from extra columns
  attr_cols <- setdiff(names(new_edges), c("i", "j"))

  # Add a default time attribute unless user provided one already
  if (!("time" %in% attr_cols)) {
    network::set.edge.attribute(net, attrname = "time", value = rep(t_k, m), e = new_eids)
  }

  # Set any additional attributes (including user-provided 't', 'weight', etc.)
  if (length(attr_cols) > 0L) {
    for (nm in attr_cols) {
      network::set.edge.attribute(net, attrname = nm, value = new_edges[[nm]], e = new_eids)
    }
  }

  net
}
