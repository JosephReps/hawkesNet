#' Initialize the network-state object
#'
#' @return a network object
#' @export
#'
#' @examples
net_init <- function() {
  net <- network::network.initialize(0, directed = FALSE)

  return(net)
}
