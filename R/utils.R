#' Title
#'
#' @param events
#' @param net_init_fun
#' @param net0
#' @param net_step_fun
#'
#' @return
#' @export
#'
#' @examples
create_net_from_events <- function(events,
                                   net_init_fun = net_init,
                                   net0 = NULL,
                                   net_step_fun = net_add_event,
                                   debug = FALSE) {
  #
  arrivals <- nodes_by_event(events)
  edges <- edges_by_event(events)

  event_times <- events$times$t
  if (length(event_times) == 0L) return(0)

  net <- if (is.null(net0)) net_init_fun() else net0

  for (idx in seq_along(event_times)) {
    t_k <- event_times[idx]
    new_nodes <- arrivals[[idx]]
    new_edges <- edges[[idx]]

    net <- net_step_fun(net, new_nodes, new_edges, t_k)

    if (debug == TRUE) {
      plot(net)
      readline("Press Enter / Return to continue to next event...")
    }
  }

  return(net)
}

# tiny helper (so we don't need rlang)
`%||%` <- function(x, y) if (is.null(x)) y else x

# Internal: transform parameters to working scale
to_working <- function(p_list, transform) {
  stopifnot(is.list(p_list))

  x <- numeric(length(p_list))
  names(x) <- names(p_list)

  for (nm in names(p_list)) {
    val <- p_list[[nm]]
    tf  <- transform[[nm]] %||% "log"

    if (tf == "log") {
      if (!is.numeric(val) || length(val) != 1L || !is.finite(val) || val <= 0) {
        stop("Parameter '", nm, "' must be > 0 for log-transform.")
      }
      x[[nm]] <- log(val)
    } else if (tf == "none") {
      if (!is.numeric(val) || length(val) != 1L || !is.finite(val)) {
        stop("Parameter '", nm, "' must be a finite scalar.")
      }
      x[[nm]] <- val
    } else {
      stop("Unknown transform '", tf, "' for parameter '", nm, "'.")
    }
  }

  x
}

# Internal: transform parameters back to original scale
from_working <- function(x, transform) {
  stopifnot(is.numeric(x))

  p <- as.list(x)
  names(p) <- names(x)

  for (nm in names(x)) {
    tf <- transform[[nm]] %||% "log"
    if (tf == "log") {
      p[[nm]] <- exp(x[[nm]])
    } else if (tf == "none") {
      p[[nm]] <- x[[nm]]
    } else {
      stop("Unknown transform '", tf, "' for parameter '", nm, "'.")
    }
  }

  p
}


