# Debug helpers -------------------------------------------------------------

.debug_prepare_full_layout <- function(full_net, niter = 500, pad = 0.05) {
  if (!requireNamespace("sna", quietly = TRUE)) {
    stop("Package 'sna' is required for debug plotting (install.packages('sna')).")
  }
  coords_full <- sna::gplot.layout.fruchtermanreingold(
    full_net,
    layout.par = list(niter = niter)
  )
  rownames(coords_full) <- network::network.vertex.names(full_net)

  xr <- range(coords_full[, 1], finite = TRUE)
  yr <- range(coords_full[, 2], finite = TRUE)
  dx <- diff(xr); dy <- diff(yr)
  if (!is.finite(dx) || dx == 0) dx <- 1
  if (!is.finite(dy) || dy == 0) dy <- 1

  xlim_full <- xr + c(-pad, pad) * dx
  ylim_full <- yr + c(-pad, pad) * dy

  list(coords_full = coords_full, xlim_full = xlim_full, ylim_full = ylim_full)
}

.debug_draw_info_panel <- function(t_k, mark_ll, lambda_g, ll, new_nodes, new_edges, edge_probs = NULL, top_n = 15) {
  graphics::plot.new()
  graphics::par(usr = c(0, 1, 0, 1))
  y <- 0.97
  line <- function(txt, adj = 0, cex = 0.95, font = 1) {
    graphics::text(0.02, y, labels = txt, adj = adj, cex = cex, font = font)
    y <<- y - 0.045
  }

  line(sprintf("t = %.6f", t_k), font = 2)
  line(sprintf("log q(mark) = %.6f", mark_ll))
  line(sprintf("log lambda_g = %.6f", log(lambda_g)))
  line(sprintf("ll (cum) = %.6f", ll), font = 2)
  y <- y - 0.02

  nn <- if (is.null(new_nodes)) 0 else nrow(new_nodes)
  ne <- if (is.null(new_edges)) 0 else nrow(new_edges)
  line(sprintf("new_nodes: %d", nn))
  line(sprintf("new_edges: %d", ne))
  y <- y - 0.02

  # Show top candidate edge probabilities if provided
  if (!is.null(edge_probs) && is.data.frame(edge_probs) && nrow(edge_probs) > 0 && FALSE) {
    # Expect columns: i, j, p (names may differ slightly)
    cn <- names(edge_probs)
    pcol <- if ("p" %in% cn) "p" else if ("prob" %in% cn) "prob" else NULL
    if (!is.null(pcol)) {
      ord <- order(edge_probs[[pcol]], decreasing = TRUE)
      edge_probs2 <- edge_probs[ord, , drop = FALSE]
      k <- min(top_n, nrow(edge_probs2))
      line(sprintf("Top %d candidate probs:", k), font = 2, cex = 0.9)
      for (r in seq_len(k)) {
        row <- edge_probs2[r, , drop = FALSE]
        i <- if ("i" %in% cn) as.character(row[["i"]]) else if ("tail" %in% cn) as.character(row[["tail"]]) else "?"
        j <- if ("j" %in% cn) as.character(row[["j"]]) else if ("head" %in% cn) as.character(row[["head"]]) else "?"
        pr <- as.numeric(row[[pcol]])
        line(sprintf("  %s -> %s : %.4g", i, j, pr), cex = 0.85)
        if (y < 0.05) break
      }
    }
  }
}

.debug_draw_net_panel <- function(net, coords_full, xlim_full, ylim_full,
                                  new_nodes, new_edges,
                                  edge_probs = NULL,
                                  highlight_prob = TRUE,
                                  prob_cex_scale = 3) {

  # vertex styling
  vcol <- rep("grey30", network::network.size(net))
  vcex <- rep(1, network::network.size(net))
  names(vcol) <- network::network.vertex.names(net)
  names(vcex) <- network::network.vertex.names(net)

  if (!is.null(new_nodes) && nrow(new_nodes) > 0 && "id" %in% names(new_nodes)) {
    ids <- as.character(new_nodes$id)
    vcol[ids] <- "forestgreen"
    vcex[ids] <- 1.5
  }

  # edge styling (highlight edges incident to new nodes)
  ne <- network::network.edgecount(net)
  ecol <- rep("grey60", ne)
  elwd <- rep(1, ne)

  if (!is.null(new_nodes) && nrow(new_nodes) > 0 && "id" %in% names(new_nodes)) {
    ids <- as.character(new_nodes$id)
    idx_new <- unlist(lapply(ids, function(id) network::get.edgeIDs(net, v = id)), use.names = FALSE)
    idx_new <- idx_new[is.finite(idx_new)]
    if (length(idx_new)) {
      ecol[idx_new] <- "forestgreen"
      elwd[idx_new] <- 5
    }
  }

  vn_now <- network::network.vertex.names(net)
  graphics::plot(
    net,
    coord = coords_full[vn_now, , drop = FALSE],
    xlim = xlim_full,
    ylim = ylim_full,
    vertex.col = vcol,
    vertex.cex = vcex,
    edge.col = ecol,
    edge.lwd = elwd,
    displaylabels = FALSE
  )

  # Optionally annotate candidate edge probs (for the current step)
  if (isTRUE(highlight_prob) && !is.null(edge_probs) && length(edge_probs) != 1 && FALSE) {
    cn <- names(edge_probs)

    # Only label a manageable number to avoid clutter: top 30
    ord <- order(edge_probs, decreasing = TRUE)
    edge_probs2 <- edge_probs[ord, drop = FALSE]
    edge_probs2 <- edge_probs2[seq_len(min(30, nrow(edge_probs2))), drop = FALSE]

    # THIS IS CURRENTLY WRONG
    # EDGE PROBS IS A NAMED VECTOR
    # IT IS THE EDGE PROBABILITIES OF THE NEW NODE CONNECTING TO EACH VECTOR
    # SO EDGE_PROBS HAS LENGTH OF NUM VERTICES - NUMBER NEW NODES
    # IT HAS NAMES OF VERTEX NAMES

    # Determine endpoint column names
    # icol <- if ("i" %in% cn) "i" else if ("tail" %in% cn) "tail" else NULL
    # jcol <- if ("j" %in% cn) "j" else if ("head" %in% cn) "head" else NULL
    # if (!is.null(icol) && !is.null(jcol)) {
    #   xy <- coords_full
    #   for (r in seq_len(nrow(edge_probs2))) {
    #     i <- as.character(edge_probs2[[icol]][r])
    #     j <- as.character(edge_probs2[[jcol]][r])
    #     pr <- as.numeric(edge_probs2[[pcol]][r])
    #     if (!is.finite(pr) || pr <= 0) next
    #     if (!(i %in% rownames(xy) && j %in% rownames(xy))) next
    #
    #     mx <- (xy[i, 1] + xy[j, 1]) / 2
    #     my <- (xy[i, 2] + xy[j, 2]) / 2
    #     graphics::text(mx, my, labels = sprintf("%.2g", pr), cex = 0.75)
    #   }
    # }

  }
}


.debug_overlay_box <- function(lines,
                               x0 = 0.02, y0 = 0.98,
                               width = 0.10,
                               line_h = 0.032,
                               padding = 0.012,
                               bg = NULL,
                               border = "grey30",
                               cex = 0.7) {
  if (is.null(lines) || length(lines) == 0) return(invisible(NULL))

  # Allow either:
  #  - character vector (all font=1)
  #  - data.frame with columns: text (or txt), font (optional), cex (optional)
  if (is.data.frame(lines)) {
    txt_col  <- if ("text" %in% names(lines)) "text" else if ("txt" %in% names(lines)) "txt" else NULL
    if (is.null(txt_col)) stop(".debug_overlay_box(): data.frame 'lines' must have column 'text' (or 'txt').", call. = FALSE)
    txt  <- as.character(lines[[txt_col]])
    fnt  <- if ("font" %in% names(lines)) as.integer(lines[["font"]]) else rep.int(1L, length(txt))
    cexv <- if ("cex"  %in% names(lines)) as.numeric(lines[["cex"]])  else rep.int(cex, length(txt))
  } else {
    txt  <- as.character(lines)
    fnt  <- rep.int(1L, length(txt))
    cexv <- rep.int(cex, length(txt))
  }

  # Try to use a semi-transparent background if supported by device
  if (is.null(bg)) {
    bg <- tryCatch(grDevices::adjustcolor("white", alpha.f = 0.80),
                   error = function(e) "white")
  }

  n <- length(txt)
  height <- padding * 2 + n * line_h
  x1 <- min(x0 + width, 0.33)
  y1 <- y0
  y2 <- max(y0 - height, 0.01)

  graphics::rect(x0, y2, x1, y1, col = bg, border = border, lwd = 1)

  y <- y1 - padding - line_h/2
  for (i in seq_len(n)) {
    if (is.na(txt[i]) || txt[i] == "") {
      y <- y - line_h
      next
    }
    graphics::text(x0 + padding, y, labels = txt[i], adj = c(0, 0.5),
                   cex = cexv[i], font = fnt[i])
    y <- y - line_h
    if (y < y2) break
  }
  invisible(NULL)
}

.debug_make_info_lines <- function(t_k, mark_ll, lambda_g, ll, new_nodes, new_edges, event_id = NULL) {
  nn <- if (is.null(new_nodes)) 0 else nrow(new_nodes)
  ne <- if (is.null(new_edges)) 0 else nrow(new_edges)

  # Helper to align "Label:" column neatly
  kv_block <- function(labels, values) {
    labels <- as.character(labels)
    values <- as.character(values)
    w <- max(nchar(labels), na.rm = TRUE)
    sprintf(paste0("%-", w, "s  %s"), paste0(labels, ":"), values)
  }

  event_labels <- c("Event ID", "Event time", "# new nodes", "# new edges")
  event_vals   <- c(
    if (is.null(event_id) || is.na(event_id)) "" else as.character(event_id),
    sprintf("%.7f", t_k),
    sprintf("%d", nn),
    sprintf("%d", ne)
  )

  ll_labels <- c("log q(mark)", "log lambda_g", "loglik (cumul)")
  ll_vals   <- c(
    sprintf("%.6f", mark_ll),
    sprintf("%.6f", log(lambda_g)),
    sprintf("%.6f", ll)
  )

  lines_txt <- c(
    "Event Info",
    kv_block(event_labels, event_vals),
    "",
    "Likelihood Info",
    kv_block(ll_labels, ll_vals)
  )

  # Style: bold section headers
  data.frame(
    text = lines_txt,
    font = ifelse(lines_txt %in% c("Event Info", "Likelihood Info"), 2L, 1L),
    stringsAsFactors = FALSE
  )
}

.debug_plot_step <- function(net, layout, t_k, mark_ll, lambda_g, ll, new_nodes, new_edges, event_id, edge_probs = NULL) {

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar), add = TRUE)

  # Single centered network plot
  graphics::par(mar = c(0.5, 0.5, 0.5, 0.5))
  .debug_draw_net_panel(
    net,
    layout$coords_full, layout$xlim_full, layout$ylim_full,
    new_nodes, new_edges,
    edge_probs = edge_probs,
    highlight_prob = FALSE
  )

  # Overlay text box in normalized device coordinates
  graphics::par(new = TRUE)
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), asp = NA)

  lines <- .debug_make_info_lines(t_k, mark_ll, lambda_g, ll, new_nodes, new_edges, event_id)
  .debug_overlay_box(lines, x0 = 0.02, y0 = 0.98, width = 0.50)

  graphics::par(new = FALSE)
}
