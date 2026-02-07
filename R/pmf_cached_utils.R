# R/pmf_cached_utils.R
#
# Shared helpers for cache-based mark PMFs.
# Keep these small + boring: pure functions, no network mutation.

.get_event_times <- function(events) {
  stopifnot(is.list(events), !is.null(events$times), is.data.frame(events$times))
  if ("t" %in% names(events$times)) return(events$times$t)
  if ("time" %in% names(events$times)) return(events$times$time)
  stop("events$times must contain column 't' or 'time'.", call. = FALSE)
}

# Clamp probabilities away from 0/1 to avoid -Inf / NaNs in log terms.
.clamp_prob <- function(p, eps = 1e-12) {
  stopifnot(is.numeric(eps), length(eps) == 1L, is.finite(eps), eps > 0, eps < 0.5)
  pmin(pmax(p, eps), 1 - eps)
}

# Bernoulli-product loglik with a boolean mask is_obs
# ll = sum(is_obs * log(p) + (1-is_obs) * log(1-p))
.bern_ll <- function(p, is_obs) {
  stopifnot(length(p) == length(is_obs))
  sum(ifelse(is_obs, log(p), log1p(-p)))
}

# Undirected edge key helper (order-invariant)
.edge_key_undirected <- function(a, b) {
  a <- as.character(a); b <- as.character(b)
  lo <- ifelse(a <= b, a, b)
  hi <- ifelse(a <= b, b, a)
  paste0(lo, "|", hi)
}

.edge_keys_from_df_undirected <- function(edges_df) {
  if (is.null(edges_df) || nrow(edges_df) == 0L) return(character(0))
  stopifnot(is.data.frame(edges_df), all(c("i","j") %in% names(edges_df)))
  .edge_key_undirected(edges_df$i, edges_df$j)
}
