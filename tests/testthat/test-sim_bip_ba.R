# test-sim_ba_bip.R

test_that("edge_probs_ba_bip only returns probs on old perps", {
  st <- state_init()
  st <- state_add_nodes(st, c("1", "2", "3", "4"), t_k = 0, role = c("perp", "event", "perp", "event"))
  st$deg[c("1","3")] <- c(2L, 5L)

  p <- edge_probs_ba_bip(st, t_k = 1, beta_edges = 0, delta = 0.1)

  expect_setequal(names(p), c("1", "3"))
  expect_true(is.numeric(p))
  expect_true(all(p >= 0 & p <= 1))
  expect_equal(sum(p), 1, tolerance = 1e-12)
})

test_that("sim_mark_ba_bip returns 1 event + K perps and deterministic edges to new perps", {
  st <- state_init()
  st <- state_add_nodes(st, c("1", "2", "3"), t_k = 0, role = c("perp", "perp", "event"))
  st <- state_add_edges(st, data.frame(i = "1", j = "3"), t_k = 0)

  # Make K deterministic by duplicating the seed/rpois call
  set.seed(123)
  K_expected <- stats::rpois(1L, 3)

  set.seed(123)
  mk <- sim_mark_ba_bip(st, t_k = 5, params = list(beta_edges = 0.2, lambda_new = 3), delta = 0.1)

  arrivals <- mk$arrivals
  edges <- mk$edges

  expect_true(is.character(arrivals))
  expect_equal(length(arrivals), 1L + K_expected)

  new_event <- arrivals[1]
  new_perps <- arrivals[-1]

  # deterministic edges must be present: (new_event, each new_perp)
  if (length(new_perps) > 0L) {
    key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "|")
    edge_keys <- key(edges$i, edges$j)
    det_keys  <- key(rep(new_event, length(new_perps)), new_perps)
    expect_true(all(det_keys %in% edge_keys))
  }

  # all edges must involve new_event (in this generator)
  if (nrow(edges) > 0L) {
    expect_true(all(edges$i == new_event | edges$j == new_event))
  }
})

test_that("sim_mark_ba_bip with lambda_new=0 and no old perps yields 0 edges", {
  st <- state_init()
  st <- state_add_nodes(st, c("2"), t_k = 0, role = "event")  # no perps exist

  set.seed(1)
  mk <- sim_mark_ba_bip(st, t_k = 1, params = list(beta_edges = 0, lambda_new = 0), delta = 0.01)

  expect_equal(length(mk$arrivals), 1L)
  expect_equal(nrow(mk$edges), 0L)
})


# test_that("sim_mark_ba_bip never attaches to old events (Bernoulli targets are old perps only)", {
#   st <- state_init()
#   st <- state_add_nodes(st, c("1", "2", "3", "4"), t_k = 0, role = c("perp", "event", "perp", "event"))
#   st$deg[c("1","3")] <- c(10L, 1L) # make old perps high-ish degree
#
#   set.seed(42)
#   mk <- sim_mark_ba_bip(st, t_k = 2, params = list(beta_edges = 0, lambda_new = 0), delta = 0.01)
#
#   new_event <- mk$arrivals[1]
#   edges <- mk$edges
#
#   if (nrow(edges) > 0L) {
#     other <- ifelse(edges$i == new_event, edges$j, edges$i)
#     expect_true(all(other %in% c("1", "3")))  # only old perps
#   }
# })
