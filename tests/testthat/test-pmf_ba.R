# tests/testthat/test-pmf_ba.R

test_that("log_pmf_ba matches Bernoulli-over-possible-edges for one arrival", {
  st <- make_state_test_helper(
    deg  = c("1" = 1L, "2" = 3L),
    born = c("1" = 0.0, "2" = 0.0)
  )

  arrival_k <- "3"
  edges_k <- data.frame(i = "3", j = "2")

  # beta_edges=0, delta=0:
  # weights: w1=1, w2=3 => p1=1/4, p2=3/4
  # candidates: (3,1) missing, (3,2) present
  # log P = log(p2) + log(1 - p1) = log(3/4) + log(3/4) = 2*log(3/4)
  lp <- log_pmf_ba(st, arrival_k, edges_k, t_k = 1.0, params = list(beta_edges = 0), delta = 0)
  expect_equal(lp, 2 * log(3/4), tolerance = 1e-12)
})

test_that("log_pmf_ba uses beta_edges ageing via born times", {
  st <- make_state_test_helper(
    deg  = c("1" = 1L, "2" = 1L),
    born = c("1" = 0.0, "2" = 0.9)
  )

  arrival_k <- "3"
  edges_k <- data.frame(i = "3", j = "2")

  t_k <- 1.0
  beta_edges <- 1.0

  # weights w_i = (deg_i + delta) * exp(-beta_edges * (t - born_i)), delta=0
  w1 <- 1 * exp(-beta_edges * (t_k - st$born["1"]))
  w2 <- 1 * exp(-beta_edges * (t_k - st$born["2"]))
  p2 <- unname(w2 / (w1 + w2))
  p1 <- 1 - p2

  # observed: (3,2), missing: (3,1)
  lp_expected <- log(p2) + log(1 - p1)

  lp <- log_pmf_ba(st, arrival_k, edges_k, t_k = t_k, params = list(beta_edges = beta_edges), delta = 0)
  expect_equal(lp, lp_expected, tolerance = 1e-12)
})

test_that("log_pmf_ba returns correct non-edge contribution for edgeless arrival (>=2 old nodes)", {
  st <- make_state_test_helper(
    deg  = c("1" = 1L, "2" = 1L),
    born = c("1" = 0.0, "2" = 0.0)
  )

  lp <- log_pmf_ba(st, arrival_k = "3", edges_k = empty_edges(), t_k = 1.0, params = list(beta_edges = 0), delta = 0)

  # p_old = (1/2, 1/2) => log(1-1/2)+log(1-1/2) = 2*log(1/2)
  expect_equal(lp, 2 * log(0.5), tolerance = 1e-12)
})

test_that("log_pmf_ba is invariant to edge endpoint order (i,j swapped)", {
  st <- make_state_test_helper(
    deg  = c("1" = 2L, "2" = 6L),
    born = c("1" = 0.0, "2" = 0.0)
  )

  arrival_k <- "3"
  edges_a <- data.frame(i = "3", j = "2")
  edges_b <- data.frame(i = "2", j = "3") # swapped

  lp_a <- log_pmf_ba(st, arrival_k, edges_a, t_k = 1.0, params = list(beta_edges = 0), delta = 0)
  lp_b <- log_pmf_ba(st, arrival_k, edges_b, t_k = 1.0, params = list(beta_edges = 0), delta = 0)
  expect_equal(lp_a, lp_b, tolerance = 1e-12)
})

test_that("log_pmf_ba is invariant to edge row order", {
  st <- make_state_test_helper(
    deg  = c("1" = 1L, "2" = 3L),
    born = c("1" = 0.0, "2" = 0.0)
  )

  arrival_k <- "3"
  edges1 <- data.frame(i = c("3", "2"), j = c("2", "3"))  # same edge twice but swapped
  edges2 <- edges1[c(2, 1), , drop = FALSE]

  # NOTE: this test assumes your validator allows swapped orientation but disallows duplicates
  # As written, edges1 contains duplicates; so we instead use two distinct edges:
  edges1 <- data.frame(i = c("3", "3"), j = c("2", "1"))
  edges2 <- edges1[c(2, 1), , drop = FALSE]

  lp1 <- log_pmf_ba(st, arrival_k, edges1, t_k = 1.0, params = list(beta_edges = 0), delta = 0)
  lp2 <- log_pmf_ba(st, arrival_k, edges2, t_k = 1.0, params = list(beta_edges = 0), delta = 0)
  expect_equal(lp1, lp2, tolerance = 1e-12)
})

test_that("log_pmf_ba delta shifts probabilities as expected", {
  st <- make_state_test_helper(
    deg  = c("1" = 0L, "2" = 2L),
    born = c("1" = 0.0, "2" = 0.0)
  )

  arrival_k <- "3"
  edges_k <- data.frame(i = "3", j = "1") # attach to node 1

  # With delta=1: w1=1, w2=3 => p1=1/4, p2=3/4
  # Observed (3,1): log(p1) = log(1/4)
  # Missing  (3,2): log(1-p2) = log(1/4)
  lp <- log_pmf_ba(st, arrival_k, edges_k, t_k = 1.0, params = list(beta_edges = 0), delta = 1)
  expect_equal(lp, 2 * log(1/4), tolerance = 1e-12)
})

test_that("log_pmf_ba errors when sum of weights is zero (delta=0 and all degrees zero)", {
  st <- make_state_test_helper(
    deg  = c("1" = 0L, "2" = 0L),
    born = c("1" = 0.0, "2" = 0.0)
  )

  expect_error(
    log_pmf_ba(st, arrival_k = "3", edges_k = empty_edges(), t_k = 1.0, params = list(beta_edges = 0), delta = 0),
    "Invalid BA weights"
  )
})

test_that("log_pmf_ba allows zero degrees if delta > 0", {
  st <- make_state_test_helper(
    deg  = c("1" = 0L, "2" = 0L),
    born = c("1" = 0.0, "2" = 0.0)
  )

  # delta=1 => uniform: p=(1/2,1/2); no edges => 2*log(1/2)
  lp <- log_pmf_ba(st, arrival_k = "3", edges_k = empty_edges(), t_k = 1.0, params = list(beta_edges = 0), delta = 1)
  expect_equal(lp, 2 * log(0.5), tolerance = 1e-12)
})

test_that("log_pmf_ba is a proper PMF over all edge subsets (1 arrival × 2 old)", {
  st <- make_state_test_helper(
    deg  = c("1" = 1L, "2" = 3L),
    born = c("1" = 0.0, "2" = 0.0)
  )

  arrival_k <- "3"
  old_nodes <- names(st$deg)
  all_edges <- enumerate_candidate_edge_subsets(arrival_k, old_nodes)

  probs <- vapply(
    all_edges,
    function(edges_k) exp(log_pmf_ba(st, arrival_k, edges_k, t_k = 1.0, params = list(beta_edges = 0), delta = 0)),
    numeric(1)
  )

  expect_equal(sum(probs), 1.0, tolerance = 1e-12)
})

test_that("paper-ratio implies forced attachment when there is exactly 1 old node", {
  st <- make_state_test_helper(
    deg  = c("1" = 5L),
    born = c("1" = 0.0)
  )

  arrival_k <- "2"

  # No-edge mark should be impossible because p=1 when only one old node
  lp_no_edge <- log_pmf_ba(st, arrival_k, empty_edges(), t_k = 1.0, params = list(beta_edges = 0), delta = 0)
  expect_true(is.infinite(lp_no_edge) && lp_no_edge < 0)

  # Edge present should have log prob 0 (prob = 1)
  lp_edge <- log_pmf_ba(st, arrival_k, data.frame(i = "2", j = "1"), t_k = 1.0, params = list(beta_edges = 0), delta = 0)
  expect_equal(lp_edge, 0.0, tolerance = 1e-12)
})
