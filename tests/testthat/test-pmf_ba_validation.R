# tests/testthat/test-pmf_ba_validation.R

test_that("log_pmf_ba errors if arrival_k is missing, empty, or length != 1", {
  st <- make_state_test_helper(
    deg  = c("1" = 1L),
    born = c("1" = 0.0)
  )

  expect_error(
    log_pmf_ba(st, arrival_k = character(0), edges_k = empty_edges(), t_k = 1.0,
               params = list(beta_edges = 0.1)),
    "length\\(arrival_k\\) == 1"
  )

  expect_error(
    log_pmf_ba(st, arrival_k = c("2", "3"), edges_k = empty_edges(), t_k = 1.0,
               params = list(beta_edges = 0.1)),
    "length\\(arrival_k\\) == 1"
  )

  expect_error(
    log_pmf_ba(st, arrival_k = "", edges_k = empty_edges(), t_k = 1.0,
               params = list(beta_edges = 0.1)),
    "nzchar"
  )
})

test_that("log_pmf_ba returns 0 for empty edge set when no old nodes exist", {
  st <- state_init()
  st$deg  <- setNames(integer(0), character(0))
  st$born <- setNames(numeric(0), character(0))
  st$adj  <- list()

  lp <- log_pmf_ba(st, arrival_k = "1", edges_k = empty_edges(), t_k = 1.0,
                   params = list(beta_edges = 0.1))
  expect_equal(lp, 0.0, tolerance = 1e-12)
})

test_that("log_pmf_ba errors if no old nodes but edges are observed", {
  st <- state_init()
  st$deg  <- setNames(integer(0), character(0))
  st$born <- setNames(numeric(0), character(0))
  st$adj  <- list()

  edges_k <- data.frame(i = "1", j = "2")
  expect_error(
    log_pmf_ba(st, arrival_k = "1", edges_k = edges_k, t_k = 1.0,
               params = list(beta_edges = 0.1)),
    "no eligible existing nodes|no eligible old nodes|no eligible"
  )
})

test_that("log_pmf_ba errors on new-new and old-old edges", {
  st <- make_state_test_helper(
    deg  = c("1" = 1L, "2" = 1L),
    born = c("1" = 0.0, "2" = 0.0)
  )

  # new-new: only allowed node is arrival_k, so (3,3) is invalid (not new-old)
  expect_error(
    log_pmf_ba(st, arrival_k = "3", edges_k = data.frame(i = "3", j = "3"), t_k = 1.0,
               params = list(beta_edges = 0.1)),
    "non-BA edges|non-BA"
  )

  # old-old: edge entirely among old nodes
  expect_error(
    log_pmf_ba(st, arrival_k = "3", edges_k = data.frame(i = "1", j = "2"), t_k = 1.0,
               params = list(beta_edges = 0.1)),
    "non-BA edges|non-BA"
  )
})

test_that("log_pmf_ba errors if edge endpoint not in {arrival_k} ∪ old_nodes", {
  st <- make_state_test_helper(
    deg  = c("1" = 1L),
    born = c("1" = 0.0)
  )

  edges_k <- data.frame(i = "2", j = "X") # X is neither old nor arrival
  expect_error(
    log_pmf_ba(st, arrival_k = "2", edges_k = edges_k, t_k = 1.0,
               params = list(beta_edges = 0.1)),
    "endpoint not in|edge endpoint"
  )
})

test_that("log_pmf_ba errors on duplicate new-old edges within the event", {
  st <- make_state_test_helper(
    deg  = c("1" = 1L),
    born = c("1" = 0.0)
  )

  arrival_k <- "2"
  edges_k <- data.frame(i = c("2", "2"), j = c("1", "1"))

  expect_error(
    log_pmf_ba(st, arrival_k, edges_k, t_k = 1.0,
               params = list(beta_edges = 0.1)),
    "Duplicate BA edges|duplicate"
  )
})
