# tests/testthat/test-sim_ba.R

test_that("sim_mark_ba returns 0-row edges when there are no old nodes", {
  st <- state_init()
  st$deg  <- setNames(integer(0), character(0))
  st$born <- setNames(numeric(0), character(0))
  st$adj  <- list()

  mk <- sim_mark_ba(st, t_k = 1.0, new_node = "1",
                       params = list(beta_edges = 0.1), delta = 0.5)

  expect_true(is.list(mk))
  expect_true(all(c("arrivals","edges") %in% names(mk)))
  expect_true(all(c("i","j") %in% names(mk$edges)))
  expect_equal(nrow(mk$edges), 0L)
})

test_that("sim_mark_ba errors if new_node is not a single nonempty ID", {
  st <- make_state_test_helper(
    deg  = c("1" = 1L),
    born = c("1" = 0.0)
  )

  expect_error(sim_mark_ba(st, t_k = 1.0, new_node = character(0),
                           params = list(beta_edges = 0.1), delta = 0.5))
  expect_error(sim_mark_ba(st, t_k = 1.0, new_node = c("2","3"),
                           params = list(beta_edges = 0.1), delta = 0.5))
  expect_error(sim_mark_ba(st, t_k = 1.0, new_node = "",
                           params = list(beta_edges = 0.1), delta = 0.5))
})

test_that("sim_mark_ba returns all candidate edges when there is exactly one old node", {
  st <- make_state_test_helper(
    deg  = c("1" = 5L),
    born = c("1" = 0.0)
  )

  # Under paper-ratio, with exactly one old node, p_old = 1
  mk <- sim_mark_ba(st, t_k = 1.0, new_node = "2",
                       params = list(beta_edges = 0.1), delta = 0.5)

  expect_equal(nrow(mk$edges), 1L)
  expect_equal(mk$edges$i, "2")
  expect_equal(mk$edges$j, "1")
})

test_that("sim_mark_ba output is BA-valid (new-old only, no duplicates) for random draws", {
  st <- make_state_test_helper(
    deg  = c("1" = 1L, "2" = 3L, "3" = 2L),
    born = c("1" = 0.0, "2" = 0.0, "3" = 0.0)
  )

  arrival <- "4"
  mk <- sim_mark_ba(st, t_k = 1.0, new_node = arrival,
                       params = list(beta_edges = 0.1), delta = 0.5)

  expect_true(is.data.frame(mk$edges))
  expect_true(all(c("i","j") %in% names(mk$edges)))

  if (nrow(mk$edges) > 0L) {
    # orientation
    expect_true(all(mk$edges$i == arrival))
    expect_true(all(mk$edges$j %in% names(st$deg)))

    # no duplicate (new,old) pairs
    keys <- paste0(mk$edges$i, "|", mk$edges$j)
    expect_equal(anyDuplicated(keys), 0L)

    # validator (new signature)
    expect_silent(validate_edges_ba(mk$edges, arrival_k = arrival, old_nodes = names(st$deg)))
  } else {
    # empty edge set is always valid; validator returns empty endpoints
    expect_silent(validate_edges_ba(mk$edges, arrival_k = arrival, old_nodes = names(st$deg)))
  }
})

test_that("sim_mark_ba empirical distribution matches log_pmf_ba on a tiny case", {
  set.seed(1)

  st <- make_state_test_helper(
    deg  = c("1" = 1L, "2" = 3L),
    born = c("1" = 0.0, "2" = 0.0)
  )

  arrival <- "3"
  old_nodes <- names(st$deg)

  # Enumerate all 4 possible edge subsets for candidates {(3,1),(3,2)}
  all_marks <- enumerate_candidate_edge_subsets(arrival, old_nodes)

  # Theoretical probs from log_pmf_ba
  p_theory <- vapply(
    all_marks,
    function(edges_k) exp(log_pmf_ba(st, arrival_k = arrival, edges_k = edges_k,
                                     t_k = 1.0, params = list(beta_edges = 0.1), delta = 0)),
    numeric(1)
  )

  # Monte Carlo sample
  B <- 20000L
  counts <- integer(length(all_marks))

  key_for_edges <- function(edges_k) {
    if (nrow(edges_k) == 0L) return("empty")
    paste(sort(edges_k$j), collapse = ",")
  }

  keys <- vapply(all_marks, key_for_edges, character(1))
  key_to_idx <- setNames(seq_along(keys), keys)

  for (b in seq_len(B)) {
    mk <- sim_mark_ba(st, t_k = 1.0, new_node = arrival,
                     params = list(beta_edges = 0.1), delta = 0)
    k <- key_for_edges(mk$edges)
    counts[[ key_to_idx[[k]] ]] <- counts[[ key_to_idx[[k]] ]] + 1L
  }

  p_emp <- counts / B

  expect_equal(p_emp, p_theory, tolerance = 0.02)
})
