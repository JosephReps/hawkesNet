# test-state_step.R

test_that("state_step adds explicit arrivals then edges", {
  s <- state_init()
  edges <- data.frame(i = "A", j = "B", stringsAsFactors = FALSE)

  s2 <- state_step(s, arrivals_k = c("A", "B"), edges_k = edges, t_k = 5, implicit_birth = TRUE)

  expect_setequal(names(s2$deg), c("A", "B"))
  expect_identical(s2$born[["A"]], 5)
  expect_identical(s2$born[["B"]], 5)
  expect_identical(s2$deg[["A"]], 1L)
  expect_identical(s2$deg[["B"]], 1L)

  # role defaults to NA when not provided
  expect_true(is.na(s2$role[["A"]]))
  expect_true(is.na(s2$role[["B"]]))
})

test_that("state_step supports role_k for arrivals", {
  s <- state_init()
  edges0 <- data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)

  s2 <- state_step(
    s,
    arrivals_k = c("E1", "P1", "P2"),
    edges_k = edges0,
    t_k = 1,
    implicit_birth = TRUE,
    role_k = c("event", "perp", "perp")
  )

  expect_identical(s2$role[["E1"]], "event")
  expect_identical(s2$role[["P1"]], "perp")
  expect_identical(s2$role[["P2"]], "perp")
})

test_that("state_step implicit_birth=TRUE births edge endpoints not in arrivals with role=NA", {
  s <- state_init()
  edges <- data.frame(i = "A", j = "B", stringsAsFactors = FALSE)

  s2 <- state_step(s, arrivals_k = "A", edges_k = edges, t_k = 2, implicit_birth = TRUE, role_k = "event")

  expect_setequal(names(s2$deg), c("A", "B"))
  expect_identical(s2$born[["A"]], 2)
  expect_identical(s2$born[["B"]], 2)
  expect_identical(s2$deg[["A"]], 1L)
  expect_identical(s2$deg[["B"]], 1L)

  expect_identical(s2$role[["A"]], "event")
  expect_true(is.na(s2$role[["B"]]))  # implicit birth => unknown role
})

test_that("state_step implicit_birth=FALSE errors if edges reference unseen endpoints (STRICT via state_add_edges)", {
  s <- state_init()
  edges <- data.frame(i = "A", j = "B", stringsAsFactors = FALSE)

  expect_error(
    state_step(s, arrivals_k = character(0), edges_k = edges, t_k = 2, implicit_birth = FALSE),
    "endpoints not in state|Missing nodes|Call state_add_nodes",
    ignore.case = TRUE
  )
})

test_that("state_step with empty arrivals and empty edges is identity", {
  s <- state_init()
  edges0 <- data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)

  s2 <- state_step(s, arrivals_k = character(0), edges_k = edges0, t_k = 1)

  expect_identical(s2, s)
})
