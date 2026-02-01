# test-make_state.R

test_that("make_state() returns an empty state when no baseline nodes/edges", {
  st <- make_state(nodes0 = NULL, edges0 = NULL, T0 = 0)

  expect_true(is.list(st))
  expect_true(all(c("deg", "born", "adj", "role") %in% names(st)))
  expect_equal(length(st$deg),  0L)
  expect_equal(length(st$born), 0L)
  expect_equal(length(st$adj),  0L)
  expect_equal(length(st$role), 0L)
})

test_that("make_state() assigns birth times correctly", {
  # vector nodes0 -> all born at T0
  st1 <- make_state(nodes0 = c(1, 2, 3), T0 = 5)
  expect_equal(unname(st1$born), rep(5, 3))
  expect_setequal(names(st1$born), as.character(c(1, 2, 3)))
  expect_true(all(is.na(unname(st1$role)))) # no roles provided

  # data.frame nodes0 with time -> uses those times (must be <= T0)
  st2 <- make_state(nodes0 = data.frame(id = c("a", "b"), time = c(0, 4)), T0 = 4)
  expect_equal(st2$born[["a"]], 0)
  expect_equal(st2$born[["b"]], 4)
})

test_that("make_state() preserves roles if nodes0 has role column", {
  nodes0 <- data.frame(
    id   = c("e1", "p1", "p2"),
    time = c(0, 0, 0),
    role = c("event", "perp", "perp"),
    stringsAsFactors = FALSE
  )
  st <- make_state(nodes0 = nodes0, T0 = 0)

  expect_identical(st$role[["e1"]], "event")
  expect_identical(st$role[["p1"]], "perp")
  expect_identical(st$role[["p2"]], "perp")
})

test_that("make_state() validates role values if provided", {
  nodes0 <- data.frame(
    id   = c("x"),
    time = c(0),
    role = c("banana"),
    stringsAsFactors = FALSE
  )

  expect_error(
    make_state(nodes0 = nodes0, T0 = 0),
    "invalid",
    ignore.case = TRUE
  )
})

test_that("make_state() validates nodes0 input strictly", {
  expect_error(make_state(T0 = NA_real_), "T0 must be a finite numeric scalar")

  expect_error(
    make_state(nodes0 = data.frame(time = 0), T0 = 0),
    "nodes0 must have column 'id'"
  )

  expect_error(
    make_state(nodes0 = data.frame(id = c("x", "x")), T0 = 0),
    "duplicates"
  )

  expect_error(
    make_state(nodes0 = data.frame(id = c("x", ""), time = c(0, 0)), T0 = 0),
    "NA/empty"
  )

  expect_error(
    make_state(nodes0 = data.frame(id = c("x"), time = c(1)), T0 = 0),
    "<= T0"
  )
})

test_that("make_state() validates edges0, and handles implicit vs explicit baseline nodes", {
  # bad edges0 schema
  expect_error(
    make_state(nodes0 = c("1"), edges0 = list(i = "1", j = "2"), T0 = 0),
    "edges0 must be a data.frame"
  )
  expect_error(
    make_state(nodes0 = c("1"), edges0 = data.frame(i = "1", k = "2"), T0 = 0),
    "must have columns 'i' and 'j'"
  )

  edges0 <- data.frame(i = c("1"), j = c("2"))

  # implicit_nodes = TRUE => edges0 endpoints may introduce nodes not in nodes0
  st_imp <- make_state(nodes0 = NULL, edges0 = edges0, T0 = 0, implicit_nodes = TRUE)
  expect_setequal(names(st_imp$born), c("1", "2"))
  expect_true(all(is.na(unname(st_imp$role)))) # no roles provided

  # implicit_nodes = FALSE => edges0 endpoints must be listed in nodes0
  expect_error(
    make_state(nodes0 = "1", edges0 = edges0, T0 = 0, implicit_nodes = FALSE),
    "edges0 references nodes not present in nodes0"
  )
})
