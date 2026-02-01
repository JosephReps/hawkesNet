# test-state_init.R

test_that("state_init returns empty deterministic state", {
  st <- state_init()

  expect_identical(names(st), c("deg", "born", "adj", "role", "net", "idx"))

  expect_identical(st$deg,  structure(integer(0),  names = character(0)))
  expect_identical(st$born, structure(numeric(0),  names = character(0)))
  expect_identical(st$adj,  structure(list(),      names = character(0)))
  expect_identical(st$role, structure(character(0), names = character(0)))
})

test_that("state_init returns correctly typed empty structures with names", {
  s <- state_init()

  expect_true(is.list(s))
  expect_true(all(c("deg", "born", "adj", "role") %in% names(s)))

  expect_type(s$deg, "integer")
  expect_type(s$born, "double")
  expect_true(is.list(s$adj))
  expect_type(s$role, "character")

  expect_identical(length(s$deg),  0L)
  expect_identical(length(s$born), 0L)
  expect_identical(length(s$adj),  0L)
  expect_identical(length(s$role), 0L)

  expect_true(!is.null(names(s$deg)))
  expect_true(!is.null(names(s$born)))
  expect_true(!is.null(names(s$adj)))
  expect_true(!is.null(names(s$role)))

  expect_identical(names(s$deg),  character(0))
  expect_identical(names(s$born), character(0))
  expect_identical(names(s$adj),  character(0))
  expect_identical(names(s$role), character(0))
})
