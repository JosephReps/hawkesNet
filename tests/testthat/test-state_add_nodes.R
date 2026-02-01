# test-state_add_nodes.R

test_that("state_add_nodes adds new nodes with deg=0, born=t_k, empty adj, role=NA by default", {
  s <- state_init()
  s2 <- state_add_nodes(s, c(1, "2"), t_k = 3.5)

  expect_setequal(names(s2$deg),  c("1", "2"))
  expect_setequal(names(s2$born), c("1", "2"))
  expect_setequal(names(s2$adj),  c("1", "2"))
  expect_setequal(names(s2$role), c("1", "2"))

  expect_identical(unname(s2$deg), c(0L, 0L))
  expect_true(all(s2$born == 3.5))

  expect_true(all(vapply(s2$adj, is.character, logical(1))))
  expect_true(all(lengths(s2$adj) == 0L))

  expect_true(all(is.na(unname(s2$role))))
})

test_that("state_add_nodes is idempotent and ignores NA/empty ids (roles align to kept ids)", {
  s <- state_init()
  s <- state_add_nodes(s, c("A", NA, "", "A"), t_k = 1)

  expect_setequal(names(s$deg), "A")
  expect_identical(s$deg[["A"]], 0L)
  expect_identical(s$born[["A"]], 1)
  expect_true(is.na(s$role[["A"]]))

  s2 <- state_add_nodes(s, c("A"), t_k = 999, role = "perp")
  # A already exists so it should not be re-born nor have role overwritten
  expect_identical(s2$born[["A"]], 1)
  expect_true(is.na(s2$role[["A"]]))
})

test_that("state_add_nodes does not alter existing nodes' deg/born/adj/role", {
  s <- state_init()
  s <- state_add_nodes(s, c("A", "B"), t_k = 1, role = c("perp", "event"))
  s$deg[["A"]]  <- 7L
  s$adj[["A"]]  <- c("B")
  s$born[["A"]] <- 0.25
  s$role[["A"]] <- "perp"

  s2 <- state_add_nodes(s, c("A", "C"), t_k = 2, role = c("event", "event"))

  expect_identical(s2$deg[["A"]], 7L)
  expect_identical(s2$born[["A"]], 0.25)
  expect_identical(s2$adj[["A"]],  c("B"))
  expect_identical(s2$role[["A"]], "perp")  # unchanged

  expect_identical(s2$deg[["C"]],  0L)
  expect_identical(s2$born[["C"]], 2)
  expect_identical(s2$adj[["C"]],  character(0))
  expect_identical(s2$role[["C"]], "event")
})

test_that("state_add_nodes accepts role length-1 and recycles", {
  s <- state_init()
  s2 <- state_add_nodes(s, c("A", "B", "C"), t_k = 0, role = "perp")

  expect_identical(unname(s2$role[c("A", "B", "C")]), rep("perp", 3))
})

test_that("state_add_nodes accepts role vector matching node_ids length (including NA), and drops roles for dropped ids", {
  s <- state_init()
  ids  <- c("A", NA, "B", "", "C")
  role <- c("perp", "event", NA, "perp", "event")

  s2 <- state_add_nodes(s, ids, t_k = 1, role = role)

  expect_setequal(names(s2$deg), c("A", "B", "C"))
  expect_identical(s2$role[["A"]], "perp")
  expect_true(is.na(s2$role[["B"]]))
  expect_identical(s2$role[["C"]], "event")
})

test_that("state_add_nodes errors on invalid role values", {
  s <- state_init()
  expect_error(
    state_add_nodes(s, c("A", "B"), t_k = 0, role = c("perp", "banana")),
    "invalid role",
    ignore.case = TRUE
  )
})

test_that("state_add_nodes errors if role is non-scalar and length mismatch", {
  s <- state_init()
  expect_error(
    state_add_nodes(s, c("A", "B", "C"), t_k = 0, role = c("perp", "event")),
    "must match length",
    ignore.case = TRUE
  )
})
