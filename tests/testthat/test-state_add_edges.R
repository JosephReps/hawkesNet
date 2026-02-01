# test-state_add_edges.R

test_that("state_add_edges updates degrees and undirected adjacency for existing nodes", {
  s <- state_init()
  s <- state_add_nodes(s, c("A", "B", "C"), t_k = 0)

  edges <- data.frame(i = c("A", "A"), j = c("B", "C"), stringsAsFactors = FALSE)
  s2 <- state_add_edges(s, edges, t_k = 1, update_adj = TRUE)

  expect_identical(s2$deg[["A"]], 2L)
  expect_identical(s2$deg[["B"]], 1L)
  expect_identical(s2$deg[["C"]], 1L)

  expect_true("B" %in% s2$adj[["A"]])
  expect_true("C" %in% s2$adj[["A"]])
  expect_true("A" %in% s2$adj[["B"]])
  expect_true("A" %in% s2$adj[["C"]])

  # role should not be touched
  expect_true(all(is.na(unname(s2$role[c("A", "B", "C")]))))
})

test_that("state_add_edges respects update_adj = FALSE", {
  s <- state_init()
  s <- state_add_nodes(s, c("A", "B"), t_k = 0)

  edges <- data.frame(i = "A", j = "B", stringsAsFactors = FALSE)
  s2 <- state_add_edges(s, edges, t_k = 1, update_adj = FALSE)

  expect_identical(s2$deg[["A"]], 1L)
  expect_identical(s2$deg[["B"]], 1L)
  expect_identical(s2$adj[["A"]], character(0))
  expect_identical(s2$adj[["B"]], character(0))
})

test_that("state_add_edges with empty edges returns identical state", {
  s <- state_init()
  s <- state_add_nodes(s, c("A", "B"), t_k = 0)

  edges0 <- data.frame(i = character(0), j = character(0), stringsAsFactors = FALSE)
  s2 <- state_add_edges(s, edges0, t_k = 1)

  expect_identical(s2, s)
})

test_that("state_add_edges errors if any endpoint is missing from state (STRICT)", {
  s <- state_init()
  s <- state_add_nodes(s, "A", t_k = 0)

  edges <- data.frame(i = "A", j = "B", stringsAsFactors = FALSE)

  expect_error(
    state_add_edges(s, edges, t_k = 1, update_adj = FALSE),
    "endpoints not in state|Missing nodes|Call state_add_nodes",
    ignore.case = TRUE
  )
})
