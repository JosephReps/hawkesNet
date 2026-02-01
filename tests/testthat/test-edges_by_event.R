test_that("edges_by_event returns one batch per event (including empty events)", {
  nodes <- data.frame(id = c("1","2","3","4"), time = c(0, 0.1, 0.2, 0.3))
  edges <- data.frame(
    time = c(0.1, 0.2, 0.4, 0.4),
    i = c("1", "2", "1", "4"),
    j = c("2", "3", "3", "1")
  )

  ev <- make_events(nodes, edges)
  batches <- edges_by_event(ev)

  expect_length(batches, 5)

  expect_identical(batches[[1]], data.frame(i = character(), j = character()))
  expect_identical(batches[[2]], data.frame(i = "1", j = "2"))
  expect_identical(batches[[3]], data.frame(i = "2", j = "3"))
  expect_identical(batches[[4]], data.frame(i = character(), j = character()))
  expect_identical(batches[[5]], data.frame(i = c("1", "1"), j = c("3", "4")))
})
