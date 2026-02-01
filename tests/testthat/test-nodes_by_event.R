test_that("nodes_by_event returns arrivals per event", {
  nodes <- data.frame(id = c("a","b","c"), time = c(0, 0.2, 0.2))
  edges <- data.frame(i = c("a"), j = c("b"), time = c(0.2))
  ev <- make_events(nodes = nodes, edges = edges)
  nb <- nodes_by_event(ev)
  expect_length(nb, nrow(ev$times))
  expect_equal(nb[[1]], "a")
  expect_equal(sort(nb[[2]]), c("b","c"))
})
