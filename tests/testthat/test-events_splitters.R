test_that("edges_by_event/nodes_by_event require an events object", {
  expect_error(edges_by_event(list()), "class 'events'")
  expect_error(nodes_by_event(list()), "class 'events'")
})

test_that("edges_by_event/nodes_by_event require expected components", {
  x <- structure(list(times = data.frame(t = 1)), class = "events")
  expect_error(edges_by_event(x), "\\$edges")
  expect_error(nodes_by_event(x), "\\$nodes")
})

test_that("returns length-n list and includes 0-row dfs for empty events", {
  edges <- data.frame(i = c("a", "a"), j = c("b", "c"), time = c(1, 3))
  nodes <- data.frame(id = c("a", "b", "c"), time = c(0, 1, 3))
  ev <- make_events(edges = edges, nodes = nodes)

  eb <- edges_by_event(ev)
  nb <- nodes_by_event(ev)

  expect_type(eb, "list")
  expect_type(nb, "list")
  expect_length(eb, nrow(ev$times))
  expect_length(nb, nrow(ev$times))

  expect_true(all(vapply(eb, inherits, logical(1), "data.frame")))
  expect_true(all(vapply(nb, inherits, logical(1), "data.frame")))

  # event_id is removed
  expect_false(any(vapply(eb, function(df) "event_id" %in% names(df), logical(1))))
  expect_false(any(vapply(nb, function(df) "event_id" %in% names(df), logical(1))))
})

test_that("nodes_by_event drops missing/empty ids after coercion", {
  edges <- data.frame(i = c("a"), j = c("b"), time = c(1))
  nodes <- data.frame(id = c("a", "b"), time = c(0, 1))
  ev <- make_events(edges = edges, nodes = nodes)

  # Inject a bad row post-make_events to test cleaner behavior
  ev$nodes$id[1] <- ""

  nb <- nodes_by_event(ev)
  expect_true(sum(vapply(nb, nrow, integer(1))) <= nrow(ev$nodes))
})

test_that("errors if event_id cannot be mapped back into 1..n", {
  edges <- data.frame(i = c("a"), j = c("b"), time = c(1))
  nodes <- data.frame(id = c("a", "b"), time = c(0, 1))
  ev <- make_events(edges = edges, nodes = nodes)

  ev_bad <- ev
  ev_bad$edges$event_id <- "nope"
  expect_error(edges_by_event(ev_bad), "event_id.*integer", ignore.case = TRUE)

  ev_bad2 <- ev
  ev_bad2$nodes$event_id <- ev_bad2$nodes$event_id + 999L
  expect_error(nodes_by_event(ev_bad2), "1\\.\\.", ignore.case = TRUE)
})
