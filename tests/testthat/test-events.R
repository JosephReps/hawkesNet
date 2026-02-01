# test-events.R

test_that("make_events() basic structure (nodes+edges)", {
  nodes <- data.frame(
    id   = c("1", "2", "3"),
    time = c(0, 1, 2)
  )
  edges <- data.frame(
    i    = c("1", "1"),
    j    = c("2", "3"),
    time = c(1, 2)
  )

  ev <- make_events(nodes = nodes, edges = edges)

  expect_s3_class(ev, "hg_events")
  expect_true(all(c("times", "nodes", "edges") %in% names(ev)))

  expect_true(is.data.frame(ev$times))
  expect_true(is.data.frame(ev$nodes))
  expect_true(is.data.frame(ev$edges))

  expect_true(all(c("event_id", "t") %in% names(ev$times)))
  expect_true(all(c("event_id", "id") %in% names(ev$nodes)))
  expect_true(all(c("event_id", "i", "j") %in% names(ev$edges)))

  # event_id grid is 1..n
  expect_identical(ev$times$event_id, seq_len(nrow(ev$times)))
  expect_true(is.unsorted(ev$times$t, strictly = TRUE) == FALSE)
})

test_that("make_events() maps times correctly", {
  nodes <- data.frame(
    id   = c("1", "2", "3"),
    time = c(0, 1, 2)
  )
  edges <- data.frame(
    i    = c("1", "1"),
    j    = c("2", "3"),
    time = c(1, 2)
  )

  ev <- make_events(nodes = nodes, edges = edges)

  # Nodes: 0,1,2 should map to event_ids 1,2,3
  expect_equal(ev$nodes$event_id[match("1", ev$nodes$id)], 1L)
  expect_equal(ev$nodes$event_id[match("2", ev$nodes$id)], 2L)
  expect_equal(ev$nodes$event_id[match("3", ev$nodes$id)], 3L)

  # Edges at time 1 -> event_id 2 ; time 2 -> event_id 3
  expect_true(all(ev$edges$event_id[ev$edges$i == "1" & ev$edges$j == "2"] == 2L))
  expect_true(all(ev$edges$event_id[ev$edges$i == "1" & ev$edges$j == "3"] == 3L))
})

test_that("make_events() preserves extra node/edge columns", {
  nodes <- data.frame(
    id    = c("1","2","3"),
    time  = c(0, 1, 2),
    role  = c("perp","event","perp"),
    score = c(0.1, 0.2, 0.3)
  )
  edges <- data.frame(
    i     = c("1"),
    j     = c("2"),
    time  = c(1),
    etype = "observed",
    w     = 5
  )

  ev <- make_events(nodes = nodes, edges = edges)

  # nodes extras preserved
  expect_true(all(c("role", "score") %in% names(ev$nodes)))
  expect_equal(ev$nodes$role[match("1", ev$nodes$id)], "perp")
  expect_equal(ev$nodes$role[match("2", ev$nodes$id)], "event")
  expect_equal(ev$nodes$role[match("3", ev$nodes$id)], "perp")
  expect_equal(ev$nodes$score[match("1", ev$nodes$id)], 0.1)
  expect_equal(ev$nodes$score[match("2", ev$nodes$id)], 0.2)
  expect_equal(ev$nodes$score[match("3", ev$nodes$id)], 0.3)

  # edges extras preserved
  expect_true(all(c("etype", "w") %in% names(ev$edges)))
  expect_equal(ev$edges$etype[1], "observed")
  expect_equal(ev$edges$w[1], 5)
})

test_that("implicit births fill NA for extra node columns", {
  nodes <- data.frame(
    id   = c("10","11"),
    time = c(1, 1),
    role = c("perp","perp")
  )
  edges <- data.frame(
    i    = c("10","10"),
    j    = c("11","12"),
    time = c(1, 2)
  )

  ev <- make_events(nodes = nodes, edges = edges, allow_implicit_birth = TRUE)

  # "12" should be created at time 2 with role NA
  idx_12 <- which(ev$nodes$id == "12")
  expect_length(idx_12, 1)
  expect_true(is.na(ev$nodes$role[idx_12]))
})

test_that("make_events() errors on duplicate node id", {
  nodes <- data.frame(id = c("1", "1"), time = c(0, 1))
  edges <- data.frame(i = character(0), j = character(0), time = numeric(0))
  expect_error(make_events(nodes = nodes, edges = edges), "at most once|Duplicate", ignore.case = TRUE)
})

test_that("make_events() errors on self-loop edge", {
  nodes <- data.frame(id = c("1"), time = c(0))
  edges <- data.frame(i = c("1"), j = c("1"), time = c(1))
  expect_error(make_events(nodes = nodes, edges = edges), "Self-loops|self", ignore.case = TRUE)
})

test_that("make_events() errors on edge before birth", {
  nodes <- data.frame(id = c("1", "2"), time = c(5, 0))
  edges <- data.frame(i = c("1"), j = c("2"), time = c(1))
  expect_error(make_events(nodes = nodes, edges = edges), "before.*born|before birth", ignore.case = TRUE)
})

test_that("validate_events(strict_cols=TRUE) rejects extra columns", {
  times <- data.frame(event_id = 1L, t = 0)
  nodes <- data.frame(event_id = 1L, id = "1", role = "perp")
  edges <- data.frame(event_id = integer(0), i = character(0), j = character(0))

  ev <- structure(list(times = times, nodes = nodes, edges = edges),
                  class = c("hg_events", "list"))

  expect_error(validate_events(ev, strict_cols = TRUE), "extra columns", ignore.case = TRUE)
})

test_that("validate_events() accepts extra columns by default", {
  times <- data.frame(event_id = 1L, t = 0)
  nodes <- data.frame(event_id = 1L, id = "1", role = "perp")
  edges <- data.frame(event_id = integer(0), i = character(0), j = character(0))

  ev <- structure(list(times = times, nodes = nodes, edges = edges),
                  class = c("hg_events", "list"))

  expect_true(isTRUE(validate_events(ev)))
})
