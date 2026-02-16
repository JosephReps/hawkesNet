test_that("make_events errors when edges is NULL or missing required cols", {
  expect_error(make_events(NULL), "`edges` must be provided")
  expect_error(make_events(data.frame(i="a", j="b")), "must have columns: i, j, time")
})

test_that("make_events validates nodes columns and uniqueness", {
  expect_error(
    make_events(data.frame(i="a", j="b", time=1), nodes=data.frame(id="a")),
    "`nodes` must have columns: id, time"
  )

  expect_error(
    make_events(data.frame(i="a", j="b", time=1),
                nodes=data.frame(id=c("a","a"), time=c(0,0.5))),
    "at most once"
  )
})

test_that("make_events validates times are finite numeric", {
  expect_error(
    make_events(data.frame(i="a", j="b", time=Inf)),
    "edges\\$time.*finite"
  )
  expect_error(
    make_events(data.frame(i="a", j="b", time=1),
                nodes=data.frame(id="a", time=NaN)),
    "nodes\\$time.*finite"
  )
})

test_that("make_events validates ids: non-missing, non-empty; no self-loops", {
  expect_error(make_events(data.frame(i=NA, j="b", time=1)), "missing endpoints")
  expect_error(make_events(data.frame(i="", j="b", time=1)), "empty endpoint")
  expect_error(make_events(data.frame(i="a", j="a", time=1)), "Self-loops")
})

test_that("make_events blocks duplicate undirected edges unless allow_multi_edges=TRUE", {
  e <- data.frame(i=c("a","b"), j=c("b","a"), time=c(1,2))
  expect_error(make_events(e), "Duplicate undirected edges")

  out <- make_events(e, allow_multi_edges = TRUE)
  expect_s3_class(out, "events")
  expect_equal(nrow(out$edges), 2)
})

test_that("make_events can do implicit births and warns if extra node cols exist", {
  # nodes has an extra column -> implicit births should warn it becomes NA
  nodes <- data.frame(id="a", time=0, foo="x")
  edges <- data.frame(i="a", j="b", time=1)

  expect_warning(
    out <- make_events(edges, nodes = nodes, allow_implicit_birth = TRUE),
    "Extra `nodes` columns filled with NA"
  )

  expect_true("b" %in% out$nodes$id)
  # implicit node b should have foo = NA
  b_row <- out$nodes[out$nodes$id == "b", , drop = FALSE]
  expect_true(is.na(b_row$foo))
})

test_that("make_events errors when allow_implicit_birth=FALSE and endpoints not present", {
  edges <- data.frame(i="a", j="b", time=1)
  nodes <- data.frame(id="a", time=0)

  expect_error(
    make_events(edges, nodes=nodes, allow_implicit_birth = FALSE),
    "allow_implicit_birth=FALSE"
  )
})

test_that("make_events errors if an edge occurs before a node birth", {
  nodes <- data.frame(id=c("a","b"), time=c(5,0))
  edges <- data.frame(i="a", j="b", time=1)
  expect_error(make_events(edges, nodes=nodes), "occur before.*born")
})

test_that("make_events returns consistent event grid and event_id mapping", {
  nodes <- data.frame(id=c("a","b"), time=c(0,2))
  edges <- data.frame(i="a", j="b", time=2)

  out <- make_events(edges, nodes=nodes)

  expect_true(all(c("times","nodes","edges") %in% names(out)))
  expect_true(all(c("event_id","t") %in% names(out$times)))
  expect_true(all(out$nodes$event_id %in% out$times$event_id))
  expect_true(all(out$edges$event_id %in% out$times$event_id))
})
