# test-pmf_ba_bip.R

test_that("log_pmf_ba_bip matches hand calculation on a tiny case", {
  # State: two old perps with degrees 1 and 3, and one old event
  st <- state_init()
  st <- state_add_nodes(st, c("A", "B", "E0"), t_k = 0, role = c("perp", "perp", "event"))
  st$deg[c("A","B")] <- c(1L, 3L)
  st$born[c("A","B")] <- c(0, 0)

  t_k <- 1
  params <- list(beta_edges = 0, lambda_new = 2)
  delta <- 1

  # Mark: one new event "X", K_new = 2 new perps "P1","P2"
  arrivals <- c("X", "P1", "P2")
  role_k   <- c("event", "perp", "perp")

  # Deterministic edges to new perps (X-P1, X-P2)
  # Plus a Bernoulli edge only to old perp A (not to B)
  edges <- data.frame(
    i = c("X", "X", "X"),
    j = c("P1", "P2", "A"),
    stringsAsFactors = FALSE
  )

  # Hand:
  # old perp probs with beta=0, age=1, delta=1:
  # wA=(1+1)=2, wB=(3+1)=4 => pA=1/3, pB=2/3
  pA <- 2 / (2 + 4)
  pB <- 4 / (2 + 4)

  lp_expected <- stats::dpois(2, 2, log = TRUE) + log(pA) + log(1 - pB)

  lp <- log_pmf_ba_bip(st, arrivals, edges, t_k, params, role_k, delta = delta)

  expect_equal(lp, lp_expected, tolerance = 1e-12)
})

test_that("log_pmf_ba_bip includes Poisson(K_new) term correctly", {
  st <- state_init()
  st <- state_add_nodes(st, c("A", "E0"), t_k = 0, role = c("perp", "event"))
  st$deg["A"] <- 1L
  st$born["A"] <- 0

  params <- list(beta_edges = 0, lambda_new = 3)
  delta <- 0.1
  t_k <- 1

  # Same old-perp edge pattern (edge to A always), different K_new
  arrivals1 <- c("X", "P1")
  role1 <- c("event", "perp")
  edges1 <- data.frame(i = c("X", "X"), j = c("P1", "A"))

  arrivals2 <- c("Y")  # K_new = 0
  role2 <- "event"
  edges2 <- data.frame(i = "Y", j = "A")  # Bernoulli edge to old perp

  lp1 <- log_pmf_ba_bip(st, arrivals1, edges1, t_k, params, role1, delta = delta)
  lp2 <- log_pmf_ba_bip(st, arrivals2, edges2, t_k, params, role2, delta = delta)

  # Difference should equal dpois(1)-dpois(0) because Bernoulli part is identical
  expect_equal(lp1 - lp2,
               stats::dpois(1, 3, log = TRUE) - stats::dpois(0, 3, log = TRUE),
               tolerance = 1e-12)
})

test_that("validate_edges_ba_bip errors if role_k is missing (strict)", {
  st <- state_init()
  st <- state_add_nodes(st, c("A", "E0"), t_k = 0, role = c("perp", "event"))

  expect_error(
    log_pmf_ba_bip(st, arrivals_k = "X", edges_k = data.frame(i=character(0), j=character(0)),
                   t_k = 1, params = list(beta_edges = 0, lambda_new = 1),
                   role_k = NULL),
    "requires `role_k`",
    ignore.case = TRUE
  )
})

test_that("validate_edges_ba_bip errors if arrivals_k does not contain exactly one event", {
  st <- state_init()
  st <- state_add_nodes(st, c("A", "E0"), t_k = 0, role = c("perp", "event"))

  expect_error(
    log_pmf_ba_bip(st, arrivals_k = c("X","Y"), edges_k = data.frame(i=character(0), j=character(0)),
                   t_k = 1, params = list(beta_edges = 0, lambda_new = 1),
                   role_k = c("event","event")),
    "exactly one",
    ignore.case = TRUE
  )

  expect_error(
    log_pmf_ba_bip(st, arrivals_k = c("X","Y"), edges_k = data.frame(i=character(0), j=character(0)),
                   t_k = 1, params = list(beta_edges = 0, lambda_new = 1),
                   role_k = c("perp","perp")),
    "exactly one",
    ignore.case = TRUE
  )
})

test_that("validate_edges_ba_bip errors on missing deterministic edges", {
  st <- state_init()
  st <- state_add_nodes(st, c("A", "E0"), t_k = 0, role = c("perp", "event"))
  st$deg["A"] <- 1L
  st$born["A"] <- 0

  arrivals <- c("X", "P1", "P2")
  role_k   <- c("event", "perp", "perp")

  # Only one of the deterministic edges is present
  edges <- data.frame(i = c("X"), j = c("P1"), stringsAsFactors = FALSE)

  expect_error(
    log_pmf_ba_bip(st, arrivals, edges, t_k = 1, params = list(beta_edges = 0, lambda_new = 1),
                   role_k = role_k, delta = 0.1),
    "missing deterministic",
    ignore.case = TRUE
  )
})

test_that("validate_edges_ba_bip errors if an edge is not incident to the new event", {
  st <- state_init()
  st <- state_add_nodes(st, c("A", "B", "E0"), t_k = 0, role = c("perp", "perp", "event"))
  st$deg[c("A","B")] <- c(1L, 1L)
  st$born[c("A","B")] <- c(0, 0)

  arrivals <- c("X")
  role_k   <- "event"

  # Edge not incident to X
  edges <- data.frame(i = "A", j = "B", stringsAsFactors = FALSE)

  expect_error(
    log_pmf_ba_bip(st, arrivals, edges, t_k = 1, params = list(beta_edges = 0, lambda_new = 1),
                   role_k = role_k, delta = 0.1),
    "incident",
    ignore.case = TRUE
  )
})

test_that("validate_edges_ba_bip errors if a non-event endpoint is not a perp", {
  st <- state_init()
  st <- state_add_nodes(st, c("A", "E0"), t_k = 0, role = c("perp", "event"))
  st$deg["A"] <- 1L
  st$born["A"] <- 0

  arrivals <- c("X")
  role_k   <- "event"

  # tries to connect new event X to old event E0 (not a perp)
  edges <- data.frame(i = "X", j = "E0", stringsAsFactors = FALSE)

  expect_error(
    log_pmf_ba_bip(st, arrivals, edges, t_k = 1, params = list(beta_edges = 0, lambda_new = 1),
                   role_k = role_k, delta = 0.1),
    "must be perps",
    ignore.case = TRUE
  )
})
