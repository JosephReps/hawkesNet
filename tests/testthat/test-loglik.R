test_that("loglik() matches hand-computed temporal Hawkes likelihood when mark pmf is zero", {
  # This is essentially just confirming the temporal part of the likelihood is
  # working fine

  # Simple event times: 1, 2, 4 (no node/edge content needed for this test)
  nodes <- data.frame(id = c("a", "b", "c"), time = c(1, 2, 4))
  edges <- data.frame(i = character(0), j = character(0), time = numeric(0))
  ev <- make_events(nodes = nodes, edges = edges)

  params <- list(mu = 10, K = 3, beta = 2)

  # mark pmf contributes 0 always; state doesn't change
  mark0 <- function(state, arrivals_k, edges_k, t_k, delta) 0
  step0 <- function(state, arrivals_k, edges_k, t_k) state

  val <- loglik(
    ev = ev,
    params = params,
    T0 = 0,
    T_end = 4,
    mark_logpmf = mark0,
    state0 = state_init(),
    state_step_fun = step0
  )

  # Hand calculation:
  # event times = c(1, 2, 4)
  t <- c(1, 2, 4)

  mu <- params$mu
  K <- params$K
  beta <- params$beta

  # recursion: S_pre(t1)=0; after event, S_post = S_pre + 1
  # Between events, S decays by exp(-beta*dt)
  S_pre_1 <- 0
  lambda1 <- mu + K * S_pre_1

  S_post_1 <- S_pre_1 + 1
  S_pre_2 <- exp(-beta * (t[2] - t[1])) * S_post_1
  lambda2 <- mu + K * S_pre_2

  S_post_2 <- S_pre_2 + 1
  S_pre_3 <- exp(-beta * (t[3] - t[2])) * S_post_2
  lambda3 <- mu + K * S_pre_3

  ll_sum <- log(lambda1) + log(lambda2) + log(lambda3)

  # compensator:
  # on (0,1): S_post_prev=0
  # on (1,2): S_post_prev = S_post_1
  # on (2,4): S_post_prev = S_post_2
  comp <- 0
  comp <- comp + (mu * (t[1] - 0) + (K / beta) * 0 * (1 - exp(-beta * (t[1] - 0))))
  comp <- comp + (mu * (t[2] - t[1]) + (K / beta) * S_post_1 * (1 - exp(-beta * (t[2] - t[1]))))
  comp <- comp + (mu * (t[3] - t[2]) + (K / beta) * S_post_2 * (1 - exp(-beta * (t[3] - t[2]))))

  expected <- ll_sum - comp

  expect_equal(val, expected, tolerance = 1e-12)
})

test_that("loglik() errors on decreasing times (in the case ev is malformed)", {
  # Construct a deliberately malformed ev (bypassing make_events() sorting/validation)
  ev_bad <- list(
    times = data.frame(event_id = 1:2, t = c(2, 1)),
    nodes = data.frame(event_id = integer(0), id = character(0)),
    edges = data.frame(event_id = integer(0), i = character(0), j = character(0))
  )
  class(ev_bad) <- c("hg_events", "list")

  params <- list(mu = 10, K = 3, beta = 2)
  mark0 <- function(state, arrivals_k, edges_k, t_k, delta) 0
  step0 <- function(state, arrivals_k, edges_k, t_k) state

  expect_error(
    loglik(ev_bad, params = params, T0 = 0, T_end = 2, mark_logpmf = mark0, state0 = state_init(), state_step_fun = step0),
    "must be strictly increasing"
  )
})

test_that("loglik() BA returns correct likelihood value (Validated on old package)", {
  # Nodes 5 through 10 arrive at corresponding times
  nodes <- data.frame(
    id   = c(5, 6, 7, 8, 9, 10),
    time = c(1, 2, 3, 4, 7, 10)
  )

  edges <- data.frame(
    i    = c(4, 4, 5, 6, 4, 9),
    j    = c(5, 6, 7, 8, 9, 10),
    time = c(1, 2, 3, 4, 7, 10)
  )

  params <- list(mu = 10,
                 beta = 0.1,
                 K = 0.1,
                 beta_edges = 0.1)

  # Initial network of a single node (id = 4)
  state0 <- make_state(nodes0 = data.frame(id = 4), T0 = 0)

  ev <- make_events(nodes = nodes, edges = edges, state0 = state0)

  # debug = TRUE prints per-event contributions
  ll_1 <- loglik(ev, params = params, state0 = state0, T0 = 0, delta = 0.5)

  expect_equal(ll_1, -99.7706689122, tolerance = 1e-8)
})

# Need to add validation for this in the likelihood function
# test_that("loglik() errors if T_end < T0", {
#   nodes <- data.frame(id = c("a"), time = c(1))
#   ev <- make_events(nodes = nodes, edges = data.frame(i=character(), j=character(), time=numeric()))
#   params <- list(mu=1, K=1, beta=1)
#   mark0 <- function(...) 0
#   step0 <- function(s, ...) s
#
#   expect_error(
#     loglik(ev, params=params, T0=2, T_end=1, mark_logpmf=mark0, state0=state_init(), state_step_fun=step0),
#     "T_end"
#   )
# })

test_that("loglik() changes with mu", {
  nodes <- data.frame(id=c("a","b"), time=c(1,2))
  ev <- make_events(nodes = nodes,
                    edges = data.frame(i=character(), j=character(), time=numeric()))
  mark0 <- function(...) 0
  step0 <- function(s, ...) s

  ll1 <- loglik(ev, params=list(mu=1, K=0, beta=2), T0=0, T_end=2,
                mark_logpmf=mark0, state0=state_init(), state_step_fun=step0)
  ll2 <- loglik(ev, params=list(mu=2, K=0, beta=2), T0=0, T_end=2,
                mark_logpmf=mark0, state0=state_init(), state_step_fun=step0)

  expect_true(is.finite(ll1))
  expect_true(is.finite(ll2))
  expect_false(isTRUE(all.equal(ll1, ll2)))
})

# Also want to add improved debugging, and returning each component of the
# likelihood per-event, with corresponding tests
