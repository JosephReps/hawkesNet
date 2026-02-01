#' Title
#'
#' @param events
#' @param params_init
#' @param fixed_params
#' @param transform
#' @param hessian
#' @param control
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
fit_hawkesNet <- function(events,
                          params_init,
                          fixed_params = NULL,
                          transform = NULL,
                          hessian = FALSE,
                          control = list(),
                          ...
                          ) {
  # 0) Input validation

  # Make sure starting values are provided for the optimizer
  stopifnot(is.list(params_init), length(params_init) > 0)

  start_time <- proc.time()

  # 1) Deal with the parameter transformations

  # If no transformation information was provided, assume we are fitting on the
  # log scale
  if (is.null(transform)) {
    transform <- rep("log", length(params_init))
    names(transform) <- names(params_init)
  }
  # Make sure that each parameter has a valid transform, default to "log"
  # and warn user in the case any are missing
  missing_tf <- setdiff(names(params_init), names(transform))
  if (length(missing_tf) > 0) {
    warning(paste0('Missing `transform` entries for: ',
                  paste(missing_tf, collapse = ", "), '. Defaulting to "log"'))

    transform[missing_tf] <- "log"
  }

  # 2) Deal with any possible fixed parameters

  # Coerce fixed_params to character(0) if NULL
  fixed_params <- fixed_params %||% character(0)
  # Make sure to warn the user if fixed_params contains names not present in
  # params_init
  unknown_fixed <- setdiff(fixed_params, names(params_init))
  if (length(unknown_fixed) > 0) {
    warning(
      "The following fixed_params are not present in params_init and will be ignored: ",
      paste(unknown_fixed, collapse = ", "))
  }
  # Derive which parameters are free using fixed_params
  fixed_params <- intersect(fixed_params, names(params_init))
  free_params <- setdiff(names(params_init), fixed_params)

  # And in the case there are no free parameters, just return the likelihood
  # at the param_init values
  if (length(free_params) == 0L) {
    ll <- loglik(ev = events, params = params_init, ...)
    return(list(
      fit = NULL,
      par = params_init,
      loglik = ll,
      convergence = 0L,
      message = "No free parameters (all fixed)."
    ))
  }

  # 4) Optimizer setup

  # Define the objective function
  fn <- function(x) {
    # Make sure to transform to back from working scale
    p <- from_working(x, transform)
    ll <- tryCatch(
      loglik(ev = events, params = p, ...),
      error = function(e) -Inf
    )
    if (!is.finite(ll)) return(1e100)
    -ll
  }

  # Transform the initial parameter values to working scale
  x0_all <- to_working(params_init, transform)
  # Make sure to only optimize the free parameters
  x0 <- x0_all[free_params]

  # 5) Fit the model
  fit <- stats::optim(
    par = x0,
    fn = fn,
    method = "BFGS",
    control = control,
    hessian = hessian
  )

  # 6) Structure output
  x_all <- x0_all
  x_all[free_params] <- fit$par

  p_hat <- from_working(x_all, transform)
  ll_hat <- loglik(ev = events, params = p_hat, ...)

  out <- list(
    fit = fit,
    par = p_hat,
    loglik = ll_hat,
    convergence = fit$convergence,
    message = fit$message
  )
  if (hessian) out$hessian <- fit$hessian

  print(paste0("Fitting took ", round((proc.time()-start_time)[3],2), " seconds"))
  out
}

