#' Create a data.frame with log-normalising constant and its derivatives for a vector of values of the tempering factor.
#'
#' @param compiled.model.prior a \code{stan_model} object containing the tempered target.
#' @param a0_grid a vector of values of the tempering factor.
#' @param stan.list list with data for the Stan model.
#' @param pars which parameters to track.
#' @param niter number of iterations, defaults to \code{5000}.
#' @param strict whether to use strict computational specs.
#' @param step_size if \code{spec=TRUE}, sets the step_size. Default is \code{step_size = 0.99}.
#' @param tdepth  if \code{spec=TRUE}, sets the tree_depth. Default is \code{tdepth = 15}.
#'
#' @return a list containing a data.frame with the grid points, log-normalising constant and derivatives.
#' If \code{pars != NA}, MCMC summaries for selected parameters are also returned.
#' @export create_lc_df
#'
create_lc_df <- function(a0_grid, compiled.model.prior, stan.list, pars = NA,
                         niter = 5000, strict = FALSE,  step_size = 0.99, tdepth = 15){
  a0_grid <- unique(c(0, a0_grid)) # adds zero and avoids double computing if it's already there
  all.outs <- lapply(a0_grid, get_deriv_and_mal, compiled.model.prior = compiled.model.prior,
                     stan.list = stan.list, pars = pars,
                     strict = strict, niter = niter, tdepth = tdepth, step_size = step_size)
  mals <- unlist(lapply(all.outs, function(x) x$lc))
  derivs <- unlist(lapply(all.outs, function(x) x$deriv_lc))
  second.derivs <- unlist(lapply(all.outs, function(x) x$second_deriv_lc))
  if(!is.na(pars[1])) summaries <- lapply(all.outs, function(x) x$summaries)
  
  res <- data.frame(
    a0 = a0_grid,
    lc_a0 =  mals,
    deriv_lc = derivs,
    second_deriv_lc = second.derivs
  )
  if(!is.na(pars[1])){
    out <- list(
      result = res,
      summaries = summaries
    )
  }else{
    out <- list(
      result = res
    )
  }
  return(out)
}