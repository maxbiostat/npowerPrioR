#'  Compute the derivative of the log-normalising constant  for a given tempering scalar.
#'
#' @param compiled.model.prior a \code{stan_model} object containing the tempered target.
#' @param a  the tempering scalar.
#' @param stan.list list with data for the Stan model.
#' @param pars which parameters to track.
#' @param niter number of iterations, defaults to \code{5000}.
#' @param strict whether to use strict computational specs.
#' @param step_size if \code{spec=TRUE}, sets the step_size. Default is \code{step_size = 0.99}.
#' @param tdepth  if \code{spec=TRUE}, sets the tree_depth. Default is \code{tdepth = 15}.
#' 
#' @return a list with the \code{stanfit} object, the first and second derivatives of the log-normalising constant.
#' @export get_deriv_only
#' @seealso \code{\link[npowerPrioR]{get_deriv_and_mal}}
#' 
get_deriv_only <- function(a, compiled.model.prior, stan.list, pars = NA, niter = 5000, strict = FALSE,  step_size = 0.99, tdepth = 15){
  stan.list$a_0 <- a
  if(strict){
    power.prior <- rstan::sampling(object = compiled.model.prior, data = stan.list, iter = 5000,
                                   control = list(adapt_delta = 0.99, max_treedepth = 15), refresh = 0)
  }else{
    power.prior <- rstan::sampling(object = compiled.model.prior, data = stan.list, refresh = 0)  
  }
  second_deriv_part <- mean(rstan::extract(power.prior, 'logL_sq')$logL_sq)
  first_deriv <- mean(rstan::extract(power.prior, 'logL')$logL)
  result <- list(
    chain = power.prior,
    deriv_lc = first_deriv,
    second_deriv_lc = second_deriv_part - (first_deriv^2)
  )
  if(!is.na(pars[1])) result$summaries <- summarise_posterior(rstan::extract(power.prior, pars = pars))
  return(result)
}