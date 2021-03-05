#' Mean and Bayesian credibility interval (BCI)
#'
#' @param x a vector of MCMC samples.
#' @param alpha a credibility level. Default is 95%. 
#'
#' @return a data frame containing the mean and alpha*100% BCI for x.
#'
mean_bci <- function(x, alpha = 0.95){
  quants <- stats::quantile(x, probs = c( (1-alpha)/2 , (1 + alpha)/2)) 
  out <- data.frame(
    lwr = quants[1],
    mean = mean(x, na.rm = TRUE),
    upr = quants[2]
  )
  return(out)
}
#' Obtain all summaries for stanfit
#'
#' @param obj a stanfit object.
#'
#' @return a data.frame with mean and BCI for each 'parameter' in \code{obj}
#' @export summarise_posterior
#'
summarise_posterior <- function(obj){
  K <- length(obj)
  obj <- lapply(obj, as.matrix) ## taking care of one-dimensional "arrays"
  dims <- lapply(obj, ncol)
  par.names <- names(obj)
  all.names <-  vector(K, mode = "list")
  for(k in 1:K){
    all.names[[k]] <- paste(par.names[k], '[', 1:dims[[k]], ']', sep = "")
  }
  all.names.flat <- unlist(all.names)
  summys <- lapply(obj, function(par.dt) apply(par.dt, 2, mean_bci))
  res <- data.frame(
    parameter = all.names.flat,
    do.call(rbind, lapply(summys, function(y) do.call(rbind, y)))
  )
  row.names(res) <- NULL
  return(res)
}
#' Plug the gaps in a vector of a0s.
#'
#' @param x a vector of a0s.
#'
#' @return a new vector a0s with "gaps" filled in.
#'
plug_gap <- function(x){
  sx <- sort(x)
  deltas <- diff(sx)
  pos <- which.max(deltas)
  new.x <- mean(sx[c(pos, pos + 1)])
  return(new.x)
}