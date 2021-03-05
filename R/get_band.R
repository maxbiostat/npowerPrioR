#' Get confidence bands from (GAM, lm) fit
#'
#' @param fit  an object containing a \code{lm} or \code{gam}, with a \code{predict} method.
#' @param xpred a vector with values at which to predict.
#' @param alpha the level of confidence to be used in a normal approximation (default is 95%).
#' @importFrom stats predict
#' @importFrom mgcv predict.gam
#' @return a data.frame with mean, lower and upper predictions.
#' @export get_band
#'
get_band <- function(fit, xpred, alpha = 0.95){
  z <- stats::qnorm((1+alpha)/2)
  inipred <- predict(fit, newdata = data.frame(a0 = xpred),  se.fit = TRUE)
  err <- z * inipred$se.fit
  ans <- data.frame(
    mean = inipred$fit,
    lwr =  inipred$fit - err,
    upr =  inipred$fit + err
  )
  return(ans)
}