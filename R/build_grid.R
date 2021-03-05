#' Build grid of values for the tempering scalar along with log-normalising constant and its derivatives.
#'
#' @param compiled.model.prior a \code{stan_model} object containing the tempered target.
#' @param eps grid increment factor. Default is \code{eps=0.01}.
#' @param M maximum value for the tempering scalar. Default is \code{M=1}.
#' @param J number of grid points. Default is \code{J=15}.
#' @param v1  scalar that controls the grid spacing in the minimum-finding phase.
#'  Default is \code{v1 = 5}.
#' @param v2 scalar that controls the grid spacing in the "gap-plugging" phase.
#' Default is \code{v2 = 10}.
#' @param stan.list list with data for the Stan model.
#' @param pars which parameters to track.
#' @param niter number of iterations, defaults to \code{5000}.
#' @param strict whether to use strict computational specs.
#' @param step_size if \code{spec=TRUE}, sets the step_size. Default is \code{step_size = 0.99}.
#' @param tdepth  if \code{spec=TRUE}, sets the tree_depth. Default is \code{tdepth = 15}.
#'
#' @return a list containing a data.frame with the grid points, log-normalising constant and derivatives.
#' If \code{pars != NA}, MCMC summaries for selected parameters are also returned.
#' @export build_grid
#' @seealso \code{\link[npowerPrioR]{build_grid_derivOnly}}
#'
build_grid <- function(compiled.model.prior, eps = .01,  M = 1, J = 15, v1 = 5, v2 = 10, stan.list, pars,
                       niter = 5000, strict = FALSE,  step_size = 0.99, tdepth = 15){
  
    f0 <- get_deriv_and_mal(a = 0.0, compiled.model.prior = compiled.model.prior, stan.list = stan.list, pars = pars,
                          strict = strict, niter = niter, tdepth = tdepth, step_size = step_size) ## need to run for a_0 = 0 to get parameter summaries.
  
  ## Preliminary checking 
  fm <- get_deriv_and_mal(a = eps, compiled.model.prior = compiled.model.prior, stan.list = stan.list, pars = pars,
                          strict = strict, niter = niter, tdepth = tdepth, step_size = step_size)
  fM <- get_deriv_and_mal(a = M, compiled.model.prior = compiled.model.prior,  stan.list = stan.list, pars = pars,
                          strict = strict, niter = niter, tdepth = tdepth, step_size = step_size)
  
  all.outs <- vector(J + 1, mode = "list")
  all.outs[[1]]$summaries <- f0$summaries
  all.outs[[2]]$summaries <- fm$summaries
  all.outs[[J + 1]]$summaries <- fM$summaries 
  
  same_sign <- sign(fm$deriv_lc) == sign(fM$deriv_lc)
  
  if(same_sign){ ## if they are the same sign, can construct a regularly spaced grid
    a0_grid <- seq(2*eps, M-eps, length.out = J - 2)
    all.outs[3:J] <- lapply(a0_grid, get_deriv_and_mal, compiled.model.prior = compiled.model.prior,
                            stan.list = stan.list, pars = pars,
                            strict = strict, niter = niter, tdepth = tdepth, step_size = step_size)
    
    inner.mals <- unlist(lapply(all.outs, function(x) x$lc ))
    inner.derivs <- unlist(lapply(all.outs, function(x) x$deriv_lc ))
    inner.2nd.derivs <- unlist(lapply(all.outs, function(x) x$second_deriv_lc ))
    
    res <- data.frame(
      a0 = c(0, eps, a0_grid, M),
      lc_a0 = c(f0$lc, ## should be zero. Maybe gives a gauge of how accurate estimation is.
                fm$lc, inner.mals, fM$lc),
      deriv_lc =  c(f0$deriv_lc, fm$deriv_lc, inner.derivs, fM$deriv_lc),
      second_deriv_lc = c(f0$second_deriv_lc, fm$second_deriv_lc, inner.2nd.derivs, fM$second_deriv_lc) 
    )
  }else{
    ## if signs change, then we need to be a bit more thoughtful
    ### Initialising
    a0s <- c(0, eps) 
    mals <- c(f0$lc, fm$lc)
    lderivs <- c(f0$deriv_lc, fm$deriv_lc)
    l2derivs <- c(f0$second_deriv_lc, fm$second_deriv_lc)
    
    ### Now start the "search"
    counter <- J - 2
    l <- eps
    u <- M
    z <- (u + l)/2
    d.old <- fm$deriv_lc
    delta <- 100 * eps
    while(counter > 0 && delta > v1*eps){
      lz <- get_deriv_and_mal(a = z, compiled.model.prior = compiled.model.prior, stan.list = stan.list,
                              pars = pars, 
                              strict = strict, niter = niter, tdepth = tdepth, step_size = step_size)
      a0s <- c(a0s, z)
      mals <- c(mals, lz$lc)
      lderivs <- c(lderivs, lz$deriv_lc)
      l2derivs <- c(l2derivs, lz$second_deriv_lc)
      all.outs[[3 + (J-2)-counter]]$summaries <- lz$summaries
      
      same_sign_below <- sign(d.old) == sign(lz$deriv_lc)
      if(same_sign_below){
        l <- z
        u <- u
      }else{ ## then it must be the case that sign(du) == sign(lz$deriv_lc)
        l <- l
        u <- z
      }
      z.old <- z
      z <- (u + l)/2
      delta <- abs(z - z.old)
      counter <- counter - 1
    }
    if(counter > 1){## in case it stopped via getting close enough to the zero
      lz <- get_deriv_and_mal(a = z, compiled.model.prior = compiled.model.prior,  stan.list = stan.list, pars = pars,
                              strict = strict, niter = niter, tdepth = tdepth, step_size = step_size)
      a0s <- c(a0s, z)
      mals <- c(mals, lz$lc)
      lderivs <- c(lderivs, lz$deriv_lc)
      l2derivs <- c(l2derivs, lz$second_deriv_lc)
      all.outs[[3 + (J-2)-counter]]$summaries <- lz$summaries
      
      counter <- counter - 1
      newrange <- c(max(0, z - v2*eps), min(z + v2 * eps, M))
      new_a0s <- c( a0s[a0s > newrange[1] & a0s < newrange[2]], newrange)
      
      while(counter > 0){
        ## above 
        z <-  plug_gap(new_a0s)
        new_a0s <- c(new_a0s, z)
        lz <- get_deriv_and_mal(a = z, compiled.model.prior = compiled.model.prior,
                                stan.list = stan.list, pars = pars,
                                strict = strict, niter = niter, tdepth = tdepth, step_size = step_size)
        lderivs <- c(lderivs, lz$deriv_lc)
        l2derivs <- c(l2derivs, lz$second_deriv_lc)
        a0s <- c(a0s, z)
        mals <- c(mals, lz$lc)
        all.outs[[3 + (J-2)-counter]]$summaries <- lz$summaries
        
        counter <- counter - 1
      }
    }
    
    a0s <- c(a0s, M)
    mals <- c(mals, fM$lc)
    lderivs <- c(lderivs, fM$deriv_lc)
    l2derivs <- c(l2derivs, fM$second_deriv_lc)
    
    res <- data.frame(a0 = a0s, lc_a0 = mals, deriv_lc = lderivs, second_deriv_lc = l2derivs)
  }
  
  if(!is.na(pars[1])) parameter.summaries <- lapply(all.outs, function(x) x$summaries)
  
  if(!is.na(pars[1])){
    out <- list(
      result = res,
      summaries = parameter.summaries
    )
  }else{
    out <- list(
      result = res
    )
  }
  return(out)
}