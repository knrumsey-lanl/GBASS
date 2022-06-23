#' T BASS Wrapper
#'
#' gbmars() wrapper for t likelihoods
#'
#' @param df degrees of freedom, default=1 (Cauchy)
#' @export
#' @examples
#' foo <- 1 + 1
#'
tbassold <- function(X, y, df=1, w_prior=list(type="GIG", p=1e-3, a=1e-3, b=0), ...){
  v_prior <- build_prior("GIG", c(-df/2, 0, df))
  beta_prior <- list(type="fixed", value=0)
  gbmars(X, y, w_prior, v_prior, beta_prior, ...)
}

#' bass_shoe or hbass Wrapper
#'
#' gbmars() wrapper for horseshoe likelihood
#'
#' @param prop_sigma proposal sd for v_i
#' @export
#' @examples
#' foo <- 1 + 1
#'
hbassold <- function(X, y, prop_sigma, w_prior=list(type="GBP", p=1, a=1/2, b=1/2, prop_sigma=0.5), control=list(), ...){
  con <- list(iter=5000, Jmax=2, Mmax=300, Pmove=c(1/3, 1/3, 1/3),
              b0=10, tau_prior=c(1/2, length(y)/2), lam_prior=c(1, 1),
              init_weights=c(1,1), c=1)
  con[(namec <- names(control))] <- control #Swap any user inputs
  con$c <- var(y)

  v_prior <- list(type="GBP", p=1, a=1/2, b=1/2, prop_sigma=prop_sigma)
  beta_prior <- list(type="fixed", value=0)
  gbmars(X, y, w_prior, v_prior, beta_prior, control=con, ...)
}


#' quantile regression bass wrapper
#'
#' gbmars() wrapper for quantile regression
#'
#' @param q quantile of interest (default 0.5)
#' @export
#' @examples
#' foo <- 1 + 1
#'
qbassold <- function(X, y, q=0.5, w_prior=list(type="GIG", p=1e-3, a=1e-3, b=0, prop_sigma=0.2), control=list(), ...){
  con <- list(iter=5000, Jmax=2, Mmax=300, Pmove=c(1/3, 1/3, 1/3),
              b0=10, tau_prior=c(1/2, length(y)/2), lam_prior=c(1, 1),
              init_weights=c(1,1), c=1)
  con[(namec <- names(control))] <- control #Swap any user inputs
  con$c <- 2/(q*(1-q))

  v_prior <- build_prior("GIG", c(1, 2, 0))
  beta_prior <- list(type="fixed", value=1/q - 1/(1-q))
  gbmars(X, y, w_prior, v_prior, beta_prior, control=con, ...)
}





