#' Random generator for the "monomial perturbation of normal" distribution
#'
#' A function to generate from the mpon distribution. Uses a uniformly bounded rejection sampling scheme, adapted from Devroye (2014).
#'
#' @param n currently, n = 1. Changing this argument has no effect on the code.
#' @param alpha parameter
#' @param gamma parameter
#' @param mu parameter
#' @param max_steps maximum number of steps for Newton Raphson
#' @param tol tolerance for Newton Raphson
#' @return a random draw from the mpon distribution
#' @details Generates a RV with (unnormalized) density function f(x) = x^(alpha-1)exp(-gamma*(x-mu)^2)
#' @export
#' @examples
#' n <- 10000
#' alpha <- 5
#' gamma <- 1
#' mu <- 2
#' y <- rep(NA, n)
#' for(i in 1:n){
#'   y[i] <- rmpon(1, alpha, gamma, mu)
#' }
#' hist(y, breaks=30, freq=F)
#' c <- integrate(dmpon, lower=0, upper=20, alpha=alpha, gamma=gamma, mu=mu)$value
#' curve(dmpon(x, alpha, gamma, mu)/c, add=TRUE, lwd=2, col='blue')
rmpon <- function(n=1, alpha, gamma, mu, max_steps=50, tol=1e-4){
  m <- mu/2 + 1/(2*gamma)*sqrt(gamma*(2*alpha+gamma*mu^2-2))
  lf <- function(xx) dmpon(xx, alpha, gamma, mu, log=TRUE)
  phi <- function(xx) exp((alpha-1)*(log(xx+m)-log(m)) - gamma*((xx+m-mu)^2 - (m-mu)^2))
  psi <- function(xx) (alpha-1)*log(xx+m) - gamma*(xx+m-mu)^2 - lf(m)
  psip <- function(xx)  (alpha-1)/(xx+m) - 2*gamma*(xx+m-mu)
  #Newton-Raphson
  s <- t <- 0.5
  #Get decent starting values
  if(-psi(t) < 1){
    while(-psi(t) < 0.1){
      t <- 3*t
      s <- 3*s
    }
  }else{
    while(-psi(t) > 10){
      t <- t/3
      s <- s/3
    }
  }
  t_flag <- s_flag <- TRUE
  cnt <- 1
  while((t_flag | s_flag) & cnt <= max_steps){
    t <- t - (psi(t)+1)/psip(t)
    s <- s + (psi(-s)+1)/psip(-s)

    if(abs(psi(t) + 1) < tol) t_flag <- FALSE
    if(abs(psi(-s) + 1) < tol) s_flag <- FALSE
    cnt <- cnt + 1
  }
  if(cnt > max_steps) warning("cnt reached max_steps -- may not have converged. psi(t) = ", round(psi(t), 7), " and psi(-s) = ", round(psi(-s),7))

  if(min(s,t) < 0) stop("Error: s and/or t is negative")
  #Devroye algorithm
  p <- 1/psip(-s)
  r <- -1/psip(t)
  tp <- t + r*psi(t)
  sp <- s + p*psi(-s)
  q <- tp + sp
  flag <- TRUE

  while(flag){
    U <- runif(1); V <- runif(1); W <- runif(1)
    if(U < q/(q+r+p)){
      X <- -sp + q*V
    }else{
      if(U < (q+r)/(q+r+p)){
        X <- tp - r*log(V)
      }else{
        X <- -sp + p*log(V)
      }
    }

    if(X > -m){
      if(X > tp){
        chi <- exp(psi(t) + psip(t)*(X-t))
      }else{
        if(X > -sp){
          chi <- 1
        }else{
          chi <- exp(psi(-s) + psip(-s)*(X+s))
        }
      }
      if(log(W) <= psi(X) - log(chi)){
        flag <- FALSE
      }
    }
  }
  return(X+m)
}


#' Un-normalized density of the "monomial perturbation of normal (mpon)" distribution
#'
#' A function to evaluate the unnormalized density of the mpon distribution.
#'
#' @param x locations for density evaluation
#' @param alpha parameter
#' @param gamma parameter
#' @param mu parameter
#' @return density of the mpon distribution
#' @details Evaluates the (unnormalized) density function f(x) = x^(alpha-1)exp(-gamma*(x-mu)^2)
#' @export
#' @examples
#' n <- 10000
#' alpha <- 5
#' gamma <- 1
#' mu <- 2
#' y <- rep(NA, n)
#' for(i in 1:n){
#'   y[i] <- rmpon(1, alpha, gamma, mu)
#' }
#' hist(y, breaks=30, freq=F)
#' c <- integrate(dmpon, lower=0, upper=20, alpha=alpha, gamma=gamma, mu=mu)$value
#' curve(dmpon(x, alpha, gamma, mu)/c, add=TRUE, lwd=2, col='blue')
dmpon <- function(x, alpha, gamma, mu, log=FALSE){
  res <- (alpha-1)*log(x) - gamma*(x-mu)^2 + log(x > 0)
  if(log){
    return(res)
  }else{
    return(exp(res))
  }
}

