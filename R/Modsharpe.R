#' Modified Sharpe Ratio function
#'
#' Function that calculates the modified version of the Sharpe ratio, more coherent
#' when the expected excess portfolio return is negative.
#'
#' @param x a vector or matrix - if it is a matrix, the function will return the value of the ModSharpe for each row.
#' @param mu vector of the expected returns of the n risky assets
#' @param rf risk-free rate
#' @param Sigma estimator of the covariance matrix
#'
#' @return ModSharpe() returns as object the value of the function modified sharpe ratio calculated in x. If x is a matrix,
#' it returns a column vector and calculates the modified sharpe for each row.
#' @export
#'
#' @examples
#'
#' r <- matrix( ,nrow =60 ,ncol=25)
#' for (i in 1:60) {
#'  r[i,] <- runif(25,-2.5,2.5)
#' }
#' mu <- rowMeans(t(r))
#' Sigma <- cov(r)
#' xx <- runif(25,0.01,0.30)
#' ModSharpe(xx,mu,0,Sigma)
#'
ModSharpe <- function(x,mu,rf,Sigma)
{
  num <- x%*%mu - rf
  if (is.vector(x))
  {
    den <- sqrt(x%*%Sigma%*%x)
  }
  if (is.matrix(x))
  {
    den <- sqrt(apply((x%*%Sigma)*x,1,sum))
  }
  modsharpe <- -num/((den)^sign(num))
  return(modsharpe)
}
