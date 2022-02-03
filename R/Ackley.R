#' Ackley function
#'
#' This function is a classic benchmark function for testing optimization algorithms.
#'
#' Function has many local minima and represents a threat for almost any optimization algorithm.
#' The recommended values for the parameters are used, with a=20, b=0.2 and c=2*pi.
#' Usually x_i is set between [-32,32] for each i, but can also be restricted to a smaller domain.
#'
#' @param x vector for calculating the Ackley function value.
#'
#' @return If \bold{x} is a vector, returns the Ackley function value corrisponding to this vector.
#' If is a matrix, returns
#' the value of the function applied on each row vector.
#' @export
#'
#' @examples
#' ##INPUT:
#' xx = runif(3)
#' g <- ackley(xx)
#'
ackley <- function(x)
{
  a=20
  b=0.2
  c=2*pi
  ##recommended values

  if(is.vector(x))
  {
    d <- length(x)
    sum1 <- sum(x^2)
    sum2 <- sum(cos(c*x))
    term1 <- -a * exp(-b*sqrt(sum1/d))
    term2 <- -exp(sum2/d)
    y <- term1 + term2 + a + exp(1)
    return(y)
  }
  if(is.matrix(x))
  {
    d <- ncol(x)
    sum1 <- apply(x^2,1,sum)
    sum2 <- apply(cos(c*x),1,sum)
    term1 <- -a * exp(-b*sqrt(sum1/d))
    term2 <- -exp(sum2/d)
    y <- term1 + term2 + a + exp(1)
    return(y)
  }
}
