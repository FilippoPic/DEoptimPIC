#' Schwefel function
#'
#' This function is a classic benchmark function for testing optimization algorithms.
#'
#' This function is complex with many local minima. The function presents a global minimum at
#' \bold{x'} = (420.9687,...,420.9687) and f(\bold{x'})=0.
#' Usually is evaluated on the hypercube denoted by: x_i in [-500,500].
#'
#'
#' @param x vector for calculating the schwefel function value.
#'
#' @return If \bold{x} is a vector, returns the schwefel function value corrisponding to this vector.
#' If is a matrix, returns
#' the value of the function applied on each row vector.
#'
#' @export
#'
#' @examples
#' ##INPUT:
#' xx = runif(3)
#' s <- schwefel(xx)
#'
schwefel <- function(x)
{
  if(is.vector(x))
  {
    d <- length(x)
    sum <- sum(x*sin(sqrt(abs(x))))
    y <- 418.9829*d - sum
    return(y)
  }
  if(is.matrix(x))
  {
    d <- ncol(x)
    sum <- apply(x*sin(sqrt(abs(x))),1,sum)
    y <- 418.9829*d - sum
    return(y)
  }
}
