#' Sphere function
#'
#' This function is a classic benchmark function for testing optimization algorithms.
#'
#' The function presents a global minimum at
#' \bold{x'} = (0,...,0) and f(\bold{x'})=0.
#' Usually is evaluated on the hypercube denoted by: x_i in  [-5.12,5.12] for each i.
#'
#' @param x vector for calculating the sphere function value.
#'
#' @return If \bold{x} is a vector, returns the sphere function value corrisponding to this vector.
#' If is a matrix, returns
#' the value of the function applied on each row vector.
#' @export
#'
#' @examples
#' ##INPUT:
#' xx = runif(3)
#' s <- sphere(xx)
#'
sphere <- function(x)
{
  if(is.vector(x))
  {
    sum <- sum(x^2)
    y <- sum
    return(y)
  }
  if(is.matrix(x))
  {
    s <- apply(x^2,1,sum)
    y <- s
    return(y)
  }
}
