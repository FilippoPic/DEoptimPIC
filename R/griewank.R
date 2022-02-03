#' Griewank function
#'
#' This function is a classic benchmark function for testing optimization algorithms.
#'
#' Complex function with many local minima, which are regularly distributed. Usually evaluated
#' with x_i in [-600,600] for each i. Global minimum f(\bold{x'})=0 in \bold{x'}=(0,...,0).
#'
#' @param x vector for calculating the Griewank function value.
#'
#' @return If \bold{x} is a vector, returns the Griewank function value corrisponding to this vector.
#' If is a matrix, returns
#' the value of the function applied on each row vector.
#' @export
#'
#' @examples
#'
#' ##INPUT:
#' xx = runif(3)
#' g <- griewank(xx)
#'
griewank <- function(x){
  if(is.vector(x))
  {
  ii <- c(1:length(x))
  sum <- sum(x^2/4000)
  prod <- prod(cos(x/sqrt(ii)))
  y <- sum - prod + 1
  return (y)
  }
  if(is.matrix(x))
  {
  i <- c(1:ncol(x))
  ii <- rbind(i)[rep(1,nrow(x)),]
  sum <- apply(x^2/4000,1,sum)
  prod <- apply(cos(x/sqrt(ii)),1,prod)
  y <- sum - prod + 1
  return (y)
  }
}
