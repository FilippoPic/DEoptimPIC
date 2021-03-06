% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Selfadaptive_MutationDE.R
\name{SAmDE}
\alias{SAmDE}
\title{Self Adaptive mutation DE (SAmDE)}
\usage{
SAmDE(funcname, lo, up, n, NP, n_gen, ...)
}
\arguments{
\item{funcname}{function to optimize}

\item{lo}{lower bound (same for each component of the vector x)}

\item{up}{upper bound (same for each component of the vector x)}

\item{n}{problem size}

\item{NP}{population size (should be ~ four times the size of the problem)}

\item{n_gen}{number of generations (~ ten times the value of NP)}

\item{...}{parameters to pass at objective function to optimize}
}
\value{
The output of the function \code{DEbase} is a list (of length 3) containing the following elements:\cr

\itemize{
\item \code{f_best}: the best value found by the algorithm for the last generation
\item \code{x_best}: the vector which corresponds to the best overall function value
\item \code{f}: the vector of the \code{n_gen} optimal values found by the algorithm at each generation
}
}
\description{
Modified version of the base DE algorithm, proposed by Silva et al., which selects from a pool the mutation scheme
selecting between \code{DE/rand/1}, \code{DE/best/1}, \code{DE/best/2} and \code{DE/current-to-rand/1}.
}
\details{
The probability for choosing a mutation scheme is updated by a \code{DE/rand/1} mutation operator
on the probability vector, starting from a vector of 4 random values.
Mutation factor \code{F} and crossover rate \code{CR} are randomly selected between [0.7,1.1], and
[0.4,0.9] respectively.
}
\examples{
##implementation of the schwefel function
schwef <- function(x)
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
##application of DEbase function
set.seed(123)
d <- SAmDE(schwef,-500,500,10,40,400)

}
\references{
R. C. P. Silva, R. A. Lopes, F.G. Guimaraes, \emph{Self adaptive mutation in the differential evolution}, Proceedings of GECCO'11, 2011, Dublin, Ireland, pp. 1939--1946.

S. Das, S. Mullick, P. N. Suganthan, \emph{Recent advances in differential evolution--
an updated survey}. Swarm and evolutionary computation, vol. 23, 2016, pp. 1--30
}
