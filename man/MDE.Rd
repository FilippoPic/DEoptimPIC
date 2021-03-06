% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModifiedDE.R
\name{MDE}
\alias{MDE}
\title{Modified DE}
\usage{
MDE(funcname, lo, up, n, NP, n_gen, ...)
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
Modified version of the DE algorithm, proposed by Zhu et al., for solving optimization problems with bound constraints, where the mutation strategy is randomly selected
between \code{DE/best/1} and \code{DE/rand/1} with a probability that favours the first with the progress of the run of the algorithm.
}
\details{
This modified version selects two mutation strategies every generation:
\itemize{
\item \code{DE/best/1} with probability \code{p1 = 1/(1+exp(-0.005*g))}
\item \code{DE/rand/1} with probability \code{1-p1}
}
Note that the probability for the strategy \code{DE/best/1} increases as the generation counter does.
This will favour this strategy with the run of the algorithm.
Mutation factor \code{F} is selected from a Gaussian distribution with mean 0.6 and standard deviation 0.2. Of course negative values
are not acceptable. The crossover rate is picked randomly from the range between 0.1 and 0.9.
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
d <- MDE(schwef,-500,500,10,40,400)


}
\references{
W. Zhu, Y. Tang, J. Fang, W. Zhang, \emph{Adaptive population tuning scheme for
differential evolution}. Information Sciences, vol. 223 n. 13, 2013, pp. 390--401.

S. Das, S. Mullick, P. N. Suganthan, \emph{Recent advances in differential evolution--
an updated survey}. Swarm and evolutionary computation, vol. 23, 2016, pp. 1--30
}
