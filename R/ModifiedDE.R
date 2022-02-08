#' Modified DE
#'
#' Modified version of the DE algorithm, proposed by Zhu et al., for solving optimization problems with bound constraints, where the mutation strategy is randomly selected
#' between \code{DE/best/1} and \code{DE/rand/1} with a probability that favours the first with the progress of the run of the algorithm.
#'
#' This modified version selects two mutation strategies every generation:
#' \itemize{
#' \item \code{DE/best/1} with probability \code{p1 = 1/(1+exp(-0.005*g))}
#' \item \code{DE/rand/1} with probability \code{1-p1}
#' }
#' Note that the probability for the strategy \code{DE/best/1} increases as the generation counter does.
#' This will favour this strategy with the run of the algorithm.
#' Mutation factor \code{F} is selected from a Gaussian distribution with mean 0.6 and standard deviation 0.2. Of course negative values
#' are not acceptable. The crossover rate is picked randomly from the range between 0.1 and 0.9.
#'
#' @param funcname function to optimize
#' @param lo lower bound (same for each component of the vector x)
#' @param up upper bound (same for each component of the vector x)
#' @param n problem size
#' @param NP population size (should be ~ four times the size of the problem)
#' @param n_gen number of generations (~ ten times the value of NP)
#' @param ... parameters to pass at objective function to optimize
#'
#' @return The output of the function \code{DEbase} is a list (of length 3) containing the following elements:\cr
#'
#' \itemize{
#' \item \code{f_best}: the best value found by the algorithm for the last generation
#' \item \code{x_best}: the vector which corresponds to the best overall function value
#' \item \code{f}: the vector of the \code{n_gen} optimal values found by the algorithm at each generation
#' }

#' @export
#'
#' @examples
#' ##implementation of the schwefel function
#'schwef <- function(x)
#'{
#'  if(is.vector(x))
#'  {
#'    d <- length(x)
#'    sum <- sum(x*sin(sqrt(abs(x))))
#'    y <- 418.9829*d - sum
#'    return(y)
#'  }
#'  if(is.matrix(x))
#'  {
#'    d <- ncol(x)
#'    sum <- apply(x*sin(sqrt(abs(x))),1,sum)
#'    y <- 418.9829*d - sum
#'    return(y)
#'  }
#'}
#'##application of DEbase function
#'set.seed(123)
#'d <- MDE(schwef,-500,500,10,40,400)
#'
#'
#'@references
#' W. Zhu, Y. Tang, J. Fang, W. Zhang, \emph{Adaptive population tuning scheme for
#' differential evolution}. Information Sciences, vol. 223 n. 13, 2013, pp. 390--401.
#'
#' S. Das, S. Mullick, P. N. Suganthan, \emph{Recent advances in differential evolution--
#' an updated survey}. Swarm and evolutionary computation, vol. 23, 2016, pp. 1--30
#'

MDE <- function(funcname,lo,up,n,NP,n_gen,...)
{

  #definizione della funzione obiettivo
  objfun <- funcname
  x_best <- matrix( ,nrow = n_gen,ncol = n)
  f_best <- vector()

  #costruzione della popolazione iniziale
  x <- matrix( ,nrow = NP,ncol = n)
  for (j in 1:NP) {
    x[j,] <- runif(n,lo,up)
  }
  #normalizzazione dei vettori
  #s<-apply(x,1,sum)
  #x=x/s
  #procedura di normalizzazione molto pi? facile
  #for (j in 1:NP) {
  #  x[j,] <- x[j,]/s[j]
  #}
  #x==x1
  f_x = objfun(x,...)
  #inizio del loop
  for (g in 1:n_gen)
  {
    #inizializzazione di F e Cr
    f = max(0,rnorm(1,0.6,0.2))
    CR = runif(1,0.1,0.9)
    #calcolo delle probabilità per le due strategie
    p1 = 1/(1+exp(-0.005*g)) #probabilit? strategia de/best/1 secondo una curva logistica
    #con k = 0.005
    if (runif(1) < p1)
    {
      strategy <- 2
    }
    else
    {
      strategy <- 1
    }
    #STRATEGIA RAND/1/BIN
    if (strategy == 1)
    {
      #inizializzazione vettore mutazione e indici per rand 1 bin
      v <- matrix( ,nrow=NP,ncol=n)
      for (i in 1:NP)
      {
        flag <- 0
        while (flag == 0)
        {
          r = sample(1:NP,1)
          if (r != i)
          {
            r1 = r
            flag <- 1
          }
        }
        flag <- 0
        while (flag == 0)
        {
          r = sample(1:NP,1)
          if (r != i && r != r1)
          {
            r2 = r
            flag <- 1
          }
        }
        flag <- 0
        while (flag == 0)
        {
          r = sample(1:NP,1)
          if (r != i && r != r2 && r != r1)
          {
            r3 = r
            flag <- 1
          }
        }
        #Operatore mutazione
        v[i,] = x[r1,] + f*(x[r2,]-x[r3,])
      }
    }
    #STRATEGIA BEST/1/BIN
    else if (strategy == 2)
    {
      r_best = which.min(f_x)
      v <- matrix( ,nrow=NP,ncol=n)
      for (i in 1:NP)
      {
        flag <- 0
        while (flag == 0)
        {
          r = sample(1:NP,1)
          if (r != i)
          {
            r2 = r
            flag <- 1
          }
        }
        flag <- 0
        while (flag == 0)
        {
          r = sample(1:NP,1)
          if (r != i && r != r2)
          {
            r3 = r
            flag <- 1
          }
        }
        #Operatore mutazione
        v[i,] = x[r_best,] + f*(x[r2,]-x[r3,])
      }
    }
    #Crossover binario
    r = matrix(runif(NP*n),nrow = NP,ncol=n)
    si = sample(1:n,NP,replace = TRUE)
    mat1 = matrix( ,nrow = NP, ncol = n)
    for (k in 1:NP)
    {
      mat1[k,] = c(1:n)
    }
    mat2 = matrix(si,NP,n)
    u = v*((r <= CR) | (mat1 == mat2)) +
      x*((r > CR) & (mat1!=mat2))
    #Procedura per la riparazione delle variabili che escono dai bordi
    L = matrix(lo,NP,n)
    U = matrix(up,NP,n)
    A = u<L
    B = u>U
    while (any(A | B))
    {
      u = A*(2*L-u) + B*(2*U-u) + (!A & !B)*u
      A = u<L
      B = u>U
    }
    #Normalizzazione dei vettori riga in u con procedura di Tarkur
    #s = apply(u,1,sum)
    #for (p in 1:NP)
    #{
    #  if (s[p]>1)
    #  {
    #    id = u[p,]>0
    #    dw = (1-sum(L[p,id]))/sum(u[p,id]-L[p,id])
    #    u[p,id] = L[p,id]+dw*(u[p,id]-L[p,id])
    #  }
    #  if (s[p]<1)
    #  {
    #    id = u[p,]>0
    #    dw = (sum(U[p,id])-1)/sum(U[p,id]-u[p,id]);
    #    u[p,id] = U[p,id]-dw*(U[p,id]-u[p,id]);
    #  }
    #}
    #Processo finale di selezione
    f_u = objfun(u,...)
    site = f_u<=f_x
    x[site,] = u[site,]
    f_x[site] = f_u[site]
    #calcolo della p measure
    #P_measure[gen] = max(norm(x-repmat(mean(x),NP,1),2))

    #Salvo alcune cose in un vettore a livello di generazione
    id_best = which.min(f_x)
    x_best[g,] = x[id_best,]
    f_best[g] = f_x[id_best]
  }

  fu_best = f_best[n_gen]
  best_x = x_best[n_gen,]
  list_de = list("x_best"=best_x,"f_best"=fu_best,"f"=f_best)
  return(list_de)
}
