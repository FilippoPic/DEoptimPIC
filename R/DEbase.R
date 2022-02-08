#' Base DE
#'
#' Differential Evolution procedure to solve minimization problems with bound constraints
#'
#' @param funcname function to optimize
#' @param lo lower bound (same for each component of the vector x)
#' @param up upper bound (same for each component of the vector x)
#' @param n problem size
#' @param NP population size (should be ~ four times the size of the problem)
#' @param n_gen number of generations (~ ten times the value of NP)
#' @param f F value (usually from 0.1 to 1.1)
#' @param CR crossover rate (it is a probability, should stay between 0 and 1)
#' @param strategy type of strategy implemented. Default to DE/RAND/1
#' @param ... parameters to pass at objective function to optimize
#'
#' @return The output of the function \code{DEbase} is a list (of length 3) containing the following elements:\cr
#'
#' \itemize{
#' \item \code{f_best}: the best value found by the algorithm for the last generation
#' \item \code{x_best}: the vector which corresponds to the best overall function value
#' \item \code{f}: the vector of the \code{n_gen} optimal values found by the algorithm at each generation
#' }
#'
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
#'d <- DEbase(schwef,-500,500,10,40,400,0.8,0.4,strategy=2)
#'
#'@references
#' S. Das, S. Mullick, P. N. Suganthan, \emph{Recent advances in differential evolution--
#' an updated survey}. Swarm and evolutionary computation, vol. 23, 2016, pp. 1--30
#'

DEbase <- function(funcname,lo,up,n,NP,n_gen,f,CR,strategy,...)
{

  #objective function definition
  objfun <- funcname
  x_best <- matrix( ,nrow = n_gen,ncol = n)
  f_best <- vector()

  #population initialization
  x <- matrix( ,nrow = NP,ncol = n)
  for (j in 1:NP) {
  x[j,] <- runif(n,lo,up)
  }
  f_x = objfun(x,...)
  #loop start
  for (g in 1:n_gen)
  {
    #RAND/1/BIN strategy

    if (strategy == 1)
    {
      #initialization of the mutation vector
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
    #STRATEGIA RAND/2/BIN
    else if (strategy == 3)
    {
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
        flag <- 0
        while (flag == 0)
        {
          r = sample(1:NP,1)
          if (r != i && r != r2 && r != r1 && r!=r3)
          {
            r4 = r
            flag <- 1
          }
        }
        flag <- 0
        while (flag == 0)
        {
          r = sample(1:NP,1)
          if (r != i && r != r2 && r != r1 && r!=r3 && r!=r4)
          {
            r5 = r
            flag <- 1
          }
        }
        #Operatore mutazione
        v[i,] = x[r1,] + f*(x[r2,]-x[r3,]) + f*(x[r4,]-x[r5,])
      }
    }
    #STRATEGIA BEST/2/BIN
    else if (strategy == 4)
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
        flag <- 0
        while (flag == 0)
        {
          r = sample(1:NP,1)
          if (r != i && r != r2 && r != r3)
          {
            r4 = r
            flag <- 1
          }
        }
        flag <- 0
        while (flag == 0)
        {
          r = sample(1:NP,1)
          if (r != i && r != r2 && r != r4 && r != r3)
          {
            r5 = r
            flag <- 1
          }
        }
        #Operatore mutazione
        v[i,] = x[r_best,] + f*(x[r2,]-x[r3,])+ f*(x[r4,]-x[r5,])
      }
    }
    #STRATEGIA RAND-TO-BEST/1/BIN
    else if (strategy == 5)
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
          if (r != i && r != r1 && r != r2)
          {
            r3 = r
            flag <- 1
          }
        }
        flag <- 0
        while (flag == 0)
        {
          r = sample(1:NP,1)
          if (r != i && r != r1 && r != r2 && r != r3)
          {
            r4 = r
            flag <- 1
          }
        }
        #Operatore mutazione
        v[i,] = x[r1,] + f*(x[r_best,]-x[r2,])+ f*(x[r3,]-x[r4,])
      }
    }
    #Crossover binario
    r = matrix(runif(NP*n),nrow = NP,ncol=n)
    si = sample(1:n,NP)
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
