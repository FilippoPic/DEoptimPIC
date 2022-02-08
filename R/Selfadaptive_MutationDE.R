#' Self Adaptive mutation DE (SAmDE)
#'
#' Modified version of the base DE algorithm, proposed by Silva et al., which selects from a pool the mutation scheme
#' selecting between \code{DE/rand/1}, \code{DE/best/1}, \code{DE/best/2} and \code{DE/current-to-rand/1}.
#'
#' The probability for choosing a mutation scheme is updated by a \code{DE/rand/1} mutation operator
#' on the probability vector, starting from a vector of 4 random values.
#' Mutation factor \code{F} and crossover rate \code{CR} are randomly selected between [0.7,1.1], and
#' [0.4,0.9] respectively.
#'
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
#'d <- SAmDE(schwef,-500,500,10,40,400)
#'
#'@references
#' R. C. P. Silva, R. A. Lopes, F.G. Guimaraes, \emph{Self adaptive mutation in the differential evolution}, Proceedings of GECCO'11, 2011, Dublin, Ireland, pp. 1939--1946.
#'
#'S. Das, S. Mullick, P. N. Suganthan, \emph{Recent advances in differential evolution--
#' an updated survey}. Swarm and evolutionary computation, vol. 23, 2016, pp. 1--30
#'
SAmDE <- function(funcname,lo,up,n,NP,n_gen,...)
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
  K = 0.8
  p <- vector()
  for (k in 1:4)
  {
    p[k] = runif(1)
  }
  sum = sum(p)
  p = p/sum
  #inizio del loop
  for (g in 1:n_gen)
  {
    #inizializzazione di F e Cr
    f = runif(1,0.7,1.1)
    CR = runif(1,0.4,0.9)
    #vettore delle probabilit? per le due strategie

    for (j in 1:4)
    {
      flag <- 0
      while (flag == 0)
      {
        b = sample(1:4,1)
        if (b != j)
        {
          b1 = b
          flag <- 1
        }
      }
      flag <- 0
      while (flag == 0)
      {
        b = sample(1:4,1)
        if (b != j && b!= b1)
        {
          b2 = b
          flag <- 1
        }
      }
      p[j] = p[j] + f*(abs(p[b1]-p[b2]))
    }
    p = p/sum(p)


    tmp = runif(1)
    if (tmp <= p[1])
    {
      strategy <- 1
    }
    else if (tmp <= p[1]+p[2] & tmp > p[1])
    {
      strategy <- 2
    }
    else if (tmp <= p[1]+p[2]+p[3] & tmp > p[1]+p[2])
    {
      strategy <- 3
    }
    else
    {
      strategy <- 4
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
    #STRATEGIA BEST/2/BIN
    else if (strategy == 3)
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
    #STRATEGIA RAND-TO-CURRENT/1 che non ha il crossover
    else if (strategy == 4)
    {
      u <- matrix( ,nrow=NP,ncol=n)
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
        #Operatore mutazione
        u[i,] = x[i,] + k*(x[r1,]-x[i,]+f*(x[r1,]-x[r2,]))
      }
    }
    if (strategy != 4)
    {
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
    }
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
