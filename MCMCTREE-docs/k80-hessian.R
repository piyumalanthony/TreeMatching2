# K80 log-likelihood
# d: distance, k: transtion (ts) / transversion (tv) ratio, ns: number of
# transition differences nv: number of transversion differences, 
# n: alignment length
k80.2slnL <- function(d, k, ns, nv, n) {
  
  p0 <- .25 + .25 * exp(-4 * d / (k + 2)) +
    .5 * exp(-2 * d * (k + 1) / (k + 2))
  
  p1 <- .25 + .25 * exp(-4 * d / (k + 2)) -
    .5 * exp(-2 * d * (k + 1) / (k + 2))
  
  p2 <- .25 - .25 * exp(-4 * d / (k + 2))
  
  return ((n - ns - nv) * log(p0 / 4) +
            ns * log(p1 / 4) + nv * log(p2 / 4))
}

# wrapper of K80 likelihood for numerical optimisation
k80.optim.wrap <- function(par, ...) {
  -k80.2slnL(d=par[1], k=par[2], ...)
}

# analytic MLE's
k80.2smle <- function(ns, nv, n) {
  dk <- numeric(2)
  S <- ns / n; V <- nv / n
  dk[1] = -.5 * log (1 - 2 * S - V) - .25 * log(1 - 2 * V)
  dk[2] = 2 * log(1 - 2 * S - V) / log(1 - 2 * V) - 1
  names(dk) <- c("d", "k")
  return(dk)
}

D0 <- expression( log( (.25 + .25 * exp(-4 * d / (k + 2)) +
                          .5 * exp(-2 * d * (k + 1) / (k + 2))) / 4 ))
D1 <- expression( log( (.25 + .25 * exp(-4 * d / (k + 2)) -
                           .5 * exp(-2 * d * (k + 1) / (k + 2))) / 4) )
D2 <- expression( log( (.25 - .25 * exp(-4 * d / (k + 2))) / 4 ))

# first derivatives of log-likelihood at d:
D0.d <- D(D0, "d"); D1.d <- D(D1, "d"); D2.d <- D(D2, "d")
# first derivatives of log-likelihood at k:
D0.k <- D(D0, "k"); D1.k <- D(D1, "k"); D2.k <- D(D2, "k")

# returns the gradient matrix (calculated analytically) of 
# the K80 log-likelihood on site patterns and on d and k
# n0 = n - ns - nv, ns, nv
k80.gradient <- function(d, k, ns, nv, n) {
  g <- matrix(nrow=2, ncol=3)
  
  g[1,1] <- eval(D0.d) * (n - ns - nv)
  g[1,2] <- eval(D1.d) * ns
  g[1,3] <- eval(D2.d) * nv
  
  g[2,1] <- eval(D0.k) * (n - ns - nv)
  g[2,2] <- eval(D1.k) * ns
  g[2,3] <- eval(D2.k) * nv
  
  rownames(g) <- c("d", "k")
  colnames(g) <- c("conserved", "ts", "tv")
  
  return (g)
}