supplc <- function(a, ca.res)
{

  evals <- ca.res$evals
  rproj <- ca.res$rproj
  cproj <- ca.res$cproj

  numa         <- nrow(a)
  mass         <- apply(a, 2, sum)
  newprojc     <- t(a) %*% rproj
  newprojc     <- apply(newprojc, 2, "/", mass)
  newprojc     <- apply(newprojc, 1, "/", sqrt(evals))

  list(proj=newprojc)

}

