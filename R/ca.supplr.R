supplr <- function(a, ca.res)
{

  evals <- ca.res$evals
  rproj <- ca.res$rproj
  cproj <- ca.res$cproj

  numa         <- nrow(a)
  mass         <- apply(a, 1, sum)
  newprojr     <- a %*% cproj
  newprojr     <- apply(newprojr, 2, "/", mass)
  newprojr     <- apply(newprojr, 1, "/", sqrt(evals))

  list(proj=newprojr)

}

