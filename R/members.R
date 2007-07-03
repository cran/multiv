members <- function(a)
{
# Hierarchical clustering assignments, using seq. of agglomerations.
#

# Assume arg. 'a' is name of structure resulting from hier. clust.
  a <- a$merge

  if (!is.matrix(a)) stop("Input argument must be a matrix.\n")
  n  <- nrow(a)
  if (ncol(a)!=2) stop("Input array must have two columns.\n")
  nplus1 <- n+1
# Remember: # aggloms. is 1 less than # obs.
# n = # aggloms.; nplus1 = # observations.
  iclass <- matrix(0.0, nplus1, n)
  storage.mode(iclass) <- "integer"
  ia <- a[,1]
  ib <- a[,2]
#
assn <- .Fortran("assgn",
                 n = as.integer(n),
                 nplus1 = as.integer(nplus1),
                 ia = ia,
                 ib = ib,
                 iclass = iclass,
                 hvals = integer(n),
                 iia = integer(n),
                 iib = integer(n),
                 PACKAGE = "multiv")
#
  nmin1 <- n-1

# Return with cluster assignments at all levels of the hierarchy
  assn$iclass[,1:nmin1]

}


