hc <- function(a, method=1, bign=F)
{

# Hierarchical clustering, on raw input data; we will use Euclidean distance.
# A range of criteria are supported; also there is a storage-economic option.
# Author: F. Murtagh, May 1992


 if (!is.matrix(a)) {
    n <- length(a)
    m <- 1
    }
 if (is.matrix(a)) {
    n <- nrow(a)
    m <- ncol(a)
    }
 storage.mode(a) <- "double"


# Here we will branch.  We either choose the general routine, `hc', which
# caters for 7 criteria, using a half dissimilarity matrix; (BTW, this uses the
# very efficient nearest neighbor chain algorithm, which makes this algorithm
# of O(n^2) computational time, and differentiates it from the less efficient
# -- i.e. O(n^3) -- implementations in all commercial statistical packages
# -- as far as I am aware -- except Clustan.)  alternatively we branch to
# the routine `hcstoreff', which implements the Ward method again in O(n^2)
# time, but without storage of dissimilarities (-- dissimilarities are det'd.
# on the fly; the reciprocal nearest neighbor algorithm is used).


 if ( method == 1 && bign)
 {

# 1st step - get sequence of agglomerations

 istat <- 0
 hcl   <- .Fortran("hcon2",
          n = as.integer(n),
          m = as.integer(m),
          a = as.matrix(a),
          ia = integer(n),
          ib = integer(n),
          crit = double(n),
          membr = double(n),
          diss = double(n),
          ichain = integer(n),
          flag = logical(n),
          istat = as.integer(istat))

 if (hcl$istat!=0) stop("Pb. with NN-chain storage mgt. in routine hcon2\n")
 }


# Other path of branching -- more general branch...

 else
 { 

 len <- n*(n-1)/2

 hcl <- .Fortran("hc",
          n = as.integer(n),
          m = as.integer(m),
          len = as.integer(len),
          method = as.integer(method),
          a = as.matrix(a),
          ia = integer(n),
          ib = integer(n),
          crit = double(n),
          membr = integer(n),
          nn = integer(n),
          disnn = double(n),
          flag = logical(n),
          diss = double(len))
 }


# Now we're back to common ground: 
# 2nd step: interpret the information that we now have, -- seq. of aggloms., --
# as merge, height, and order lists.


 hcass <- .Fortran("hcass2",
          n = as.integer(n),
          ia = as.integer(hcl$ia),
          ib = as.integer(hcl$ib),
          order = integer(n),
          iia = integer(n),
          iib = integer(n))
 merge <- cbind(hcass$iia[1:n-1],hcass$iib[1:n-1])


 hhh <- list(merge = merge, height = hcl$crit[1:n-1], order = hcass$order)


 hhh

}



