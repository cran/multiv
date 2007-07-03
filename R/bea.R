bea <- function(a, istart=0, jstart=0)
#          Permute rows and colums, using "bond energy algorithm".
{
  if (!is.matrix(a)) stop("First input argument must be a matrix.\n")
           n <- nrow(a)
           m <- ncol(a)
           b    <- matrix(0.0, n, m)
           storage.mode(a) <- "double"
           storage.mode(b) <- "double"
           ib   <- integer(n)
           jb   <- integer(m)
           ifin <- integer(n)
           jfin <- integer(m)
           ener <- 0.0
#
           if (istart == 0) istart <- floor(runif(1,1,n))
           if (jstart == 0) jstart <- floor(runif(1,1,m))
#
           bea1  <- .Fortran("rbea",
                       n = as.integer(n),
                       m = as.integer(m),
                       a = as.matrix(a),             # input data
                       istart = as.integer(istart),  # 1st row placement
                       b = as.matrix(b),             # permuted array
                       ib = as.integer(ib),          # permuted order of rows
                       ifin = as.integer(ifin),      # for book-keeping
                             PACKAGE = "multiv")
#
           a <- bea1$b
           bea2  <- .Fortran("cbea",
                       n = as.integer(n),
                       m = as.integer(m),
                       a = as.matrix(a),             # input data
                       jstart = as.integer(jstart),  # 1st col. placement
                       b = as.matrix(b),             # permuted array
                       jb = as.integer(jb),          # permuted order of cols.
                       jfin = as.integer(jfin),      # for book-keeping
                             PACKAGE = "multiv")
#
           energ  <- .Fortran("energy",
                       n = as.integer(n),
                       m = as.integer(m),
                       b = as.matrix(bea2$b),
                       ener = as.double(ener),
                              PACKAGE = "multiv")
#
           list(b = bea2$b, ib = bea1$ib, jb = bea2$jb, e = energ$ener)
}


