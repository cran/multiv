sammon <- function(a, p=2, maxit=100, tol=0.05, alpha=0.3, diagnostics=FALSE)
{

# Sammon (non-metric multidimensional scaling) mapping.  
# Author: F. Murtagh, May 1992

   n <- nrow(a)
   m <- ncol(a)
   storage.mode(a) <- "double"
   b <- matrix(runif(n*p), n, p)
   storage.mode(b) <- "double"
   dstar <- as.vector(dist(a))
   dd    <- as.vector(dist(b))
   ndis  <- length(dstar)
   iter  <- 0
   err   <- 1.e+5
   diag  <- 0
   if (diagnostics) {
#     Note that most diagnostics are from within the Fortran program!
      cat("Mapping error\n")
      diag <- 1
   }

   redmap <- .Fortran("sammon",
                      n     = as.integer(n),
                      m     = as.integer(m),
                      p     = as.integer(p),
                      a     = as.matrix(a),
                      b     = as.matrix(b),
                      ndis  = as.integer(ndis),
                      dstar = as.double(dstar),
                      dd    = as.double(dd),
                      alpha = as.double(alpha),
                      maxit = as.integer(maxit),
                      diag  = as.integer(diag),
                      iter  = as.integer(iter),
                      tol   = as.double(tol),
                      err   = as.double(err),
                      PACKAGE = "multiv")

   cat("Number of iterations: ",redmap$iter,"\n")
   cat("Mapping error:        ",redmap$err,"\n")

   rproj <- redmap$b
   res   <- list(rproj = rproj)

   res

}


