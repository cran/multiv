ca <- function(a, nf=7)
{
 if (!is.matrix(a)) stop("Input argument must be a matrix.\n")
          n   <- nrow(a)
          m   <- ncol(a)
          b   <- matrix(0.0, m, m)
          z   <- matrix(0.0, m, m)
          rcn <- matrix(0.0, n, m)
          ccn <- matrix(0.0, m, m)
          storage.mode(a)   <- "double"
          storage.mode(b)   <- "double"
          storage.mode(z)   <- "double"
          storage.mode(rcn) <- "double"
          storage.mode(ccn) <- "double"
#
          ierr   <- 0        # error indicator
#         Up to 'nf' prin. axes are determined in this implementation.  
          nf <- min(n, m, nf)
#
          corresp <- .Fortran("ca",
               n      = as.integer(n),
               m      = as.integer(m),
               a      = as.matrix(a),
               b      = as.matrix(b),
               v1     = double(m),
               v2     = double(m),
               z      = as.matrix(z),
               rmass  = double(n),
               cmass  = double(m),
               rcn    = as.matrix(rcn),
               ccn    = as.matrix(ccn),
               nf     = as.integer(nf),
               ierr   = as.integer(ierr))
    if (corresp$ierr==2) 
       stop("All vals. in a row or col. are zero; pls. remove.")
    if (corresp$ierr > 0) stop("No convergence for eigenvalue: ",ierr)

    clbls  <- as.list(paste("Factor",1:m,sep=""))

    dimnames(corresp$rcn) <- list(NULL, clbls)
    dimnames(corresp$ccn) <- list(NULL, clbls)
    dimnames(corresp$a)   <- list(NULL, clbls)
    dimnames(corresp$b)   <- list(NULL, clbls)


    res <- list(evals = (rev(corresp$v1))[2:(nf+1)],
                rproj = corresp$a[,1:nf], 
                cproj = corresp$b[,1:nf],
                rcntr = corresp$rcn[,1:nf],
                ccntr = corresp$ccn[,1:nf]
               )

    res
}




