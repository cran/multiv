pca <- function(a, method=3)

# PCA driver, F. Murtagh (fmurtagh@eso.org), August 1994.

# Method = 1: PCA of sums of squares and cross products;
#        = 2: PCA of covariances; 
#        = 3: PCA of correlations; default;
#        = 4: PCA of covariances of range-normalized data;
#        = 5: PCA of Kendall rank-order correlations;
#        = 6: PCA of Spearman rank-order correlations;
#        = 7: PCA of sample covariances;
#        = 8: PCA of sample correlations.

{ 

        if (!is.matrix(a)) stop("First input argument must be a matrix.\n")

        n  <- nrow(a)
        m  <- ncol(a)
        if (m > 100) 
        cat("Warning - # variables is large - possible memory problems...\n")
        b  <- matrix(0.0, m, m)
        z  <- matrix(0.0, m, m)
        storage.mode(a)  <- "double"
        storage.mode(b)  <- "double"
        storage.mode(z)  <- "double"

        if (method > 8 || method < 1) 
           stop("method must be 1, 2, 3 (default), 4, 5, 6, 7, or 8.")

        ierr   <- 0        # error indicator

        princomp <- .Fortran("pca",
             n       = as.integer(n),
             m       = as.integer(m),
             a       = as.matrix(a),
             method  = as.integer(method),
             b       = as.matrix(b),
             v1      = double(m),
             v2      = double(m),
             w1      = double(n),
             w2      = double(n),
             z       = as.matrix(z),
             ierr    = as.integer(ierr))

        if (princomp$ierr > 0) stop("No convergence for eigenvalue: ",ierr)

        inhdim  <- min(m, 7)
        inhdm0  <- min(n, m, 7)

        lim <- max(m-6,1)

        clbls <-c("Comp1","Comp2","Comp3","Comp4","Comp5","Comp6","Comp7")
        clbl  <- clbls[1:inhdim]
        vlbls <- as.list(dimnames(a)[[2]])
        olbls <- as.list(dimnames(a)[[1]])

        rproj <- princomp$a[,1:inhdim]
        if (length(olbls)!=0) dimnames(rproj) <- list(olbls,clbl)
        if (length(olbls)==0) dimnames(rproj) <- list(NULL,clbl)

        cprj  <- princomp$b[,1:inhdim]
        if (length(vlbls)!=0) dimnames(cprj) <- list(vlbls,clbl)
        if (length(vlbls)==0) dimnames(cprj) <- list(NULL,clbl)

        evcs  <- princomp$z[,m:lim]
        if (length(vlbls)!=0) dimnames(evcs) <- list(vlbls,clbl)
        if (length(vlbls)==0) dimnames(evcs) <- list(NULL,clbl)
          
        evls  <- (rev(princomp$v1))[1:inhdim]
  

        rproj <- rproj[,1:inhdm0]
        cproj <- cprj[,1:inhdm0]
        evals <- evls[1:inhdm0]
        evecs <- evcs[,1:inhdm0]
        rlbls <- NULL


      ret <- list(rproj = rproj,
                cproj = cproj,
                evals = evals,
                evecs = evecs,
                rlbls = rlbls)                      

      ret
   }
