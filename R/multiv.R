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
                       ifin = as.integer(ifin))      # for book-keeping
#
           a <- bea1$b
           bea2  <- .Fortran("cbea",
                       n = as.integer(n),
                       m = as.integer(m),
                       a = as.matrix(a),             # input data
                       jstart = as.integer(jstart),  # 1st col. placement
                       b = as.matrix(b),             # permuted array
                       jb = as.integer(jb),          # permuted order of cols.
                       jfin = as.integer(jfin))      # for book-keeping
#
           energ  <- .Fortran("energy",
                       n = as.integer(n),
                       m = as.integer(m),
                       b = as.matrix(bea2$b),
                       ener = as.double(ener))
#
           list(b = bea2$b, ib = bea1$ib, jb = bea2$jb, e = energ$ener)
}


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

flou <- function(x)
 {

   if (!is.vector(x)) stop("Only vectors handled currently as input.\n")

   log.x <- cbind(x, x)

   lo    <- quantile(x, .33)
   hi    <- quantile(x, .67)

   log.x[x >= hi, 1] <- 1   
   log.x[x >= hi, 2] <- 0   

   log.x[x <= lo, 1] <- 0
   log.x[x <= lo, 2] <- 1

   log.x[x > lo & x < hi, 1] <- (x[x > lo & x < hi] - lo) / (hi - lo)
   log.x[x > lo & x < hi, 2] <- (hi - x[x > lo & x < hi]) / (hi - lo)

   log.x

 }

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



hierclust    <- function(a, method=1, bign=F, diagnostics=F, sim=" ",
                    movie=F, option="thru",alpha=1.0,repres="chull",
                   show="change")
{


# Function to manage the calls to 'hc' and 'hclus', and to establish 
# appropriate attributes assoc. with resulting hierarchical clusterings.
# 'hclus': already-present S-Plus function; performs hier. clust. of 
#          distance or similarity input.  
# 'hc:     new function, performs hier. clust. of raw input matrix.  Eucl. 
#          dist. is used. Many criteria and options are available.  One 
#          criterion (Ward's) is supported with storage-efficient algorithm, 
#          i.e. just orig. matrix stored, and not dist. matrix.

# Author:  F. Murtagh, May 1992

#    Check whether we're dealing with input dissimilarities or an array
     if (is.matrix(a)) dat  <- "matdat"
     if (!is.matrix(a)) dat <- "dissdat"
     

     if (!movie && dat == "dissdat") {

#       Distances as input; use function 'hclust'
        
#       'hclust' assumes 'method' = 'connected', 'compact', or 'average';
#       establish appropriate default:

        if (method==1) method <- "compact"

        hier <- hclust(a, method)

#       Note: there is a bug in S+ vis-a-vis plotting of sim.-based result.
#       We need to do the following:
#           hier <- hclust(sim=a, method) 
#           hier$height <- max(hier$height) - hier$height  }
#       Here, we will ignore this.


     }

    else if (!movie)  {

#       Non-distances as input, - i.e. raw data
#       In the help file, warn the user to normalize (etc.) data first!

        if (method==1) 
           meth   <- "Ward's min. variance or error sum of squares method."
        if (method==2) 
           meth   <- "Single linkage or nearest neighbor method."
        if (method=="connected") {
           meth   <- "Single linkage or nearest neighbor method."
           method <- 2 }
        if (method==3) 
           meth   <- "Complete linkage or diameter method."
        if (method=="compact") { 
           meth   <- "Complete linkage or diameter method."
           method <- 3 }
        if (method==4) 
           meth   <- "Average linkage (group average or UPGMA) method."
        if (method=="average") { 
           meth   <- "Average linkage (group average or UPGMA) method."
           method <- 4 }
        if (method==5)  
           meth   <- "McQuitty's or WPGMA method."
        if (method==6) 
           meth   <- "Median (Gower's or WPGMC) method."
        if (method==7) 
           meth   <- "Centroid or UPGMC method."
        if (diagnostics == TRUE && method >= 1 && method <= 7) cat(meth,"\n")
        if (method >7 || method <1) 
           stop("method must be an integer between 1 and 7.") 


        hier <- hc(a, method, bign)    


    }

    else if (movie) {

#   Hierarchical clustering with real-time representation of agglomerations.


  if (!is.matrix(a)) stop("First input argument must be a matrix.\n")
  n                    <- nrow(a)
  m                    <- ncol(a)
#-----------------------------------------------------------------------------
# Input parameters: description and preliminary tests on validity.
#
#
# a is input data matrix, - no NAs allowed.
#
# method = 1 (default), 2, 3, or 4:
# 1 = minimum variance method (with the further possibility to prioritize
# linearity through use of alpha < 1.0); 2 = single linkage method; 3 =
# complete linkage method; 4 = average linkege method.
  if (!(method==1 || method==2 || method==3 || method==4))
     stop("method must be 1 (default), 2, 3, or 4.")
#
# option = "thru" (default) or "prompt":
# Go right thru all n-1 agglomerations, or ask the user whether to continue
# at each agglomeration?  Here is what we recommend: On the first occassion,
# go right through.  Check out the cluster criterion values reported on in the
# command window.  Then using the 'prompt' option, go as far as the desired
# partition, and print it out.
#
  if (!(option=="thru"||option=="prompt"))
     stop("option must be `thru` (default) or `prompt`.")
#
# alpha = 1.0 (default) or 0.0 < alpha <= 1.0:
# Only used by the minimum variance agglomerative criterion. Alpha downweights
# the variance component on the prinipal axis of the cluster, and thereby
# redefines the cluster criterion such that linear clusters are preferred.  A
# value of alpha = 1.0 gives exactly the traditional Ward's minimum variance
# method.  See Murtagh and Raftery (1984) for more on this altered criterion.
#
  if (alpha <= 0.0 || alpha > 1.0)
     stop("coefficient alpha must be > 0, <= 1. (Default 1.0).")
#
# repres = "chull" (default), "lines", "gauss":
# Representation of clusters: For input data which has dimensionality greater
# than 2, only "chull" is allowed (since the representations were found to be
# too problematic, given that the clustering takes place in m-dimensions, and
# the the representation is necessarily a 2-dimensional projection of this).
# Representation "chull" = convex hull, or straight line in case
# of colinear points.  (Note that we test for collinearity, using principal
# components analysis, since otherwise the S "hull" routine blows up.)
# Representation "lines" = lines of best fit (least squares fit is used,
# and the extremities are defined by the intersection of the min. and max. x

# values on this LS fit line). "gauss" = 1-sigma circles.  Side-effect of
# latter: warning messages indicate whenever circles overlap plot boundaries.
#
  if (!(repres=="chull"||repres=="lines"||repres=="gauss"))
     stop("repres must be one of: chull, lines, gauss.")
  if (m>2 && repres!="chull")
  warning("Only `chull' representation permitted for dim. > 2. Using this.\n")
#
# show = "change" or "all":
# In this method's present implementation, the clusters at each stage are
# 'highlighted' by means of a different symbol.  Option show = "change" has
# this highlighting changed at each agglomeration, so that only the most recent
# agglomerands are shown.  Option show = "all" does not turn off the highlight-
# ing in this way.
#
  if (!(show=="change"||show=="all"))
     stop("show option must be `change` (default) or `all`.")
#-----------------------------------------------------------------------------
#
  storage.mode(a)      <- "double"
  iklass               <- matrix(0.0, n, n) # Will store cluster assgmnts.
  storage.mode(iklass) <- "integer"
  membr                <- double(n)         # Singleton, cluster cardinalities
  flag                 <- logical(n)        # Bookkeeping: agglomerated yet?
  ia                   <- integer(n)        # Subnodes of a given node in the
  ib                   <- integer(n)        # ... hierarchy.
  crit                 <- double(n)         # Criterion val. at each level.
  potcl                <- integer(n)        # Bookkeeping: potential clus. mems
  potcl1               <- integer(n)        # Ditto
  potcl2               <- integer(n)        # Ditto
  vect1                <- double(m)
  vect2                <- double(m)
  array1               <- matrix(0.0, m, m)
  storage.mode(array1) <- "double"
  array2               <- matrix(0.0, m, m)
  storage.mode(array2) <- "double"
  centr                <- double(m)
  dat2                 <- double(n*m)
  ileast               <- 0                 # Agglomerand
  jleast               <- 0                 # Agglomerand
  dleast               <- 0.0               # Criterion value
  classx               <- 0                 # X-values corresp. to a cluster
  classy               <- 0                 # Y-values corresp. to a cluster
#  if (method==1) cat("   Alpha: ",round(as.single(alpha),5),"\n")
#
  if (method==1) meth <- "minimum variance"
  if (method==2) meth <- "single linkage"
  if (method==3) meth <- "complete linkage"
  if (method==4) meth <- "average linkage"
#
# x and y labels for plot; use principal components if dim. > 2
  if (m == 2)
     {
               xvalues              <- a[,1]
               yvalues              <- a[,2]
               xlabel               <- "coordinate 1"
               ylabel               <- "coordinate 2"
               if (method==1) {
                   titl                 <- paste("Hier. clust. ",
                                             "- crit.: ",meth,
                                             "- alpha: ",
                                              round(as.double(alpha),5) ) }
                   else       {
                   titl                 <- paste("Hier. clust. ",
                                             "- crit.: ",meth)  }
      }
  else
     {
               prco                 <- pca(a, 2)        # PCA on covariances
               xvalues              <- prco$rproj[,1]
               yvalues              <- prco$rproj[,2]
               pc1                  <- prco$evals[1]*100.0/sum(prco$evals)
               pc2                  <- prco$evals[2]*100.0/sum(prco$evals)
               xlabel               <- paste("principal component 1 (",
                                             round(as.double(pc1),2),
                                             "% of variance)")
               ylabel               <- paste("principal component 2 (",
                                             round(as.double(pc2),2),
                                             "% of variance)")
               if (method==1) {
                   titl                 <- paste("Hier. clust. ",
                                             "- crit.: ",meth,
                                             "- alpha: ",
                                              round(as.double(alpha),5),
                                              " - dim.: ",
                                              m)  }
                   else       {
                   titl                 <- paste("Hier. clust. ",
                                             "- crit.: ",meth,
                                              " - dim.: ",
                                              m)  }

     }
#
#  Set up plot
   plot(xvalues, yvalues, type="n", add=F,xlab=xlabel,ylab=ylabel,main=titl)
   points(xvalues, yvalues, type="p")
#
#  init: Initialize
   initx <- .Fortran("init",
             membr  = as.double(membr),
             flag   = as.logical(flag),
             iklass = as.matrix(iklass),
             n      = as.integer(n),
             m      = as.integer(m)
             )
#
  membr   <- initx$membr
  flag    <- initx$flag
  iklass  <- initx$iklass
#
  nmin1   <- n-1
  for (ncl in nmin1:1)
            {
#
             if (method == 1)
                {
#                gbd: Get best dissimilarity (method==1 case)
                 gbdx <- .Fortran("gbd",
                                  membr  = as.double(membr),
                                  flag   = as.logical(flag),
                                  dleast = as.double(dleast),
                                  ileast = as.integer(ileast),
                                  jleast = as.integer(jleast),
                                  ncl    = as.integer(ncl),
                                  ia     = as.integer(ia),
                                  ib     = as.integer(ib),
                                  iklass = as.matrix(iklass),
                                  potcl  = as.integer(potcl),
                                  array1 = as.matrix(array1),
                                  array2 = as.matrix(array2),
                                  vect1  = as.double(vect1),
                                  vect2  = as.double(vect2),
                                  dat2   = as.double(dat2),
                                  n      = as.integer(n),
                                  m      = as.integer(m),
                                  alpha  = as.double(alpha),
                                  a      = as.matrix(a),
                                  centr  = as.double(centr)
                                  )
                  }
             else
                 {
#                 gbd2: Get best dissimilarity (method!=1 case)
                  gbdx <- .Fortran("gbd2",
                                  membr  = as.double(membr),
                                  flag   = as.logical(flag),
                                  dleast = as.double(dleast),
                                  ileast = as.integer(ileast),
                                  jleast = as.integer(jleast),
                                  ncl    = as.integer(ncl),
                                  ia     = as.integer(ia),
                                  ib     = as.integer(ib),
                                  iklass = as.matrix(iklass),
                                  potcl1 = as.integer(potcl1),
                                  potcl2 = as.integer(potcl2),
                                  array1 = as.matrix(array1),
                                  array2 = as.matrix(array2),
                                  vect1  = as.double(vect1),
                                  vect2  = as.double(vect2),
                                  dat2   = as.double(dat2),
                                  n      = as.integer(n),
                                  m      = as.integer(m),
                                  a      = as.matrix(a),
                                  method = as.integer(method)
                                  )
                 }
#
       ileast <- gbdx$ileast
       jleast <- gbdx$jleast
       dleast <- gbdx$dleast
       membr  <- gbdx$membr
       flag   <- gbdx$flag
       ia     <- gbdx$ia
       ib     <- gbdx$ib
       iklass <- gbdx$iklass
       array1 <- gbdx$array1
       array2 <- gbdx$array2
       vect1  <- gbdx$vect1
       vect2  <- gbdx$vect2
       dat2   <- gbdx$dat2
       if (method == 1) {
           potcl  <- gbdx$potcl
           centr  <- gbdx$centr  }
       else  {
           potcl1 <- gbdx$potcl1
           potcl2 <- gbdx$potcl2 }
#
       lev    <- n-ncl
       cat("   Cluster ",lev," criterion value: ",round(dleast,5),"\n")
#
#      agg: Agglomerate
       aggx <- .Fortran("agg",
                 ileast = as.integer(ileast),
                 jleast = as.integer(jleast),
                 dleast = as.double(dleast),
                 ncl    = as.integer(ncl),
                 ia     = as.integer(ia),
                 ib     = as.integer(ib),
                 crit   = as.double(crit),
                 membr  = as.double(membr),
                 flag   = as.logical(flag),
                 n      = as.integer(n)
                 )
#
       ileast <- aggx$ileast
       jleast <- aggx$jleast
       dleast <- aggx$dleast
       ia     <- aggx$ia
       ib     <- aggx$ib
       crit   <- aggx$crit
       membr  <- aggx$membr
       flag   <- aggx$flag
#
#      gncm: Get new cluster memberships
       gncmx <- .Fortran("gncm",
                 ileast = as.integer(ileast),
                 jleast = as.integer(jleast),
                 ncl    = as.integer(ncl),
                 ia     = as.integer(ia),
                 ib     = as.integer(ib),
                 membr  = as.double(membr),
                 iklass = as.matrix(iklass),
                 flag   = as.logical(flag),
                 n      = as.integer(n)
                 )
#
       ileast <- gncmx$ileast
       jleast <- gncmx$jleast
       ia     <- gncmx$ia
       ib     <- gncmx$ib
       membr  <- gncmx$membr
       iklass <- gncmx$iklass
       flag   <- gncmx$flag
#
       if (ncl < nmin1 && show == "change")   {
#              Remove highlighting on last-agglomerated points:
               if (length(classx) > 1) {
                  lines(classx, classy, type="p", pch=4, col=0) }
#              In the case of "lines", remove them - else they cloud up things:
               if (length(classx) > 1 && repres=="lines") {
                  lines(c(xmin,xmax),c(ymin,ymax),col=0) }
#              In the case of "gauss", remove circles:
               if (length(classx) > 1 && repres=="gauss") {
                  symbols(mean(classx),mean(classy),
                  circles=sigma,add=T,inches=F,col=0)    }    }
       classvect <- iklass[n-ncl,]
       classx    <- xvalues[classvect==1]
       classy    <- yvalues[classvect==1]
       lines(classx, classy, type="p", pch=4)
#
       if (repres=="chull")  {
          if (length(classx)==2) lines(classx, classy, type="l")
          else if (length(classx)>2)  {
#            Assess linearity using PCA on covariances
             prco <- pca(cbind(classx,classy), 2)
#            Collinearity (arbitrarily!) defined as three thousandths of 2nd
#            eigenvalue (would be zero if data precisely collinear; 0.003
#            allows for small tolerance.)
             collin <- prco$evals[2]*1000.0
             if (collin <= 3)  {
                xmn <- min(classx)
                xmx <- max(classx)
                ymn <- min(classy)
                ymx <- max(classy)
                lines(c(xmn,xmx),c(ymn,ymx),type="l")    }
             else  {
                hull <- chull(classx, classy)
                polygon(classx[hull], classy[hull])      }
        }  }
        if (m==2 && repres=="gauss") {
           classn <- length(classx)
#          Use var in population terms - looks better for small 'classn':
           sigma <- sqrt(var(classx)*(classn-1)/classn+
                         var(classy)*(classn-1)/classn)
           symbols(mean(classx),mean(classy),circles=sigma,add=T,inches=F)  }
        if (m==2 && repres=="lines")  {
           z<-lsfit(classx,classy)
           xmin   <- min(classx)
           interc <- z$coef[1]
           slope  <- z$coef[2]
           ymin   <- slope*xmin+interc
           xmax   <- max(classx)
           ymax   <- slope*xmax+interc
           lines(c(xmin,xmax),c(ymin,ymax))    }
#
        if (ncl > 1 && option == "prompt") {
           cat("Proceed? [y]\n")
           resp <- scan("",what=character(),1)
        if (length(resp)==0) resp <- "y"
        if (resp=="n" || resp=="N") {
           list(lev=lev,class=iklass[1:lev,])
           return()   }
        }
#
# End of "for (ncl in nmin1:1)" loop:
   }
#
#
    part <- list(lev=lev,class=iklass[1:lev,])

# End of movie case:
  }

#   Finally return with what we came for.
#   For movie option, don't bother generating any output.

    if (movie) invisible(part)
    if (!movie) hier

}
logique <- function(x)
 {

   if (!is.vector(x)) stop("Only vectors handled currently as input.\n")

   log.x <- cbind(x, x)
   log.x[x >= median(x), 1] <- 1   
   log.x[x <  median(x), 1] <- 0
   log.x[x >= median(x), 2] <- 0   
   log.x[x <  median(x), 2] <- 1

   log.x

 }

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
          iib = integer(n))
#
  nmin1 <- n-1

# Return with cluster assignments at all levels of the hierarchy
  assn$iclass[,1:nmin1]

}


partition <- function(a, g, iter.max=15, option="mindst", diagnostics=FALSE)
{
#
#          Carry out partitioning
#          Author: F. Murtagh, May 1992

#          `g' is either: (i) number of clusters, in which case `option'
#          will be used (default option: "mindst"; alternative: "exch"; for
#          small datasets, former is better behaved - latter gives empty
#          clusters, and Inf of NaN center and compactness values);
#          or (ii) `g' is set of initial guesses for the cluster centers.


# HERE IS THE BIG DIVIDE... between current S+ routine 'kmeans', and 
# alternative routine 'exch' and 'mindst'.  


           if (is.matrix(g))
{
# HARTIGAN & WONG ALGORITHM
           partition <- kmeans(a, g, iter.max)

}
           else
{
# SPAETH ALGORITHMS (MURTAGH IMPLEMENTATION)
           ng <- g
           n <- nrow(a)
           m <- ncol(a)
           gpcen <- matrix(0.0, ng, m)
           storage.mode(gpcen) <- "double"
#          Check for valid number of groups
           if (ng <= 1) stop("ng (number of groups) must be > 1.")
           if (ng >= n) 
                    stop("ng (number of groups) must be < number of objects.")
#          Assign objects arbitrarily to groups.  Generate n uniform values,
#          and scale these up by just less than ng (= number of groups).  
#          'Floor'ing this gives an assignment in the range 0..ng-1 to each
#          object, so add 1 to this.
           memgp <- floor(runif(n)*(ng-0.00001))+1
#          Lowest acceptable number of objects assigned to a group:
           ng0 <- 1
#          Initial value of error return indicator:
           ierr <- 0
#          Now, go for it ...
           if (option == "exch") 
              partit <- .Fortran("exch",
                              a = as.double(a),      # input data
                              n = as.integer(n),
                              m = as.integer(m),
                              memgp = as.integer(memgp),  # group memberships
                              ngp0 = integer(1),     # smallest accept. cardin.
                              numgp = integer(ng),   # num. objects per group
                              gpcen = as.matrix(gpcen),# gp. cent. coords.
                              ng = as.integer(ng),   # number of gps. requested
                              comp = double(ng),     # group compactness vals.
                              ctot = double(1),      # sum of comp values
                              iter = integer(1),     # iterns. to converge
                              iter.max = as.integer(iter.max),
                              ierr = integer(1))     # error indicator
        else   
           partit <- .Fortran("mindst",
                              a = as.double(a),      # input data
                              n = as.integer(n),
                              m = as.integer(m),
                              memgp = as.integer(memgp),  # group memberships
                              ngp0 = integer(1),     # smallest accept. cardin.
                              numgp = integer(ng),   # num. objects per group
                              gpcen = as.matrix(gpcen),# gp. cen. coords.
                              ng = as.integer(ng),   # number of gps. requested
                              comp = double(ng),     # group compactness vals.
                              ctot = double(1),      # sum of comp values
                              iter = integer(1),     # iterns. to converge
                              iter.max = as.integer(iter.max),
                              ierr = integer(1))     # error indicator
           if (partit$ierr==1) stop("Invalid group number.")
#          Above error due to an object's group assignment being <1 or >ng,
#          in routine gmeans, called from mindst or from exch, which are the 
#          two principal routines stored in file part.f.  Could something
#          have gone wrong with the initial, arbitrary assignment of group
#          memberships at the beginning of this routine?  Hard to imagine...
           if (partit$ierr==2) 
            stop("A gp. has too few mems. - reduce # of gps. and try again.")


           if (partit$numgp[1] == 0) {
            cat("Warning: 1 or more classes are empty.\n")
            cat("Note that $centers values consequently contain NaN values.\n")
           } 

           partition <- partit$memgp

#          Resequence the cluster nos.
#          print(partition)
           tmp0 <- unique(partition)
           tmp  <- partition
           for (ich in 1:length(tmp0)) {
               tmp[partition==tmp0[ich]] <- ich
           }
           partition <- tmp
#          print(partition)
           if (diagnostics) {
           cat("Within sum of squares:", partit$comp,"\n")
           cat("Cluster cardinalities:", partit$numgp,"\n")
           cat("Cluster centers (displayed horizontally for each cluster):", 
              "\n")
              prmatrix(partit$gpcen)    
           }

}


   partition

}



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
plaxes <- function(a,b) {

       segments(   min(a),0,  max(a),0   )
       segments(   0,min(b),  0,max(b)   )

}

sammon <- function(a, p=2, maxit=100, tol=0.05, alpha=0.3, diagnostics=F)
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
             err   = as.double(err))

   cat("Number of iterations: ",redmap$iter,"\n")
   cat("Mapping error:        ",redmap$err,"\n")

   rproj <- redmap$b
   res   <- list(rproj = rproj)

   res

}


