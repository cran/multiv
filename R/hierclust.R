hierclust    <- function(a, method=1, bign=FALSE, diagnostics=FALSE, sim=" ",
                    movie=FALSE, option="thru",alpha=1.0,repres="chull",
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
   plot(xvalues, yvalues, type="n", xlab=xlabel,ylab=ylabel,main=titl)
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
                  circles=sigma,add=TRUE,inches=FALSE,col=0)    }    }
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
           symbols(mean(classx),mean(classy),circles=sigma,add=TRUE,inches=FALSE)  }
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
