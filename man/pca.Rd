\name{pca}
\alias{pca}
\title{
Principal Components Analysis
}
\usage{
pca(a, method=3)
}
\description{
  Finds a new coordinate system for multivariate data such that the first
coordinate has maximal variance, the second coordinate has maximal variance
subject to being orthogonal to the first, etc.
}
\arguments{
\item{a}{
data matrix to be decomposed, the rows representing observations and the
columns variables.  Missing values are not supported.
}
\item{method}{
integer taking values between 1 and 8.  \code{Method} = 1 implies no 
transformation of data matrix.  Hence the singular value decomposition (SVD) 
is carried out on a sums of squares and 
cross-products matrix.  \code{Method} = 2 implies that the observations are centered
to zero mean.  Hence the SVD is carried out on a variance-covariance matrix.
\code{Method} = 3 (default) implies that the observations are centered to zero 
mean, and additionally reduced to unit standard deviation.  In this case the 
observations are standardized.  Hence the SVD is carried out on a correlation
matrix. \code{Method} = 4 implies that the observations are normalized by being 
range-divided, and then the variance-covariance matrix is used.  \code{Method} = 5
implies that the SVD is carried out on a Kendall (rank-order) correlation
matrix.  \code{Method} = 6 implies that the SVD is carried out on a Spearman 
(rank-order) correlation matrix.  \code{Method} = 7 implies that the SVD is carried
out on the sample covariance matrix.  \code{Method} = 8 implies that the SVD is
carried out on the sample correlation matrix.
}}
\value{
list describing the principal components analysis:

\item{rproj}{
projections of row points on the new axes.
}
\item{cproj}{
projections of column points on the new axes.
}
\item{evals}{
eigenvalues associated with the new axes. These provide figures of merit for
the \code{variance explained} by the new axes.  They are usually quoted in terms
of percentage of the total, or in terms of cumulative percentage of the total.
}
\item{evecs}{
eigenvectors associated with the new axes. This orthogonal matrix describes
the rotation.  The first column is the linear combination of columns of 
\code{a} defining the first principal component, etc.  
}}
\section{Side Effects}{
When carrying out a PCA of a hierarchy object, the partition is specified
bt \code{lev}.  The level plus the associated number of groups equals the number
of observations, at all times.
}
\section{NOTE}{
In the case of \code{method} = 3, if any column point has zero standard deviation,
then a value of 1 is substituted for the standard deviation.


Up to 7 principal axes are determined.  The inherent dimensionality of either
of the dual spaces is ordinarily \code{min(n,m)} where \code{n} and \code{m} are respectively
the numbers of rows and columns of \code{a}.  The centering transformation which is
part of \code{method}s 2 and 3 introduces a linear dependency causing the inherent
dimensionality to be \code{min(n-1,m)}.  Hence the number of columns returned in
\code{rproj}, \code{cproj}, and \code{evecs} will be the lesser of this inherent 
dimensionality and 7.


In the case of \code{methods} 1 to 4, very small negative eigenvalues, if they 
arise, are an artifact of the SVD algorithm used, and may be treated as zero. 
In the case of PCA using rank-order correlations (\code{methods} 5 and 6), negative
eigenvalues indicate that a Euclidean representation of the data is not
possible.  The approximate Euclidean representation given by the axes 
associated with the positive eigenvalues can often be quite adequate for
practical interpretation of the data.


Routine \code{prcomp} is identical, to within small numerical
precision differences, to \code{method} = 7 here.  The examples below show
how to transform the outputs of the present implementation (\code{pca}) onto 
outputs of \code{prcomp}.


Note that a very large number of columns in the input data matrix will cause
dynamic memory problems: the matrix to be diagonalized requires O(m^2) storage
m is the number of variables.
}
\section{METHOD}{
A singular value decomposition is carried out. 
}
\section{BACKGROUND}{
Principal components analysis defines the axis which provides the best fit
to both the row points and the column points.  A second axis is determined
which best fits the data subject to being orthogonal to the first.  Third and
subsequent axes are similarly found.  Best fit is in the least squares sense.
The criterion which optimizes the fit of the axes to the points is, by virtue
of Pythagoras' theorem, simultaneously a criterion which optimizes the variance
of projections on the axes.


Principal components analysis is often used as a data reduction technique.
In the pattern recognition field, it is often termed the Karhunen-Loeve
expansion since the data matrix \code{a} may be written as a series expansion using
the eigenvectors and eigenvalues found.
}
\references{
Many multivariate statistics and data analysis books include a discussion of
principal components analysis.  Below are a few examples: 


C. Chatfield and A.J. Collins, \code{Introduction to Multivariate Analysis}, 
Chapman and Hall, 1980 (a good, all-round introduction); 


M. Kendall, \code{Multivariate Analysis}, Griffin, 1980 (dated in relation to 
computing techniques, but exceptionally clear and concise in the treatment of 
many practical aspects); 


F.H.C. Marriott, \code{The Interpretation of Multiple Observations}, Academic, 
1974 (a short, very readable textbook); 


L. Lebart, A. Morineau, and K.M. Warwick, `Multivariate Descriptive 
Statistical Analysis', Wiley, 1984 (an excellent geometric treatment of PCA); 


I.T. Joliffe, \code{Principal Component Analysis}, Springer, 1980.
}
\seealso{
\code{\link{svd}}, \code{\link{prcomp}}, \code{\link{cancor}}.
}
\examples{
# principal components of the prim4 data
pcprim <- pca(prim4)
# plot of first and second principal components
plot(pcprim$rproj[,1], pcprim$rproj[,2])
# To label the points, uses `plot' with parameter `type="n"', followed by
# `text': cf. examples below.
# Place additional axes through x=0 and y=0:
plaxes(pcprim$rproj[,1], pcprim$rproj[,2])
# variance explained by the principal components
pcprim$evals*100.0/sum(pcprim$evals)
# In the implementation of the S function `prcomp', different results are
# produced.  Here is how to obtain these results, using the function `pca'.
# Consider the following result of `prcomp':
old <- prcomp(prim4)
# With `pca', one would do the following:
new <- pca(prim4, method=7)
# Data structures of `prcomp' are defined thus:
n <- nrow(prim4)
old$sdev = sqrt(new$evals/(n-1))
old$rotation = new$evec
center <- apply(old$x, 2, mean)
new$rproj[1,] <- old$x[1,] - center[1]
# One remark: the rotation matrix satisfies: 
# old$x == prim4 \%*\% old$rotation
# up to numerical precision.  However, up to 7 principal components only
# are now determined.  
#


}
\keyword{multivariate}
\keyword{algebra}
% Converted by Sd2Rd version 0.2-a5.
