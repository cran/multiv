\name{supplr}
\alias{supplr}
\title{
Supplementary Rows in Correspondence Analysis
}
\description{
  Using the results of a correspondence analysis, project new rows into the
  factor space.
}
\usage{
supplr(a, ca.res)
}
\arguments{
\item{a}{
data matrix to be projected.  Must have same number of columns as matrix which
was initially input to the correspondence analysis.
}
\item{ca.res}{
the output of a correspondence analysis.  The following components of this
object are used: \code{evals}, \code{rproj} and \code{cproj}.
}}
\value{
a list containing the matrix \code{proj}, 
projections of the rows of \code{a} on the correspondence analysis factors.
}
\references{
See function \code{ca}.
}
\seealso{
Correspondence analysis: \code{\link{ca}}. 
Supplementary rows and columns: \code{\link{supplr}}, \code{\link{supplc}}.  Initial data coding:
\code{\link{flou}}, \code{\link{logique}}.  Other functions producing objects of class \code{reddim}:
\code{\link{pca}}, \code{\link{sammon}}.  Other related functions: \code{\link{prcomp}}, \code{\link{cancor}}, \code{\link{cmdscale}}.
}
\examples{
data(USArrests)
USArrests <- as.matrix(USArrests)
corr <- ca(USArrests[1:25,])
newproj <- supplc(USArrests[26:50,], corr)
# plot of first and second factors, and of supplementary columns:
plot(corr$cproj[,1], corr$cproj[,2],type="n")
text(corr$cproj[,1], corr$cproj[,2])
points(newproj$proj[,1], newproj$proj[,2], col=2) 
}
\keyword{multivariate}
\keyword{algebra}
% Converted by Sd2Rd version 0.2-a5.


