\name{members}
\alias{members}
\title{
Cluster Memberships 
}
\usage{
members(a)
}
\description{
  From a succession of agglomerations produced by a hierarchical routine
in function `hierclust' or `hclust', determine cluster assignments of all 
objects, at all levels of the hierarchy.
}
\arguments{
\item{a}{
output produced by functions \code{hierclust} (or 'hclust' or \code{mclust}).
}}
\value{
matrix of dimensions \code{n} by \code{n-2} giving cluster assignments to the 2, 3, \dots{}
n-1 cluster partitions of the hierarchy.  This corresponds to levels \code{n-2}
to 2.  The observations correspond to the rows, and are in
sequence, \code{1, 2, \dots{}}.  Column \code{j} specifies which of the \code{2, 3, \dots{} j-1}
clusters each observation is associated with at the level \code{j-1} of the 
agglomerative process.  The two clusters which merge in going from level 
\code{j-1} to \code{j-2} are replaced in the \code{j-2}nd column with the lower of the
two cluster sequence numbers.  
}
\section{NOTE}{


The time requirement of 'members' is \code{O(n^3)}.  Sample time for Sun 
SPARCstation 1, with 600 observations: 303 secs.
}
\seealso{
Functions \code{\link{hierclust}}, \code{\link{hclust}}, or \code{\link{mclust}} produce the 
array containing the 
sequence of agglomerations which can be used as input for \code{\link{members}}.
}
\examples{
h <- hierclust(a, method=2)
k <- members(h)
}
\keyword{cluster}
\keyword{multivariate}
% Converted by Sd2Rd version 0.2-a5.
