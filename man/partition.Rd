\name{partition}
\alias{partition}
\title{
Partitioning by Iterative Optimization
}
\usage{
partition(a, g, iter.max = 15, option = "mindst", diagnostics = FALSE)
}
\description{
  Returns cluster memberships.  In the case of a data matrix, the user 
specifies either a requested number of clusters; or - implicitly - as many 
clusters as there are rows in an initial guess at the cluster centers.  
Incorporates S function `kmeans'.
}
\arguments{
\item{a}{
matrix of multivariate data.
Each row corresponds to an observation, and each column corresponds to 
a variable.  
}
\item{g}{
Number of groups or clusters.
}
\item{centers}{
matrix of initial guesses for the cluster centers.  Each row represents a
cluster center, and thus \code{centers} must have the same number of columns as
\code{a}.  The number of rows in \code{centers} is the number of clusters that will be
formed.
}
\item{iter.max}{
maximum number of iterations.
}
\item{option}{
Either \code{mindst} (default) or \code{exch}.  Options for the Spaeth (1985) 
algorithm. In the former case, the variances of
the groups are optimized by assigning the objects to groups such that they
are minimally distant from group centers.  In the latter case, the variances
are optimized by exchanging the objects between groups such that they are
minimally distant from group centers. 
}
\item{diagnostics}{
FALSE (default) implies that cluster cardinalities, cluster center coordinates,
and cluster compactness values will not be output to the command window.  
}}
\value{
in the case of the Hartigan algorithm (where initial cluster centers 
are provided) a list with the following components:

\item{cluster}{
vector of integers, ranging from \code{1} to \code{ng} (or number of rows in 
\code{centers}), with length the same as the number of rows of \code{a}.  
The \code{i}th value indicates the group in which the
\code{i}th observation belongs.
}
\item{centers}{
matrix containing the coordinates of the final group centers.
}
\item{withinss}{
vector of length \code{ng} or number of rows in \code{centers}.  The \code{i}th value 
gives the within group sum of squares for the \code{i}th group.
}
\item{size}{
vector of length \code{ng}, or the number of rows in \code{centers}.  
The \code{i}th value gives the number of observations in group \code{i}.


In the case of the Spaeth algorithm (where the number of clusters is given),
the first of these, only, is returned.  
}}
\section{Method}{
Consider first the clustering of coordinate data.
The object is to find a partition of the observations with \code{ng} (or the 
number of rows in \code{centers}) groups that
minimizes the sum of \code{withinss} values.  To actually guarantee the 
minimum would be computationally
infeasible in many settings; this function finds a local minimum, that is,
a solution such that there is no single switch of an observation from one
group to another group that will decrease the objective criterion.  In the
case of the Spaeth (1985) algorithms, the local
minimum arrived at is dependent on the initial arbitrary assignment of 
observations to groups.  This arbitrary assignment is carried out using a
random number generator.  Subsequent executions of the \code{partition} function
will make use of different initial arbitrary assignments.  The function
\code{partition} should therefore be run a number of times, and the best result
kept.


The initial data is not standardized: it may be advisable to divide all column
values by the range of the column's values, or to standardize in some other 
way.


Note that this routine is best used by specifying a small number of groups;
specifying a large number of groups may cause many empty classes.  When this
happens, group centroids and compactness values may be undefined.  This may
also happen when there are very large values in the input data set.  
Cluster memberships, though, can still be used.  


Sample timings for the Spaeth (1985) algorithms (accessed by specifying the 
desired number of groups): 33000 rows, 4 columns, 3 groups took about 40 
secs. on a Sun SPARCstation 1; 50 groups took about 140 secs.  


When deciding on the number of clusters, Hartigan (1975, pp. 90-91) suggests
the following rough rule of thumb.  If \code{k} is the result of \code{partition} with
\code{k} groups, and \code{kplus1} is the result with \code{k+1} groups, then it is 
justifiable to add the extra group when 


(sum(k\$withinss)/sum(kplus1\$withinss)-1)*(nrow(a)-k+1)


is greater than 10.


When a partition is being obtained from a dendrogram, the procedure is
straightforward: a slice is made, orthogonal to the cluster criterion value
axis, and the clusters read off.
}
\references{
Spaeth, H. (1985).  
Cluster Dissection and Analysis: Theory, Fortran Programs, Examples.
Ellis Horwood, Chichester.


Hartigan, J.A. (1975).
Clustering Algorithms.
Wiley, New York.


Hartigan, J.A. and Wong, M.A. (1979).
A k-means clustering algorithm,
Applied Statistics,
vol. 28, pp. 100-108.
}
\seealso{
\code{\link{kmeans}} (incorporated into this function); \code{\link{hierclust}}
(hierarchical clustering routines, incorporating \code{\link{hclust}}).  
}
\examples{
data(iris3)
iris <- rbind(iris3[,,1], iris3[,,2], iris3[,,3])
pp <- partition(iris,3)
# Plot the results in the plane comprising the first two columns of `a'
x <- iris[,1]       # for convenience
y <- iris[,2]
plot(x, y, type="n")     # set up plot scales, etc.
points(x[pp==1], y[pp==1],pch="*")
points(x[pp==2], y[pp==2],pch="+")
points(x[pp==3], y[pp==3],pch="o")


# Secondly, specifying guesses at cluster centers.
irismean <- t(apply(iris3, c(2, 3), 'mean'))
library(mva)
km <- partition(iris, irismean)
wrong <- km$cluster!=rep(1:3, c(50, 50, 50))
plot(iris[,2], iris[,3], type="n")
text(iris[!wrong, 2], iris[!wrong, 3], km$cluster)
# identify cluster membership that is correct
points(iris[wrong, 2], iris[wrong, 3], pch=15, col=2)
# boxes for points in error
title(main="K-Means Clustering of the Iris Data")
}
\keyword{cluster}
\keyword{multivariate}
% Converted by Sd2Rd version 0.2-a5.


