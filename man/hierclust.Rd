\name{hierclust}
\alias{hierclust}
\alias{hc}
\title{
Hierarchical Clustering
}
\usage{
hierclust(a, method = 1, bign = FALSE, diagnostics = FALSE, sim = " ",
             movie = FALSE, option = "thru", alpha = 1.0, repres = "chull", 
             show = "change")

hierclust(dst, method="compact", sim=)


}
\description{
  Agglomerates the closest pair of observations, and replaces them with
a cluster.  Agglomerates the next closest pair of observations or clusters.
Continues through `n-1' agglomerations, assuming there are `n' observations.
}
\arguments{
\item{a}{
data matrix to be analyzed, the rows representing observations and the
columns variables; or, alternatively, lower-diagonal dissimilarity matrix in
vector or list format.  In latter case, the number of values is \code{n(n-1)/2} 
where \code{n} is the number of observations.


(See below for second variant of usage of this function.)


First variant:


}
\item{method}{
First possibility: integer taking the value between 1 and 7.  These methods 
correspond to:
Ward's minimum variance or error sum of squares method (1); single linkage
or nearest neighbor method (2); complete linkage or diameter (3); average
linkage, group average, or UPGMA method (4); McQuitty's or WPGMA method
(5); median, Gower's or WPGMC method (6); and centroid or UPGMC method (7).
If the \code{movie} option is used, only \code{method}s 1 to 4 are supported.
Second possibility: string \code{compact} (default; complete linkage), 
\code{average} (average linkage), or \code{connected} (single linkage).
}
\item{sim}{
structure giving similarities rather than dissimilarities.  This can either
be a symmetric matrix or a vector with a \code{Size} attribute.  Only one of 
dissimilarities or similarities may be specified (or, of course, as an
alternative, a data matrix.
}
\item{bign}{
Is \code{n} big?  If storage is problemsome, a different implementation of the 
Ward criterion may be tried.  This determines dissimilarities on the fly, and
hence requires \code{O(n)} storage.
}
\item{diagnostics}{
show, or suppress, some information directed to the command window. 
Suppression is useful if functional programming capabilities are availed of,
leading to lazy evaluation; this often results in multiple evaluations, 
with consequent repetition of diagnostic information.
}
\item{movie}{
whether interactive display of construction of hierarchy is wanted.  This
necessitates an open plot window.  It does not produce an output object.  
All remaining arguments apply only in the case of the \code{movie=TRUE}.
}
\item{option}{
(Only for \code{movie}.) "thru" (default) or "prompt": Go right through all n-1 
agglomerations, or ask
the user whether to continue at each agglomeration?  Here is what we
recommend: On the first occasion, go right through.  Check out the cluster
criterion values reported on in the command window.  Then using the 'prompt'
option, go as far as the desired partition, and print it out.
}
\item{alpha}{
(Only for \code{movie}.) 1.0 (default) or 0.0 < alpha <= 1.0.
Only used by the minimum variance agglomerative criterion. Alpha downweights
the variance component on the principal axis of the cluster, and thereby
redefines the cluster criterion such that linear clusters are preferred.  A
value of alpha = 1.0 gives exactly the traditional Ward's minimum variance
method.  See Murtagh and Raftery (1984) for more on this altered criterion.
}
\item{repres}{
(Only for \code{movie}.) repres = "chull" (default), "lines", "gauss":
Representation of clusters: For input data which has dimensionality greater
than 2, only "chull" is allowed (since the representations were found to be
too problematic, given that the clustering takes place in m-dimensions, and
the the representation is necessarily a 2-dimensional projection of this).
Representation "chull" = convex hull, or straight line in case
of collinear points.  (We test for collinearity using principal
components analysis.)
Representation "lines" = lines of best fit (least squares fit is used,
and the extremities are defined by the intersection of the min. and max. x
values on this LS fit line). "gauss" = 1-sigma circles.  Side-effect of
latter: warning messages indicate whenever circles overlap plot boundaries.
}
\item{show}{
(Only for \code{movie}.) "change" or "all":
In this method's present implementation, the clusters at each stage are
\code{highlighted} by means of a different symbol.  Option show = "change" has
this highlighting changed at each agglomeration, so that only the most recent
agglomerands are shown.  Option show = "all" does not turn off the
highlighting in this way.


Second variant (which re-wraps S function \code{hclust}: cf. help file on \code{hclust}):




Exactly one of dst or sim must be specified.


}
\item{dst}{
a distance structure or distance matrix.  Normally this
will be the result of the function \code{dist}, but it can be any
data of the form returned by \code{dist}, or  a  full,  symmetric
matrix.  Missing values are not allowed.
}
\item{method}{
character string giving  the  clustering  method.   The
three  methods  currently implemented are "average", "con-
nected" (single linkage) and "compact" (complete linkage).
(The first three characters of the method are sufficient.)
}
\item{sim=}{
structure giving similarities  rather  than  distances.
This  can  either be a symmetric matrix or a vector with a
"Size" attribute.  Missing values are not allowed.
}}
\value{
When \code{movie} is not specified in the first variant, this is a 
list describing the hierarchical clustering of the observations.
The \code{movie} option is considerably slower than alternative options, due
to the processing involved.  Hence its use is recommended only on small
data sets (less than 100 observations).

\item{merge}{
an \code{n-1} by 2 matrix. Row \code{i} of \code{merge} describes the merging of clusters
at step \code{i} of the clustering.  If an element \code{j} in the row is negative,
then observation \code{-j} was merged at this stage.  If \code{j} is positive 
then the merge
was with the cluster formed at the (earlier) stage \code{j} of the algorithm. Thus
negative entries in \code{merge} indicate agglomerations of singletons, and 
positive entries indicate agglomerations of non-singletons.
}
\item{height}{
a set of \code{n-1} non-decreasing real values.  The
clustering "height": that is, the value of the criterion associated with
the clustering \code{method} for the particular agglomeration.
}
\item{order}{
a vector giving the permutation of the original observations suitable for 
plotting,
in the sense that a cluster plot using this ordering and matrix \code{merge} will
not have crossings of the branches.


In hierarchical cluster displays, a decision is needed at each merge to
specify which subtree should go on the left and which on the right.  Since, 
for 
\code{n} observations there are \code{n-1} merges, there are \code{2^(n-1)} possible
orderings for the leaves in a cluster tree, or dendrogram.  The default
algorithm in \code{hc} is to order the subtree so that the tighter cluster is on
the left (the last, i.e. most recent, merge of the left subtree is at a lower 
value than the last
merge of the right subtree).  Observations are the tightest clusters possible,
and merges involving two observations place them in order by their 
observation sequence number.
}}
\section{Note}{
The squared Euclidean distance is used in the \code{movie} case, and in the case
when a data matrix is input. Note that the \code{movie} option returns 
invisibly.   
}
\section{Method}{
In the case of all methods when \code{bign = FALSE} (default), 
the Lance-Williams dissimilarity update formula is used to allow the 7 methods
to be implemented.   


In the direct data matrix case, 
a list of nearest neighbors is maintained, which 
allows the closest pair of clusters or observations to be determined by
scanning this list to find the pair of agglomerands on each iteration.
Hence \code{O(n)} operations are required for each of \code{n-1} iterations, i.e.
\code{O(n^2)} computation time in total.  The 
initial setting up of the \code{n} nearest neighbors also requires \code{O(n^2)} time.
The processing to be carried out following each agglomeration is \code{O(n)}.
In total, all methods implemented in \code{hc} are of time complexity \code{O(n^2)}.
Unless \code{bign = TRUE} (see below), storage 
complexity is \code{O(n^2)}.  Typical times on a Sun SPARCstation
1 for 400x4 and 800x4 arrays: 12 and 48 secs.


Ward's minimum variance method aims at finding compact, spherical clusters.
The complete link method finds similar clusters. The single link method,
which is closely related to the minimal spanning tree, adopts a 
\code{friends of friends} clustering strategy.  
The other methods can be regarded as aiming
for clusters with characteristics somewhere between the single and complete
link methods. 


In the case of the Ward method, and \code{bign = TRUE}, the nearest-neighbor 
chain algorithm is used, which gives an \code{O(n^2)}
computation time algorithm. The input matrix is stored, but not the
dissimilarities produced.  Computation time is marginally inferior to the
implementation of Ward's method when \code{bign = FALSE} (default).  Sample 
times for a 400x4 and an 800x4
array on a Sun SPARCstation 1: 12.5 and 50 secs.
}
\references{
A small sample: 


B. Everitt, 
Cluster Analysis
Heinemann Educ. Books, London 1974; 


P.H.A. Sneath and R.R. Sokal, 
Numerical Taxonomy
Freeman, San Francisco, 1973; 


M.R. Anderberg, 
Cluster Analysis for Applications
Academic Press, New York, 1973;


A.D. Gordon, 
Classification
Chapman and Hall, London, 1981; and


F. Murtagh, 
Multidimensional Clustering Algorithms
COMPSTAT Lectures 4,
Physica-Verlag, Wuerzburg, 1985 (for algorithmic details of algorithms used).
}
\seealso{
For finding assignments of objects to clusters, use \code{\link{members}}.
'Dist' is used for creating (dis)similarities.  Functions
\code{\link{cutree}}, and \code{\link{dist}} are already available in R.  Function \code{\link{members}}
is available in this function suite.
}
\examples{
data(iris)
iris <- as.matrix(iris[,1:4])
# Default Ward criterion
hierclust(iris, movie = TRUE)
h <- hierclust(iris)
k <- members(h)
}
\keyword{cluster}
\keyword{multivariate}
% Converted by Sd2Rd version 0.2-a5.
