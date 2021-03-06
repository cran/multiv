\name{logique}
\alias{logique}
\title{
Logical Coding
}
\usage{
logique(x)
}
\description{
  Simple logical coding of a vector: each value in the 
vector is replaced by a 1 (if it is above or equal to the median),
by a 0 (if it is below the median).
}
\arguments{
\item{x}{
real-valued vector, with no missing values.
}}
\value{
matrix of \code{length(a)} rows, and two columns.  The first column contains the
logically coded values of \code{a}, and the second column contains their 
complements.  Hence each row of this returned matrix necessarily sums to 1.
}
\section{Background}{
This form of coding is suitable for a subsequent correspondence analysis.
When all variable have been logically (or fuzzily) coded, the row masses 
(proportional to the row sums) are identical.  Logical coding results in
the input being in complete disjunctive form.
}
\references{
J.-P. Benzecri
Correspondence Analysis Handbook
Marcel Dekker, Basel, 1992.
}
\seealso{
\code{\link{flou}}, \code{\link{ca}}, \code{\link{supplr}}, \code{\link{supplc}}.  
}
\examples{
data(USArrests)
USArrests <- as.matrix(USArrests)
# Fuzzy coding of input variables
a.log <- logique(USArrests[,1])
b.log <- logique(USArrests[,2])
c.log <- logique(USArrests[,3])
d.log <- logique(USArrests[,4])
newdata <- cbind(a.log, b.log, c.log, d.log)
ca.newdata <- ca(newdata)
}
\keyword{multivariate}
\keyword{algebra}
% Converted by Sd2Rd version 0.2-a5.


