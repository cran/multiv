\name{flou}
\alias{flou}
\title{
Fuzzy Coding (3-Way)
}
\description{
  Simple fuzzy, or piecewise linear, coding of a vector: each value in the 
  vector is replaced by a 1 (if it is above or equal to the 67th quantile),
  by a 0 (if it is below or equal to the 33rd quantile), and by a linearly
  interpolated value between 0 and 1 (if it lies between the 33rd and 67th
  quantiles).  
}
\usage{
flou(a)
}
\arguments{
\item{a}{
real-valued vector, with no missing values.
}}
\value{
matrix of \code{length(a)} rows, and two columns.  The first column contains the
fuzzily coded values of \code{a}, and the second column contains their 
complements.  Hence each row of this returned matrix necessarily sums to 1.
}
\section{Background}{
This form of coding is suitable for a subsequent correspondence analysis.
When all variable have been fuzzily (or logically) coded, the row masses 
(proportional to the row sums) are identical.  
}
\references{
J.-P. Benzecri
Correspondence Analysis Handbook
Marcel Dekker, Basel, 1992.


F.J. Gallego,
Codage flou en analyse des correspondances,
Les Cahiers de l'Analyse des Donnees
vol. VII, 413-430, 1982
}
\seealso{
\code{\link{logique}}, \code{\link{ca}}, \code{\link{supplr}}, \code{\link{supplc}}.  
}
\examples{
data(USArrests)
USArrests <- as.matrix(USArrests)
# Fuzzy coding of input variables
a.fuzz <- flou(USArrests[,1])
b.fuzz <- flou(USArrests[,2])
c.fuzz <- flou(USArrests[,3])
d.fuzz <- flou(USArrests[,4])
newdata <- cbind(a.fuzz, b.fuzz, c.fuzz, d.fuzz)
ca.newdata <- ca(newdata)
}
\keyword{multivariate}
\keyword{algebra}
% Converted by Sd2Rd version 0.2-a5.


