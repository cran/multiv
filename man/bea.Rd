\name{bea}
\alias{bea}
\title{
Bond Energy Algorithm 
}
\usage{
bea(a, istart = 0, jstart = 0)
}
\description{
  Permutes  rows  and columns of an array, in order to maximize
  proximity of large-valued array elements. 
}
\arguments{
\item{a}{
data matrix to be analyzed.  The rows could represent observations and the
columns variables (a two-way, two-mode array); or the rows and columns could
both represent observations (a two-way, one-mode array).  Missing values are
not supported.
}
\item{istart}{
sequence number of first row to be placed; if not specified, then a row is 
arbitrarily chosen.
}
\item{jstart}{
sequence number of first column to be placed; if not specified, then a row is
arbitrarily chosen.
}}
\value{
list describing the reordering carried out:

\item{b}{
row and column permuted array.
}
\item{ib}{
permutation of rows.
}
\item{jb}{
permutation of columns.
}
\item{e}{
\code{bond energy} of the permuted matrix, \code{b}.
}}
\section{Note}{
This is a \code{non-destructive} approach to analyzing data, insofar as the 
original data is not altered; it is only rearranged, in order to highlight
potentially interesting aspects of the data.  Subsequent runs of this 
routine may improve the bond energy (see example below).  
}
\section{Method}{
A row is arbitrarily placed; then rows are positioned one by one.  When this
is completed, the columns are treated similarly.  The overall procedure 
amounts to two approximate traveling salesman problems, - one on the rows and 
one on the columns.  The so-called \code{best insertion} strategy is used: rows 
(or columns) are inserted into the current permuted list of rows (or columns).


To have repeatable results, control the random number generator used for the
initial selection of rows and columns with the \code{set.seed} command,
before calling the \code{bea} function.  


No methods are currently defined for the objects produced by the \code{bea}
function.
}
\section{Background}{
This simple method has been used in operations research, production 
engineering, marketing, and various other fields.  Arabie and Hubert (1990)
recommend that it be used with ratio scale data; and they question its use 
with non-binary data if the objective is to find a seriation or 
one-dimensional ordering of rows and columns.
}
\references{


W.T. McCormick, P.J. Schweitzer and T.W. White, 
\code{Problem decomposition and data reorganization by a clustering technique}, 
Operations Research, 
vol. 20, pp. 993-1009, Sept./Oct. 1972.


P. Arabie and L.J. Hubert, 
\code{The bond energy algorithm revisited}, 
IEEE Transactions on Systems, Man, and Cybernetics, 
vol. 20, pp. 268-274, 1990.


P. Arabie, S. Schleutermann, J. Daws and L. Hubert, 
`Marketing applications of sequencing and partitioning of nonsymmetric 
and/or two-mode matrices', 
in W. Gaul and M. Schader, Eds., 
Data Analysis, Decision Support, \& Expert Knowledge Representation in 
Marketing,
Springer Verlag, 1988, pp. 215-224.
}
\examples{
data(USArrests)
USArrests <- as.matrix(USArrests)
run1 <- bea(USArrests)
# look at energy:
run1$e
# Redo, on basis of first run, to see if energy increases:
run2 <- bea(run1$b)
run2$e
# Remark: to have repeatable results, issue the `set.seed' before the `bea'
# command.  Now, reorder the rows and the columns once more:
run3 <- bea(run2$b)
run3$e
# Of course, sequencing of rows and columns of each run is relative to 
# the array which was input.  Get `net' ordering of rows after third run:
run1$ib[run2$ib[run3$ib]]
# Plot the output as images:
image(1:nrow(USArrests), 1:ncol(USArrests), USArrests)
image(1:nrow(USArrests), 1:ncol(USArrests), bea(USArrests)$b)
# Remark: one could split the screen, and use color in such displays of
# the data.
}
\keyword{multivariate}
\keyword{cluster}
% Converted by Sd2Rd version 0.2-a5.
