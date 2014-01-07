\name{apply}
\alias{rowapply_simple_triplet_matrix}
\alias{colapply_simple_triplet_matrix}
\alias{crossapply_simple_triplet_matrix}
\alias{tcrossapply_simple_triplet_matrix}
\title{
Apply Functions Over Sparse Matrix Margins
}
\description{
Apply functions to (the cross-pairs of) the rows or columns of a sparse
matrix.
}
\usage{
rowapply_simple_triplet_matrix(x, FUN, ...)
colapply_simple_triplet_matrix(x, FUN, ...)

crossapply_simple_triplet_matrix(x, y = NULL, FUN, ...) 
tcrossapply_simple_triplet_matrix(x, y = NULL, FUN, ...)
}
\arguments{
\item{x, y}{a matrix in \code{simple_triplet_matrix}-form or, one of
\code{x} and \code{y}, of class \code{matrix}.
}
  \item{FUN}{the name of the function to be applied.
}
  \item{\dots}{optional arguments to \code{FUN}.
}
}
\details{
\code{colapply_simple_triplet_matrix} temporarily expands each column of
\code{x} to dense \code{vector} representation and applies the function
specified in \code{FUN}.

\code{crossapply_simple_triplet_matrix} temporarily expands each cross-pair
of columns of \code{x} (and \code{y}) to dense \code{vector} representation
and applies the function specified in \code{FUN}.

Note that if \code{y = NULL} then only the entries in the lower triangle
and the diagonal are computed, assuming that \code{FUN} is symmetric.
}
\value{
A \code{vector} (\code{matrix}) of length (dimensionality) of the margin(s)
used. The type depends on the result of \code{FUN}.

Note that the result of \code{colapply_simple_triplet_matrix} is never
simplified to \code{matrix}.
}
\author{
Christian Buchta
}
%\note{
%}
\seealso{
\code{\link[base]{apply}} for dense-on-dense computations.
}
\examples{
## 
x <- matrix(c(1, 0, 0, 2, 1, 0), nrow = 3, 
    dimnames = list(1:3, LETTERS[1:2]))
x
s <- as.simple_triplet_matrix(x)
colapply_simple_triplet_matrix(s, FUN = var)
##
simplify2array(colapply_simple_triplet_matrix(s, identity))
##
crossapply_simple_triplet_matrix(s, FUN = var)
}
\keyword{algebra}
\keyword{array}
