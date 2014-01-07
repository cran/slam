\name{crossprod}
\alias{tcrossprod_simple_triplet_matrix}
\alias{crossprod_simple_triplet_matrix}
\alias{matprod_simple_triplet_matrix}
\title{
Matrix Crossproduct
}
\description{
Compute the matrix cross-product of a sparse and a dense or sparse matrix.
}
\usage{
tcrossprod_simple_triplet_matrix(x, y = NULL)
##
crossprod_simple_triplet_matrix(x, y = NULL)
matprod_simple_triplet_matrix(x, y)
}
\arguments{
  \item{x, y}{a matrix in \code{simple_triplet_matrix}-form and/or a dense
matrix, where one must be of the form indicated by the suffix.
}
}
\details{
Function \code{tcrossprod_simple_triplet_matrix} implements fast computation 
of \code{x \%*\% t(x)} and \code{x \%*\% t(y)} (\code{tcrossprod}). The
remaining functions are (optimized) wrappers.
}
\value{
A double matrix, with appropriate \code{dimnames} taken from \code{x}
and \code{y}.
}
\author{
Christian Buchta
}
\note{
The computation is delegated to \code{tcrossprod} if \code{y} (or \code{x}
if \code{y == NULL}) contains any of the special values \code{NA},
\code{NaN}, or \code{Inf}.
}
\seealso{
\code{\link[base]{crossprod}} for dense-on-dense computations.
}
\examples{
## 
x <- matrix(c(1, 0, 0, 2, 1, 0), nrow = 3)
x
s <- as.simple_triplet_matrix(x)
tcrossprod_simple_triplet_matrix(s, x)
##
tcrossprod_simple_triplet_matrix(s)
##
tcrossprod_simple_triplet_matrix(s[1L, ], s[2:3, ])
}
\keyword{algebra}
\keyword{array}
