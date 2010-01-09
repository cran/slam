\name{crossprod}
\alias{tcrossprod_simple_triplet_matrix}
\title{
Matrix Crossproduct
}
\description{
Compute the matrix cross-product of a sparse and dense matrix.
}
\usage{
tcrossprod_simple_triplet_matrix(x, y = NULL)
}
\arguments{
  \item{x}{a matrix in \code{simple_triplet_matrix}-form.
}
  \item{y}{a numeric matrix.
}
}
\details{
Provides fast computation of \code{x \%*\% t(y)} (\code{tcrossprod}).
}
\value{
A double matrix, with appropriate \code{dimnames} taken from \code{x}
and \code{y}.
}
\author{
Christian Buchta
}
\note{
If \code{y = NULL} or contains any of the special values \code{NA},
\code{NaN}, or \code{Inf} then \code{x} is coerced to \code{matrix} and
the computation is delegated to \code{tcrossprod}.
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
## see note
tcrossprod_simple_triplet_matrix(s)
}
\keyword{algebra}
\keyword{array}
