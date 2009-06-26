\name{simple_sparse_array}
\alias{simple_sparse_array}
\alias{as.simple_sparse_array}
\alias{is.simple_sparse_array}

\title{Simple Sparse Arrays}
\description{Data structures and operators for sparse arrays based on
  a representation by index matrix and value vector.}

\usage{
simple_sparse_array(i, v, dim = NULL, dimnames = NULL)

as.simple_sparse_array(x)
is.simple_sparse_array(x)
}
\arguments{
\item{i}{Integer matrix of array indices.}
\item{v}{Vector of values.}
\item{dim}{Integer vector specifying the size of the dimensions.}
\item{dimnames}{either \code{NULL} or the names for the dimensions.
          This is a list
          with one component for each dimension, either \code{NULL} or a
          character vector of the length given by \code{dim} for that
          dimension.  The list can be named, and the list names will be
          used as names for the dimensions.  If the list is shorter
          than the number of dimensions, it is extended by \code{NULL}'s to
          the length required.}
\item{x}{An \R object.}
}
\details{
  \code{simple_sparse_array} is a generator for
  a class of \dQuote{lightweight} sparse arrays, 
  represented by index matrices and value vectors.
  Currently, only methods for indexing and coercion are implemented.
}
\seealso{
  \code{\link{simple_sparse_array}} for sparse arrays.
}
\examples{
x <- array(c(1, 0, 0, 2, 0, 0, 0, 3), dim = c(2, 2, 2))
s <- as.simple_sparse_array(x)
identical(x, as.array(s))

simple_sparse_array(matrix(c(1, 1, 1, 3, 3, 3), nrow = 2), c(1, 2))
}
\keyword{math}
