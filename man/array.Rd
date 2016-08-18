\name{simple_sparse_array}
\alias{simple_sparse_array}
\alias{simple_sparse_zero_array}
\alias{as.simple_sparse_array}
\alias{is.simple_sparse_array}

\alias{simplify_simple_sparse_array}
\alias{reduce_simple_sparse_array}
\alias{drop_simple_sparse_array}

\title{Simple Sparse Arrays}
\description{Data structures and operators for sparse arrays based on
  a representation by index matrix and value vector.}

\usage{
simple_sparse_array(i, v, dim = NULL, dimnames = NULL)

as.simple_sparse_array(x)
is.simple_sparse_array(x)

simplify_simple_sparse_array(x, higher = TRUE)
reduce_simple_sparse_array(x, strict = FALSE, order = FALSE)
drop_simple_sparse_array(x)
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
\item{x}{An \R object; an object of class \code{simple_sparse_array}
	  (see Note).}
\item{higher}{Option to use the dimensions of the values (see Note).}
\item{strict}{Option to treat violations of sparse representation as error
	  (see Note).}
\item{order}{Option to reorder elements (see Note).}
}
\details{
  \code{simple_sparse_array} is a generator for
  a class of \dQuote{lightweight} sparse arrays, 
  represented by index matrices and value vectors.
  Currently, only methods for indexing and coercion are implemented.
}
\note{
  The \emph{zero} element is defined as \code{vector(typeof(v), 1L)},
  for example, \code{FALSE} for \code{logical} values (see
  \code{\link{vector}}). Clearly, sparse arrays should not contain
  \emph{zero} elements, however, for performance reasons the class
  generator does not remove them.

  If \code{strict = FALSE} (default) \code{reduce_simple_sparse_array} tries
  to repair violations of sparse representation (\emph{zero, multiple}
  elements), otherwise it stops. If \code{order = TRUE} the elements are
  further reordered (see \code{\link{array}}).

  \code{simplify_simple_sparse_array} tries to reduce \code{v}. If
  \code{higher = TRUE} (default) augments \code{x} by the common dimensions
  of \code{v} (from the left), or the common length. Note that \emph{scalar}
  elements are never extended and unused dimensions never dropped.

  \code{drop_simple_sparse_array} drops unused dimensions.

  If \code{prod(dim(x)) > slam_options("max_dense")} empty and negative
  indexing are disabled
  for \code{[} and \code{[<-}. Further, non-negative single (vector)
  indexing is limited to 52 bits of representation.
}
\seealso{
  \code{\link{simple_triplet_matrix}} for sparse matrices.

  \code{\link{slam_options}} for options.
}
\examples{
x <- array(c(1, 0, 0, 2, 0, 0, 0, 3), dim = c(2, 2, 2))
s <- as.simple_sparse_array(x)
identical(x, as.array(s))

simple_sparse_array(matrix(c(1, 3, 1, 3, 1, 3), nrow = 2), c(1, 2))
}
\keyword{math}
