\name{abind}
\alias{abind_simple_sparse_array}
\alias{extend_simple_sparse_array}

\title{Combine Sparse Arrays}
\description{Combine a sequence of (sparse) arrays, matrices, or vectors
  into a single sparse array of the same or higher dimension.}

\usage{

abind_simple_sparse_array(..., MARGIN = 1L)
extend_simple_sparse_array(x, MARGIN = 0L)
}
\arguments{
\item{\dots}{\R objects of (or coercible to) class \code{simple_sparse_array}.}
\item{MARGIN}{The dimension along which to bind the arrays.}
\item{x}{An object of class \code{simple_sparse_array}.}
}
\details{
  \code{abind_simple_sparse_array} automatically extends the dimensions
  of the elements of \sQuote{\dots} before it combines them along the
  dimension specified in \code{MARGIN}. If a negative value is specified
  first all elements are extended left of the target dimension.
  
  \code{extend_simple_sparse_array} inserts one (or more) one-level
  dimension(s) into \code{x} to the right of the position(s) specified in
  \code{MARGIN}, or to the left if specified in negative terms. Note that
  the target positions must all be in the range of the dimensions of \code{x}
  (see Examples).
}
\value{
  An object of class \code{simple_sparse_array} where the \code{dimnames}
  are taken from the elements of \sQuote{\dots}.
}
\author{
Christian Buchta
}
\seealso{
  \code{\link{simple_sparse_array}} for sparse arrays.
}
\examples{
## automatic
abind_simple_sparse_array(1:3, array(4:6, c(1,3)))
abind_simple_sparse_array(1:3, array(4:6, c(3,1)), MARGIN = 2L)

## manual
abind_simple_sparse_array(1:3, 4:6)
abind_simple_sparse_array(1:3, 4:6, MARGIN = -2L)   ## by columns
abind_simple_sparse_array(1:3, 4:6, MARGIN = -1L)   ## by rows

##
a <- as.simple_sparse_array(1:3)
a
extend_simple_sparse_array(a, c( 0L, 1L))
extend_simple_sparse_array(a, c(-1L,-2L))   ## the same
extend_simple_sparse_array(a, c( 1L, 1L))
}
\keyword{array}
