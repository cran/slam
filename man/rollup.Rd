\name{rollup}
% NOTE need generics
\alias{rollup}
\alias{rollup.array}
\alias{rollup.matrix}
\alias{rollup.simple_triplet_matrix}
\alias{rollup.simple_sparse_array}
\title{
Rollup Sparse Arrays
}
\description{
Rollup (aggregate) sparse arrays along arbitrary dimensions.
}
\usage{
rollup(x, MARGIN, INDEX, FUN, ...)

\method{rollup}{simple_triplet_matrix}(x, MARGIN, INDEX, FUN, \dots)
\method{rollup}{simple_sparse_array}(x, MARGIN, INDEX, FUN, \dots)
\method{rollup}{matrix}(x, MARGIN, INDEX, FUN, \dots)
\method{rollup}{array}(x, MARGIN, INDEX, FUN, \dots)
}
\arguments{
  \item{x}{a sparse array containing numeric, integer, or logical values.
}
  \item{MARGIN}{a vector giving the subscripts of the dimensions to be rolled up.
}
  \item{INDEX}{a corresponding (\code{list} of) \code{factor} (components) in the sense that \code{as.factor(f)} defines the grouping.}
  \item{FUN}{the function to be applied. Note that only functions that return a scalar value are applicable.
}
  \item{\dots}{optional arguments to \code{FUN}.
}
}
\details{
Provides fast summation over the rows or columns of sparse matrices in
\code{simple_triplet}-form.

The default method tries to coerce \code{x} to \code{array}.
}
\value{
An object of the same class as \code{x} where for class
\code{simple_triplet_matrix} the values are always of type \code{double}.
The \code{dimnames} are taken from (the components of) \code{INDEX}.
}
\author{
Christian Buchta
}
\note{
If (a component of) \code{INDEX} contains \code{NA} values the corresponding
positions are ignored.
}
\seealso{
\code{simple_triplet_matrix} and \code{simple_sparse_array} for sparse
arrays.
}
\examples{
##
x <- matrix(c(1, 0, 0, 2, 1, NA), nrow = 2, 
	    dimnames = list(A = 1:2, B = 1:3))
x
zz <- rollup(x, 2L, c(1,2,1), sum, na.rm = TRUE)
zz
##
s <- as.simple_triplet_matrix(x)
z <- rollup(s, 2L, c(1,2,1), sum, na.rm = TRUE)
identical(zz, as.matrix(z))
##
a <- as.simple_sparse_array(s)
z <- rollup(a, 2L, c(1,2,1), sum, na.rm = TRUE)
identical(as.array(zz), as.array(z))
##
zz <- rollup(x, 2L, c(1,NA,1), sum, na.rm = TRUE)
zz
##
z  <- rollup(s, 2L, c(1,NA,1), sum, na.rm = TRUE)
identical(zz, as.matrix(z))
##
z  <- rollup(a, 2L, c(1,NA,1), sum, na.rm = TRUE)
identical(as.array(zz), as.array(z))
}
\keyword{array}
\keyword{algebra}
\keyword{arith}
