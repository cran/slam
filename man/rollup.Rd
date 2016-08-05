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

\method{rollup}{simple_triplet_matrix}(x, MARGIN, INDEX = NULL, FUN = sum, \dots,
       REDUCE = FALSE)
\method{rollup}{simple_sparse_array}(x, MARGIN, INDEX = NULL, FUN = sum, \dots, 
       DROP = FALSE, EXPAND = c("none", "sparse", "dense", "all"), 
       MODE = "double")
\method{rollup}{matrix}(x, MARGIN, INDEX = NULL, FUN = sum, \dots,
       DROP = FALSE, MODE = "double")
\method{rollup}{array}(x, MARGIN, INDEX = NULL, FUN = sum, \dots, 
       DROP = FALSE, MODE = "double")
}
\arguments{
  \item{x}{a sparse (or dense) array, typically of numeric or logical values.
}
  \item{MARGIN}{a vector giving the subscripts (names) of the dimensions to be rolled up.
}
  \item{INDEX}{a corresponding (\code{list} of) \code{factor} (components) in the sense that \code{as.factor} defines the grouping(s). If \code{NULL} collapses the corresponding dimension(s) (default).}
  \item{FUN}{the name of the function to be applied.
}
  \item{REDUCE}{option to remove zeros from the result.}
  \item{DROP}{option to delete the dimensions of the result which have only one level.}
  \item{EXPAND}{the cell expansion method to use (see Details).
}
  \item{MODE}{the type to use for the values if the result is empty.
}
  \item{\dots}{optional arguments to \code{FUN}.
}
}
\details{
Provides aggregation of sparse and dense arrays, in particular fast summation over the rows or columns of sparse matrices in
\code{simple_triplet}-form.

If (a component of) \code{INDEX} contains \code{NA} values the
corresponding parts of \code{x} are omitted.

For \code{simple_sparse_array} the following cell expansion methods are
provided:
\describe{
    \item{\code{none}:}{The \emph{non-zero} entries of a cell, if any, are
	supplied to \code{FUN} as a \code{vector}.}
    \item{\code{sparse}:}{The number of \emph{zero} entries of a cell is 
	supplied in addition to above, as a second argument.}
    \item{\code{dense}:}{Cells with \emph{non-zero} entries are expanded to
	a dense \code{array} and supplied to \code{FUN}.}
    \item{\code{all}:}{All cells are expanded to a dense \code{array}
	and supplied to \code{FUN}.}
}
Note that the memory and time consumption increases with the level
of expansion.

Note that the default method tries to coerce \code{x} to \code{array}.
}
\value{
An object of the same class as \code{x} where for class
\code{simple_triplet_matrix} the values are always of type \code{double}
if \code{FUN = sum} (default).

The \code{dimnames} corresponding to \code{MARGIN} are based on (the 
components of) \code{INDEX}.
}
\author{
Christian Buchta
}
\note{
Currently most of the code is written in R and, therefore, the memory
and time it consumes is not optimal.
}
\seealso{
\code{\link{simple_triplet_matrix}} and \code{\link{simple_sparse_array}}
for sparse arrays.

\code{\link{apply}} for dense arrays.
}
\examples{
##
x <- matrix(c(1, 0, 0, 2, 1, NA), nrow = 2, 
	    dimnames = list(A = 1:2, B = 1:3))
x
apply(x, 1L, sum, na.rm = TRUE)
##
rollup(x, 2L, na.rm = TRUE)
##
rollup(x, 2L, c(1,2,1), na.rm = TRUE)
## omit
rollup(x, 2L, c(1,NA,1), na.rm = TRUE)
## expand
a <- as.simple_sparse_array(x)
a
r <- rollup(a, 1L, FUN = mean, na.rm = TRUE, EXPAND = "dense")
as.array(r)
## 
r <- rollup(a, 1L, FUN = function(x, nz) 
	length(x) / (length(x) + nz), 
    EXPAND = "sparse"
)
as.array(r)
}
\keyword{array}
\keyword{algebra}
\keyword{arith}
