\name{tsums}
% NOTE need generics
\alias{row_tsums}
\alias{row_tsums.matrix}
\alias{row_tsums.simple_triplet_matrix}
\alias{col_tsums}
\alias{col_tsums.matrix}
\alias{col_tsums.simple_triplet_matrix}
\title{
Form Ragged Row and Column Sums
}
\description{
Form ragged row and column sums for sparse arrays (currently
\code{simple_triplet_matrix} and \code{matrix} only).
}
\usage{
row_tsums(x, INDEX, na.rm = FALSE, dims = 1, \dots)
col_tsums(x, INDEX, na.rm = FALSE, dims = 1, \dots)

\method{row_tsums}{matrix}(x, INDEX, na.rm = FALSE, dims = 1, \dots)
\method{col_tsums}{matrix}(x, INDEX, na.rm = FALSE, dims = 1, \dots)
\method{row_tsums}{simple_triplet_matrix}(x, INDEX, na.rm = FALSE, dims = 1, \dots)
\method{col_tsums}{simple_triplet_matrix}(x, INDEX, na.rm = FALSE, dims = 1, \dots)
}
\arguments{
  \item{x}{a sparse array containing numeric, integer, or logical values.
}
  \item{INDEX}{a \code{factor} in the sense that \code{as.factor(f)} defines
    the grouping.}
  \item{na.rm}{logical. Should missing values (including \code{NaN})
    be omitted from the calculations?
}
  \item{dims}{currently not used for sparse arrays.
}
  \item{\dots}{currently not used for sparse arrays.
}
}
\details{
Provides fast ragged summation over the rows or columns of sparse matrices in
\code{simple_triplet}-form.
}
\value{
An object of the same class as \code{x} where for class
\code{simple_triplet_matrix} the values are always of type \code{double}.
The \code{rownames} (or \code{colnames}) are taken from \code{INDEX}.
}
\author{
Christian Buchta
}
\note{
If \code{INDEX} contains \code{NA} values the corresponding columns are
ignored.
}
\seealso{
\code{simple_triplet_matrix}, \code{\link[base]{colSums}} for dense
numeric arrays.
}
\examples{
##
x <- matrix(c(1, 0, 0, 2, 1, NA), nrow = 3)
x
s <- as.simple_triplet_matrix(x)
col_tsums(s, c(1,2,1))
z <- col_tsums(s, c(1,2,1), na.rm = TRUE)
z
## check
zz <- col_tsums(as.matrix(s), c(1,2,1), na.rm = TRUE)
identical(zz, as.matrix(z))
}
\keyword{array}
\keyword{algebra}
\keyword{arith}
