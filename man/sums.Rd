\name{sums}
% NOTE need generics
\alias{rowSums}
\alias{rowSums.simple_triplet_matrix}
\alias{colSums}
\alias{colSums.simple_triplet_matrix}
\alias{rowMeans}
\alias{rowMeans.simple_triplet_matrix}
\alias{colMeans}
\alias{colMeans.simple_triplet_matrix}
\title{
Form Row and Column Sums and Means
}
\description{
Form row and column sums and means for sparse arrays (currently
\code{simple_triplet_matrix} only).
}
\usage{
rowSums(x, na.rm = FALSE, dims = 1, \dots)
colSums(x, na.rm = FALSE, dims = 1, \dots)
rowMeans(x, na.rm = FALSE, dims = 1, \dots)
colMeans(x, na.rm = FALSE, dims = 1, \dots)

\method{rowSums}{simple_triplet_matrix}(x, na.rm = FALSE, dims = 1, \dots)
\method{colSums}{simple_triplet_matrix}(x, na.rm = FALSE, dims = 1, \dots)
\method{rowMeans}{simple_triplet_matrix}(x, na.rm = FALSE, dims = 1, \dots)
\method{colMeans}{simple_triplet_matrix}(x, na.rm = FALSE, dims = 1, \dots)
}
\arguments{
  \item{x}{a sparse array containing numeric, integer, or logical values.
}
  \item{na.rm}{logical. Should missing values (including \code{NaN})
    be omitted from the calculations?
}
  \item{dims}{currently not used for sparse arrays.
}
  \item{\dots}{currently not used for sparse arrays.
}
}
\details{
Provides fast summation over the rows or columns of sparse matrices in
\code{simple_triplet}-form.
}
\value{
A numeric (double) array of suitable size, or a vector if the result
is one-dimensional. The \code{dimnames} (or \code{names} for a vector
result) are taken from the original array.
}
\author{
Christian Buchta
}
\note{
Results are always of storage type \code{double} to avoid (integer)
overflows.
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
rowSums(s)
rowSums(s, na.rm = TRUE)
colSums(s)
colSums(s, na.rm = TRUE)
}
\keyword{array}
\keyword{algebra}
\keyword{arith}
