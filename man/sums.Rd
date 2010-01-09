\name{sums}
% NOTE need generics
\alias{row_sums}
\alias{row_sums.simple_triplet_matrix}
\alias{col_sums}
\alias{col_sums.simple_triplet_matrix}
\alias{row_means}
\alias{row_means.simple_triplet_matrix}
\alias{col_means}
\alias{col_means.simple_triplet_matrix}
\title{
Form Row and Column Sums and Means
}
\description{
Form row and column sums and means for sparse arrays (currently
\code{simple_triplet_matrix} only).
}
\usage{
row_sums(x, na.rm = FALSE, dims = 1, \dots)
col_sums(x, na.rm = FALSE, dims = 1, \dots)
row_means(x, na.rm = FALSE, dims = 1, \dots)
col_means(x, na.rm = FALSE, dims = 1, \dots)

\method{row_sums}{simple_triplet_matrix}(x, na.rm = FALSE, dims = 1, \dots)
\method{col_sums}{simple_triplet_matrix}(x, na.rm = FALSE, dims = 1, \dots)
\method{row_means}{simple_triplet_matrix}(x, na.rm = FALSE, dims = 1, \dots)
\method{col_means}{simple_triplet_matrix}(x, na.rm = FALSE, dims = 1, \dots)
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
row_sums(s)
row_sums(s, na.rm = TRUE)
col_sums(s)
col_sums(s, na.rm = TRUE)
}
\keyword{array}
\keyword{algebra}
\keyword{arith}
