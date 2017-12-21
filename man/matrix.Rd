\name{simple_triplet_matrix}
\alias{simple_triplet_matrix}
\alias{simple_triplet_zero_matrix}
\alias{simple_triplet_diag_matrix}
\alias{as.simple_triplet_matrix}
\alias{is.simple_triplet_matrix}

\title{Simple Triplet Matrix}
\description{Data structures and operators for sparse matrices based on
  simple triplet representation.}

\usage{
simple_triplet_matrix(i, j, v, nrow = max(i), ncol = max(j),
                      dimnames = NULL)
simple_triplet_zero_matrix(nrow, ncol = nrow, mode = "double")
simple_triplet_diag_matrix(v, nrow = length(v))

as.simple_triplet_matrix(x)
is.simple_triplet_matrix(x)
}
\arguments{
\item{i, j}{Integer vectors of row and column indices, respectively.}
\item{v}{Vector of values.}
\item{nrow, ncol}{Integer values specifying the number of rows and
columns, respectively. Defaults are the maximum row and column indices, respectively.}
\item{dimnames}{A \code{dimnames} attribute for the matrix:
                \code{NULL} or a \code{list} of
                length 2 giving the row and column names respectively.  An
                empty list is treated as \code{NULL},
                and a list of length one as
                row names.  The list can be named, and the list names will be
                used as names for the dimensions.}
\item{mode}{Character string specifying the mode of the values.}
\item{x}{An \R object.}
}
\details{
  \code{simple_triplet_matrix} is a generator for
  a class of \dQuote{lightweight} sparse matrices, \dQuote{simply}
  represented by triplets \code{(i, j, v)} of row indices \code{i}, column
  indices \code{j}, and values \code{v}, respectively.
  \code{simple_triplet_zero_matrix} and
  \code{simple_triplet_diag_matrix} are convenience functions for the
  creation of empty and diagonal matrices.

  Currently implemented operations include the addition, subtraction,
  multiplication and division of compatible simple triplet matrices,
  as well as the multiplication and division of a simple triplet matrix
  and a vector.  Comparisons of the elements of a simple triplet
  matrices with a number are also provided. In addition, 
  methods for indexing, combining by rows
  (\code{rbind}) and columns (\code{cbind}), transposing (\code{t}),
  concatenating (\code{c}), and detecting/extracting duplicated and
  unique rows are implemented.
}
\seealso{
  \code{\link{simple_sparse_array}} for sparse arrays.
}
\examples{
x <- matrix(c(1, 0, 0, 2), nrow = 2)
s <- as.simple_triplet_matrix(x)
identical(x, as.matrix(s))

simple_triplet_matrix(c(1, 4), c(1, 2), c(1, 2))
simple_triplet_zero_matrix(3)
simple_triplet_diag_matrix(1:3)

cbind(rbind(s, t(s)), rbind(s, s))
\dontrun{
## map to default Matrix class
stopifnot(require("Matrix"))
sparseMatrix(i = s$i, j = s$j, x = s$v, dims = dim(s), 
	     dimnames = dimnames(s))
}
}
\keyword{math}
