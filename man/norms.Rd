\name{norms}
\alias{row_norms}
\alias{col_norms}
\title{Row and Column Norms}
\description{
  Compute row and column \eqn{p}-norms.
}
\usage{
row_norms(x, p = 2)
col_norms(x, p = 2)
}
\arguments{
  \item{x}{a sparse \code{\link{simple_triplet_matrix}}, or a dense
    matrix.}
  \item{p}{a numeric at least one.  Using \code{Inf} gives the maximum
    norm.}
}
\value{
  A vector with the row or column \eqn{p}-norms for the given matrix.
}
\examples{
x <- matrix(1 : 9, 3L)
## Row lengths:
row_norms(x)
## Column maxima:
col_norms(x, Inf)
}


