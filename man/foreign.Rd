\name{foreign}
\alias{read_stm_CLUTO}
\alias{write_stm_CLUTO}
\title{Read and Write Sparse Matrix Format Files}
\description{
  Read and write CLUTO sparse matrix format files.
}
\usage{
read_stm_CLUTO(file)
write_stm_CLUTO(x, file)
}
\arguments{
  \item{file}{a character string with the name of the file to read or
    write.}
  \item{x}{a matrix object.}
}
\details{
  Documentation for CLUTO including its sparse matrix format is
  available from \url{http://www-users.cs.umn.edu/~karypis/cluto/}.

  \code{read_stm_CLUTO} reads CLUTO sparse matrices, returning a
  \link[slam:matrix]{simple triplet matrix}.

  \code{write_stm_CLUTO} writes CLUTO sparse matrices.
  Argument \code{x} must be coercible to a simple triplet matrix via
  \code{\link[slam:matrix]{as.simple_triplet_matrix}}.
}
\keyword{IO}
