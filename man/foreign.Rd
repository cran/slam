\name{foreign}
\alias{read_stm_CLUTO}
\alias{write_stm_CLUTO}
\alias{read_stm_MC}
\alias{write_stm_MC}
\title{Read and Write Sparse Matrix Format Files}
\description{
  Read and write CLUTO sparse matrix format files, or the CCS format
  variant employed by the MC toolkit.
}
\usage{
read_stm_CLUTO(file)
write_stm_CLUTO(x, file)
read_stm_MC(file, scalingtype = NULL)
write_stm_MC(x, file)
}
\arguments{
  \item{file}{a character string with the name of the file to read or
    write.}
  \item{x}{a matrix object.}
  \item{scalingtype}{a character string specifying the type of scaling
    to be used, or \code{NULL} (default), in which case the scaling will
    be inferred from the names of the files with non-zero entries found
    (see \bold{Details}).
  }
}
\details{
  Documentation for CLUTO including its sparse matrix format is
  available from \url{https://www-users.cs.umn.edu/~karypis/cluto/}.

  \code{read_stm_CLUTO} reads CLUTO sparse matrices, returning a
  \link[slam:matrix]{simple triplet matrix}.

  \code{write_stm_CLUTO} writes CLUTO sparse matrices.
  Argument \code{x} must be coercible to a simple triplet matrix via
  \code{\link[slam:matrix]{as.simple_triplet_matrix}}.

  MC is a toolkit for creating vector models from text documents (see 
  \url{https://www.cs.utexas.edu/users/dml/software/mc/}).  It employs a
  variant of Compressed Column Storage (CCS) sparse matrix format,
  writing data into several files with suitable names: e.g., a file with
  \file{_dim} appended to the base file name stores the matrix
  dimensions.  The non-zero entries are stored in a file the name of
  which indicates the scaling type used: e.g., \file{_tfx_nz} indicates
  scaling by term frequency (\samp{t}), inverse document frequency
  (\samp{f}) and no normalization (\samp{x}).  See \file{README} in the
  MC sources for more information.

  \code{read_stm_MC} reads such sparse matrix information with argument
  \code{file} giving the path with the base file name, and returns a
  \link[slam:matrix]{simple triplet matrix}.

  \code{write_stm_MC} writes matrices in MC CCS sparse matrix format.
  Argument \code{x} must be coercible to a simple triplet matrix via
  \code{\link[slam:matrix]{as.simple_triplet_matrix}}.
}
\keyword{IO}
