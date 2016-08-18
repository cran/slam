\name{options}
\alias{slam_options}
\title{Options for the `slam' package}
\description{Function for getting and setting options for the
  \pkg{slam} package.}
\usage{
slam_options(option, value)
}
\arguments{
  \item{option}{character string indicating the option to get
    or set (see details). If missing, all options are returned as a list.}
  \item{value}{Value to be set. If omitted, the current value is returned.}
}
\details{
  Currently, the following options are available:

  \describe{
    \item{\code{"max_dense"}:}{numeric specifying the maximum \code{length}
      of dense vectors (default: \code{2^24}).
    }
  }
}
\seealso{
  \code{\link{simple_sparse_array}}
}
\examples{
## save defaults
.slam_options <- slam_options()
.slam_options

slam_options("max_dense", 2^25)
slam_options("max_dense")

## reset
slam_options("max_dense", .slam_options$max_dense)
}
\keyword{math}
