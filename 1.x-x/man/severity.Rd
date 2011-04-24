\name{severity}
\alias{severity}
\alias{severity.default}
\title{Manipulation of Individual Claim Amounts}
\description{
  \code{severity} is a generic function created to manipulate individual
  claim amounts. The function invokes particular \emph{methods} which
  depend on the \code{\link[base]{class}} of the first argument.
}
\usage{
severity(x, ...)

\method{severity}{default}(x, bycol = FALSE, drop = TRUE, \dots)
}
\arguments{
  \item{x}{an \R object.}
  \item{bycol}{logical; whether to \dQuote{unroll} horizontally
    (\code{FALSE}) or vertically (\code{TRUE})}
  \item{\dots}{further arguments to be passed to or from other methods.}
  \item{drop}{logical; if \code{TRUE}, the result is coerced to the
    lowest possible dimension.}
}
\details{
  Currently, the default method is equivalent to
  \code{\link{unroll}}. This is liable to change since the link between
  the name and the use of the function is rather weak.
}
\value{
  A vector or matrix.
}
\seealso{
  \code{\link{severity.portfolio}} for the original motivation of these
  functions.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Louis-Philippe Pouliot
}
\examples{
x <- list(c(1:3), c(1:8), c(1:4), c(1:3))
(mat <- matrix(x, 2, 2))
severity(mat)
severity(mat, bycol = TRUE)
}
\keyword{datagen}
\keyword{manip}
