\name{unroll}
\alias{unroll}
\title{Display a Two-Dimension Version of a Matrix of Vectors}
\description{
  Displays all values of a matrix of vectors by \dQuote{unrolling} the
  object vertically or horizontally.
}
\usage{
unroll(x, bycol = FALSE, drop = TRUE)
}
\arguments{
  \item{x}{a list of vectors with a \code{\link[base]{dim}} attribute of
    length 0, 1 or 2.}
  \item{bycol}{logical; whether to unroll horizontally
    (\code{FALSE}) or vertically (\code{TRUE}).}
  \item{drop}{logical; if \code{TRUE}, the result is coerced to the
    lowest possible dimension.}
}
\details{
  \code{unroll} returns a matrix where elements of \code{x}
  are concatenated (\dQuote{unrolled}) by row (\code{bycol = FALSE}) or
  by column (\code{bycol = TRUE}). \code{NA} is used to make
  rows/columns of equal length.

  Vectors and one dimensional arrays are coerced to \strong{row}
  matrices.
}
\value{
  A vector or matrix.
}
\seealso{
  This function was originally written for use in
  \code{\link{severity.portfolio}}.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Louis-Philippe Pouliot
}
\examples{
x <- list(c(1:3), c(1:8), c(1:4), c(1:3))
(mat <- matrix(x, 2, 2))

unroll(mat)
unroll(mat, bycol = TRUE)

unroll(mat[1, ])
unroll(mat[1, ], drop = FALSE)
}
\keyword{manip}
