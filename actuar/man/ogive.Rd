\name{ogive}
\alias{ogive}
\alias{Ogive}
\alias{plot.ogive}
\alias{lines.ogive}
\title{Empirical distribution function for grouped data}
\description{
  For grouped data, the ogive is obtained by connecting the values of
  the empirical distribution function at the group boundaries with
  straight lines. 
}
\usage{
ogive(x, lines = FALSE, xlim = NULL, ylim = NULL, xlab = "group boundaries", ylab = "F(x)", col = 1)

\method{plot}{ogive}(x, y = NULL, xlim = NULL, ylim = NULL, xlab =
"group boundaries", ylab = "F(x)", col = 1, ...)
\method{lines}{ogive}(x, xlim = NULL, ylim = NULL, col = 1, ...)
}
\arguments{
  \item{x, y}{grouped data with finite group boundaries.}
  \item{lines}{logical; if \code{FALSE}, generic function \code{plot} is
  called, if \code{TRUE}, generic function \code{lines} is called.}
  \item{xlim}{the x limits (x1,x2) of the plot.  Note that \code{x1 > x2} is
          allowed and leads to a "reversed axis".}
  \item{ylim}{the y limits of the plot.}
  \item{xlab}{a label for the x axis}
  \item{ylab}{a label for y axis.}
  \item{col}{the colors for lines and points.}
  \item{\dots}{potential further arguments required by generic.}
}
\details{
  The formula for ogive is \deqn{F_n(x) = \frac{c_j - x}{c_j - c_{j -
	1}}F_n(c_{j-1}) + \frac{x - c_{j-1}}{c_j -
      c_{j-1}}F_n(c_j)}{(c_j - x)/(c_j)- c_(j-1) F_n(c_(j-1)) + (x -
    c_(j-1))/(c_j - c_(j-1)) F_n(c_j)}
  for \eqn{c_{j-1} < x < c_j}.
}
\value{
  Returns an ogive for grouped data.
  
}
\references{Klugman, Panjer & Willmot, Loss Models, Second
  Edition, Wiley, 2004.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
ogive(data1, lines = FALSE, xlim = c(0,1000000))
ogive(data2, lines = TRUE, col = 3)
}
\keyword{hplot}
