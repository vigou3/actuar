\name{histog}
\alias{histog}
\alias{Histog}
\alias{plot.histog}
\title{Empirical density function for grouped data}
\description{
  For grouped data, the histogram is obtained by differentiating the
  ogive (see ogive). 
}
\usage{
histog(x, xlim = NULL, ylim = NULL, xlab = "group boundaries", ylab =
"f(x)")

\method{plot}{histog}(x, y = NULL, xlim = NULL, ylim = NULL, xlab =
"group boundaries", ylab = "f(x)", ...)
}
\arguments{
  \item{x, y}{grouped data with finite group boundaries.}
  \item{xlim}{the x limits (x1,x2) of the plot.  Note that \code{x1 > x2} is
          allowed and leads to a "reversed axis".}
  \item{ylim}{the y limits of the plot.}
  \item{xlab}{a label for the x axis}
  \item{ylab}{a label for y axis.}
  \item{\dots}{potential further arguments required by generic.}
}
\details{
  The formula for histogram is \deqn{f_n(x) = \frac{F_n(c_j) - F_n(c_{j-1})}{c_j - c_{j -
	1}}}{(F_n(c_j) - F_n(c_(j-1)))/(c_j - c_(j-1))}
  for \eqn{c_{j-1} <= x < c_j}.
}
\value{
  Returns a histogram for grouped data.
  
}
\references{Klugman, Panjer & Willmot, Loss Models, Second
  Edition, Wiley, 2004.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
histog(data1, xlim = c(0,1000000))
}
\keyword{hplot}
