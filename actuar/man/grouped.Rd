\name{grouped}
\alias{grouped}
\alias{Grouped}
\alias{grouped.data}
\alias{hist.grouped.data}
\alias{print.grouped.data}
\title{Grouped data}
\description{
  Compute or plot a histogram for grouped data.
  For grouped data, the histogram is obtained by differentiating the
  ogive. 
}
\usage{
grouped(x, y = NULL, digits = 2)
\method{hist}{grouped.data}(x, y = NULL, main = "Histogram", xlim = NULL, ylim = NULL, xlab =
"boundaries", ylab = "f(x)", plot = TRUE, ...)
\method{print}{grouped.data}(x, ...)
}
\arguments{
  \item{x, y}{the \code{x} and \code{y} arguments provide the x and y coordinates for
    the plot. \code{Length(x)} must be equal to \code{length(y) + 1}. If
    \code{length(x)} = \code{length(y)}, \code{y[1]} is ignored.
    Also, \code{x} may contains an object of class \code{data.frame}}
  \item{main}{title of the histogram}
  \item{xlim}{the x limits (x1,x2) of the plot.  Note that \code{x1 > x2} is
          allowed and leads to a "reversed axis".}
  \item{ylim}{the y limits of the plot.}
  \item{xlab}{a label for the x axis}
  \item{ylab}{a label for y axis.}
  \item{plot}{logical.  If \code{TRUE} (default), a histogram is plotted,
          otherwise a list of \code{cj}, \code{nj} and \code{density} is returned.}
  \item{digits}{the desired number of digits after the decimal point.}
  \item{\dots}{potential further arguments required by generic.}
}
\details{
  The theorical formula for histogram is \deqn{f_n(x) = \frac{F_n(c_j) - F_n(c_{j-1})}{c_j - c_{j -
	1}}}{(F_n(c_j) - F_n(c_(j-1)))/(c_j - c_(j-1))}
  for \eqn{c_{j-1} \leq x < c_j}.
}
\value{
  Returns an object of \code{class} \code{"grouped.data"}. Function \code{hist} can be used to
  construct histogram and function \code{ogive} can be used to calculate the ogive.
  An object of class \code{"grouped.data"} is a list with the following
  components:
    \item{cj}{boundaries.}
    \item{nj}{}
}
\references{Klugman, Panjer & Willmot, Loss Models, Second
  Edition, Wiley, 2004.
}
\note{
  See also \code{\link{ogive}}. 

}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
cj <- c(0, 10, 20, 30, 40)
nj <- c(4, 12, 11, 2)
data <- grouped(cj, nj)
hist(data, plot = FALSE)
plot(ogive(data))
}
\keyword{hplot}
