\name{grouped}
\alias{grouped}
\alias{Grouped}
\alias{groupedData}
\alias{plot.groupedData}
\alias{hist.groupedData}
\title{Grouped data}
\description{ 
  For grouped data, the histogram is obtained by differentiating the
  ogive. 
}
\usage{
grouped(x, y = NULL)
\method{plot}{groupedData}(x, y = NULL, xlim = NULL, ylim = NULL, xlab =
"group boundaries", ylab = "F(x)", col = 1, ...)
\method{hist}{groupedData}(x, xlim = NULL, ylim = NULL, xlab = "group boundaries", ylab = "f(x)", ...)
}
\arguments{
  \item{x, y}{the \code{x} and \code{y} arguments provide the x and y coordinates for
          the plot.  Any reasonable way of defining the coordinates is
          acceptable.  See the function \code{xy.coords} for
	  details. Also, \code{x} may contains an object of class \code{data.frame}}
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
  for \eqn{c_{j-1} \leq x \leq c_j}.
  The formula for histogram is \deqn{f_n(x) = \frac{F_n(c_j) - F_n(c_{j-1})}{c_j - c_{j -
	1}}}{(F_n(c_j) - F_n(c_(j-1)))/(c_j - c_(j-1))}
  for \eqn{c_{j-1} \leq x < c_j}.
}
\value{
  Returns an object of \code{class} \code{"groupedData"}. Function \code{plot} can
     be used to draw the ogive and function \code{hist} can be used to
     construct histogram.
  An object of class \code{"groupedData"} is a list with the following
  components:
    \item{cj}{group boundaries.}
    \item{Fnt}{empirical distribution function evaluated at group boundaries.}
    \item{fnt}{ density function evaluated at group boundaries.}
}
\references{Klugman, Panjer & Willmot, Loss Models, Second
  Edition, Wiley, 2004.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
t<-grouped(data1)
plot(t)
hist(grouped(data1))
}
\keyword{hplot}
