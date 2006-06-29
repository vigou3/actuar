\name{ogive}
\alias{ogive}
\alias{Ogive}
\alias{plot.ogive}
\alias{print.ogive}
\title{Empirical distribution function for grouped data}
\description{
  Compute or plot an empirical distribution function for grouped data.
  The ogive is obtained by connecting the values of the empirical distribution function at the group boundaries with
  straight lines. 
}
\usage{
ogive(x, y = NULL)
\method{plot}{ogive}(x, y = NULL, xlim = NULL, ylim = NULL, xlab =
"boundaries", ylab = "F(x)", col = 1, ...)
\method{print}{ogive}(x, digits= getOption("digits") - 2, ...)
}
\arguments{
  \item{x, y}{grouped data with finite group boundaries.}
   \item{xlim}{the x limits (x1,x2) of the plot.  Note that \code{x1 > x2} is
          allowed and leads to a "reversed axis".}
  \item{ylim}{the y limits of the plot.}
  \item{xlab}{a label for the x axis.}
  \item{ylab}{a label for y axis.}
  \item{col}{the colors for lines and points.}
  \item{digits}{the desired number of digits after the decimal point.}
  \item{\dots}{potential further arguments required by generic.}
}
\details{
  This function uses \code{approxfun} to evaluate empirical distribution
  function.
  The theorical formula for ogive is \deqn{F_n(x) = \frac{c_j - x}{c_j - c_{j -
	1}}F_n(c_{j-1}) + \frac{x - c_{j-1}}{c_j -
      c_{j-1}}F_n(c_j)}{(c_j - x)/(c_j)- c_(j-1) F_n(c_(j-1)) + (x -
    c_(j-1))/(c_j - c_(j-1)) F_n(c_j)}
  for \eqn{c_{j-1} \leq x \leq c_j}.
}
\value{
  Returns an object of \code{class} \code{ogive}. Function
  \code{plot} can be used to draw the ogive. An object of \code{class}
  \code{ogive} is a function.
}
\references{Klugman, Panjer & Willmot, Loss Models, Second
  Edition, Wiley, 2004.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
cj <- c(0, 10, 20, 30, 40)
nj <- c(4, 9, 0, 12)
data <- grouped(cj, nj)
fnt <- ogive(data)
plot(fnt)
fnt(7)
}
\keyword{hplot}
