\name{ogive}
\alias{ogive}
\alias{Ogive}
\title{Empirical distribution function for grouped data}
\description{
  For grouped data, the ogive is obtained by connecting the values of
  the empirical distribution function at the group boundaries with
  straight lines. 
}
\usage{
ogive(x, y = NULL, xlim = NULL, ylim = NULL, xlab = "group boundaries", ylab = "F(x)", col = 1)
}
\arguments{
  \item{x, y}{grouped data with finite group boundaries.}
  \item{xlim}{the x limits (x1,x2) of the plot.  Note that \code{x1 > x2} is
          allowed and leads to a "reversed axis".}
  \item{ylim}{the y limits of the plot.}
  \item{xlab}{a label for the x axis}
  \item{ylab}{a label for y axis.}
  \item{col}{the colors for lines and points.}
}
\details{
  This function calls function \code{grouped} and
  method \code{plot.groupedData}. 
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
ogive(data1, xlim = c(0,1000000), col = 3)
}
\keyword{hplot}
