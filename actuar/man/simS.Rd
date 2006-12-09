\name{simS}
\alias{simS}
\title{Simulation of the Aggregate Claim Amount Dstribution}
\description{
  Approximation of the aggregate claim amount cumulative distribution
  function of a portfolio over a period by simulation.
}
\usage{simS(n, model.freq, model.sev)}
\arguments{
  \item{n}{size of the sample to simulate.}
  \item{model.freq}{named list containing the frequency model (see
    \code{\link{simpf}} for details).}
  \item{model.sev}{named list containing the frequency model (see
    \code{\link{simpf}} for details).}
}
\details{
  \code{simS} uses \code{\link{simpf}} for the simulation of claim
  amounts. 
}
\value{
  An object of class \code{"aggregateDist"}. See
  \code{\link{aggregateDist}} for details on such objects.
}
\seealso{
  \code{\link{simpf}},
  \code{\link{aggregateDist}}
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Louis-Philippe Pouliot
}
\examples{
model.freq <- list(dist1 = "pois", par1 = c(lambda = 4))
model.sev <- list(dist1 = "norm", par1 = c(mean = 10, sd = 1.5))
simS(1e+04, model.freq, model.sev)
}
\keyword{distribution}
\keyword{models}
