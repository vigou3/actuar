\name{simS}
\alias{simS}
\title{Simulation of a Total Amount of Claims Distribution}
\description{
  \code{simS} simulates a distribution of the total amount of claims.
}
\usage{
simS(x, model.freq, model.sev)
}

\arguments{
  \item{x}{ number of observations to be simulated}
  \item{model.freq}{ named list containing the frequency model (see
    details below)}
  \item{model.sev}{ named list containing the frequency model (see
    details below)}
}
\details{
  \code{model.freq} and \code{model.sev} are named lists composed of:

  \item{dist}{base name of the distribution model}
  \item{par}{named vector of the parameters of \code{dist} as they are
    defined in \code{pdist}} 
}
\value{
 An object of class \code{AggregateDist}. For details on the elements of
 such objects, refer to the more general function \code{\link{ADgen}}. 
}
\seealso{
  \code{\link{ADgen}} \code{\link{panjer}} \code{\link{exact}}
  \code{\link{normal}} \code{\link{np2}}
  }
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Louis-Philippe Pouliot}

\examples{
model.freq <- list(dist = "pois", par = c(lambda = 4))
model.sev <- list(dist = "norm", par = c(mean = 10, sd = 1.5))
simS(1e+05, model.freq, model.sev)

}
\keyword{distribution}
\keyword{models}