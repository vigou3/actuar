\name{simS}
\alias{simS}
\title{Simulation of a Total Amount of Claims Distribution}
\description{
  Approximate the total amount of claims distribution by simulation.
}
\usage{
simS(n, model.freq, model.sev)
}

\arguments{
  \item{n}{ number of observations to be simulated}
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
 An empirical distribution function of class \code{aggregateDist}. For details on
 such objects, refer to the more general function \code{\link{aggregateDist}}. 
}
\seealso{
  \code{\link{aggregateDist}} \code{\link{panjer}} \code{\link{exact}}
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