\name{simpf}
\alias{simpf}
\title{Simulation of a portfolio of data}
\description{
  \code{simpf} simulates a portfolio of data for insurance
  applications. Both frequency and severity distributions can have an
  unknown risk parameter --- that is, they can each be mixtures of
  models.
}
\usage{
simpf(contracts, years, model.freq, model.sev, weights)
}
\arguments{
  \item{contracts}{the number of contracts in the porfolio}
  \item{years}{the number of years of experience for each contract}
  \item{model.freq}{named list containing the frequency model (see details
    below); if NULL, only claim amounts are simulated}
  \item{model.sev}{named list containing the severity model (see details
    below); if NULL, only claim numbers are simulated}
  \item{weights}{a matrix of weights (one per contract and per year) to
    be used in the simulation of frequencies}
}
\details{
  The function allows for compound models for both frequency and
  severity of losses. The compounding (or risk) parameter is called
  \code{Lambda} in the frequency model and \code{Theta} in the severity
  model.
  
  Distribution assumptions are specified using the base name of random
  number generation functions, e.g. "pois" for the Poisson distribution
  or "lnorm" for the lognormal.
  
  \code{model.freq} and \code{model.sev} are NULL or named lists
  composed of:
  \describe{
    \item{\code{dist1}}{base name of the distribution for a simple model,
      or of the conditional distribution for a compound model.}
    \item{\code{par1}}{named list of the parameters of \code{dist1} as
      they are defined in \code{rdist1}. If needed, the compounding parameter
      is identified by an unevaluated expression in \code{Lambda} and
      \code{weights} for \code{model.freq}, or \code{Theta} for
      \code{model.sev}.}
    \item{\code{dist2}}{base name of the compounding distribution, if
      any.}
    \item{\code{par2}}{named list of the parameters of \code{dist2} as
      they are defined in "rdist2".}
  }
}
\value{
  A list with two components:
    \item{\code{data}}{a two dimension (\code{contracts} rows and
      \code{years} columns) list where each element is a vector of
      losses, or a matrix if each element has length 1}
    \item{\code{weights}}{the matrix of weights given in argument}
}
\references{Goulet, V. (2006), \emph{Credibility for severity
    revisited}, North American Actuarial Journal}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Sebastien Auclair
}
\seealso{\code{\link{rearrangepf}}}
\examples{
## Portfolio where both frequency and severity models are compounded.
modelfreq <- list(dist1 = "pois",
                  par1 = list(lambda = quote(Lambda * weights)),
                  dist2 = "gamma",
                  par2 = c(shape = 2, rate = 1))
modelsev<-list(dist1 = "lnorm",
               par1 = list(meanlog = quote(Theta), sdlog = 1),
               dist2 = "norm",
               par2 = c(mean = 5, sd = 1))
data(hachemeister)
weights <- hachemeister$weights/mean(hachemeister$weights)
simpf(5, 12, modelfreq, modelsev, weights)

## Portfolio where the frequency model is compounded, but not the
## severity model.
modelsev <- list(dist1 = "lnorm", par1 = list(meanlog = 7, sd = 1))  
simpf(5, 12, modelfreq, modelsev)

## Portofolio with a severity model only and a user function for the
## simulation of observations.
rpareto <- function(n, alpha, lambda)
    lambda * (runif(n)^(-1/alpha) - 1)
modelsev <- list(dist1 = "pareto", par1 = list(alpha = 3, lambda = 8000))
simpf(5, 12, model.freq = NULL, modelsev)
}
\keyword{file}