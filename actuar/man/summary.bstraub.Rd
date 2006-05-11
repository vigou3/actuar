\name{summary.bstraub}
\alias{summary.bstraub}
\alias{print.summary.bstraub}
\title{Summary Method for BSTRAUB Objects}
\description{
  Summarizes a \code{\link{bstraub}} object and returns an object of
  class \code{summary.bstraub}. There is a \code{print} method for the latter.
}
\usage{
\method{summary}{bstraub}(object, \dots)

\method{print}{summary.bstraub}(x, \dots)
}  
\arguments{
  \item{x, object}{a \code{\link{bstraub}} object.}
  \item{\dots}{arguments passed to or from other methods.}
}
\details{
  This function is a method for the generic function summary() for class
  \code{bstraub}.It can be invoked by calling summary(x) for an object x of the
  appropriate class, or directly by calling summary.bstraub(x).
}
\value{
  Returns an object of class \code{summary.bstraub}. 
  A \code{print} method is defined to display such objects. It first
  returns the structure parameters estimated to compute the credibility
  premium. Four characteristics proper to each contract are subsequently
  displayed in a table.  
  \item{Collective premium}{collective premium estimator}
  \item{Average contract variance}{estimator of the within contract
    heterogeneity parameter}
  \item{Portefolio heterogeneity}{both unbiased and iterative estimators
    of the between contract heterogeneity parameter. The iterative
    estimator is NULL if argument \code{heterogeneity = "unbiased"} is
    passed to \code{bstraub}.}
  \item{Credibility constant}{estimator of the credibility constant,
    given by \eqn{K = \frac{s^2}/{a}}{K = s^2 / a}, where \eqn{s^2} is the average contract
    variance and \eqn{a} the portfolio heterogeneity.}
  {The following characteristics for each contract are returned in a
  table.}
  \item{Ind. premium}{individual experience weighted average}
  \item{Weight}{contract total weight, as used in credibility factors}
  \item{Cred. factor}{factor used with its complement in the linear
    expression of the credibility premium}
  \item{Cred. premium}{credibility premium for each contract}
  }
\author{
	Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and 
	Louis-Philippe Pouliot
}
\examples{
data(hachemeister)

## Credibility premiums calculated with the unbiased estimator
(x <- bstraub(hachemeister$claims,
              hachemeister$weights,
              heterogeneity ="unbiased"))
summary(x)
}

\keyword{file}

