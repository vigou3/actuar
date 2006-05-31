\name{summary.simpf}
\alias{summary.simpf}
\alias{print.summary.simpf}
\title{Summary Method for SIMPF Objects}
\description{
  Summarizes a \code{\link{simpf}} object by displaying the individual
  claim amounts for the desired years and contracts. The function returns an object of
  class \code{summary.simpf}. There is a \code{print} method for the latter.
}
\usage{
\method{summary}{simpf}(object, contracts = seq(length = dim(object$data)[1]),
       years = seq(length = dim(object$data)[2]), \dots)

\method{print}{summary.simpf}(x, contracts = x$contracts,
       years = x$years, \dots)
}  
\arguments{
  \item{x, object}{a \code{\link{simpf}} object.}
  \item{contracts}{a vector of contracts for which claim amounts should be
    displayed. All contracts are displayed by default.}
  \item{years}{a vector of years for which claim amounts should be
    displayed. All years observed are displayed by default}
  \item{\dots}{arguments passed to or from other methods.}
}
\details{
  This function is a method for the generic function summary() for class
  \code{simpf}.It can be invoked by calling summary(x) for an object x of the
  appropriate class, or directly by calling summary.simpf(x).
}
\value{
  Returns an object of class \code{summary.simpf}. 
  A \code{print} method is defined to display such objects. The latter
  displays a list of the chosen contracts (passed on to the argument
  \code{contracts}), in which each element contains itself a list of
  claim amounts for each year passed on to the argument \code{years}. 
  }
\author{
	Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and 
	Louis-Philippe Pouliot
}
\examples{
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
(x <- simpf(5, 12, model.freq = modelfreq, modelsev))
summary(x, contracts = 1:2, years = 10:12)
}

\keyword{internal}

