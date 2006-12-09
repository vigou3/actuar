\name{aggregateDist}
\alias{aggregateDist}
\alias{print.aggregateDist}
\alias{plot.aggregateDist}
\alias{summary.aggregateDist}
\alias{mean.aggregateDist}
\title{Aggregate Claim Amount Distribution}
\description{
  Compute the aggregate claim amount cumulative distribution function of
  a portfolio over a period using one of five methods.
}
\usage{
aggregateDist(method = c("recursive", "exact", "normal", "np2", "simulation"),
              model.sev, model.freq, p0 = NULL, x.scale = 1, moments, nb.simul,
              TOL = 1e-06, echo = FALSE, \dots)

\method{print}{aggregateDist}(x, \dots)

\method{plot}{aggregateDist}(x, xlim, \dots)

\method{summary}{aggregateDist}(object, \dots)

\method{mean}{aggregateDist}(x, \dots)
}
\arguments{
  \item{method}{method to be used}
  \item{model.sev}{severity model (see details below) for
    \code{"simulation"} method or a vector of claim amount probabilities for
    \code{"recursive"} and \code{"exact"} methods. Ignored with
    \code{normal} and \code{np2} methods.}
  \item{model.freq}{frequency model (see details below) for
    \code{"simulation"} and \code{"recursive"} methods or a vector of
    the claim number probabilities for \code{"exact"} method. Ignored
    with \code{normal} and \code{np2} methods.}
  \item{p0}{arbitrary probability at zero for the frequency
    distribution. Creates a zero-modified or zero-truncated
    distribution. Used only with \code{"recursive"} method.}
  \item{x.scale}{value of an amount of 1 in the severity model (monetary
    unit). Used only with \code{"recursive"} and \code{"exact"} methods.}
  \item{moments}{vector of the moments required by the \code{"normal"}
    or \code{"np2"} methods.}
  \item{nb.simul}{number of simulations in the \code{"simulation"} method.} 
  \item{TOL}{the recursion in the \code{recursive} method stops when the
    cumulative distribution function is less than \code{TOL} away from 1.}
  \item{echo}{echo the recursions to screen in the \code{recursive} method.}
  \item{x, object}{an object of class \code{"aggregateDist"}.}
  \item{xlim}{the x limits (x1, x2) of the plot.} 
  \item{\dots}{arguments potentially required by additionnal methods and
    additionnal plotting options in the case of \code{plot.aggregateDist}}   	  
}
\details{
  \code{aggregateDist} essentially calls functions creating
  \code{"aggregateDist"} objects. Each method returns a function to
  compute the aggregate claim amount cumulative distribution function (cdf).

  The \code{"recursive"} method computes the cdf using the Panjer
  algorithm; the \code{"exact"} method using convolutions; the
  \code{"normal"} method using a normal approximation; the \code{"np2"}
  method using the normal power approximation; the \code{"simulation"}
  method using simulations.
}
\value{
  A function of class \code{"aggregateDist"}, inheriting from the
  \code{"function"} class when using Normal and Normal Power
  approximations and additionally inheriting from the \code{"ecdf"} and
  \code{"stepfun"} classes when other methods are used.  
}
\seealso{
  \code{"\link{panjer}"} for the \code{"recursive"}  method,
  \code{"\link{exact}"}  for the \code{"exact"}      method,
  \code{"\link{normal}"} for the \code{"normal"}     method,
  \code{"\link{np2}"}    for the \code{"np2"}        method,
  \code{"\link{simS}"}   for the \code{"simulation"} method,
  \code{"\link{discretize}"}.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2004),
  \emph{Loss Models, From Data to Decisions, Second Edition}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Louis-Philippe Pouliot
}
\examples{
## Exact method (example 6.6 of Klugman et al. (2004))
fx <- c(0, 0.15, 0.2, 0.25, 0.125, 0.075, 0.05, 0.05, 0.05, 0.025, 0.025)
pn <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.15, 0.06, 0.03, 0.01)
Fs <- aggregateDist("exact", model.sev = fx, model.freq = pn, x.scale = 25)
summary(Fs)
c(Fs(0), diff(Fs(25 * 0:21))) # prob. mass function
plot(Fs)

## Recursive method
#  model.sev <- c(0.1, 0.1, 0.2, 0.3, 0.2, 0.1)
model.freq <- list(dist = "pois", par = list(lambda = 3))
Fs <- aggregateDist("recursive", model.sev = fx, model.freq = model.freq,
                    x.scale = 25)
plot(Fs)

## Normal Power approximation
Fs <- aggregateDist("np2", moments = c(200, 200, 0.5))
Fs(210)

## Simulation method
model.freq <- list(dist1 = "pois", par1 = list(lambda = 3))
model.sev <- list(dist1 = "gamma", par1 = list(shape = 100, rate = 2))
Fs <- aggregateDist("simulation", model.sev = model.sev,
                     model.freq = model.freq, nb.simul = 10000)
mean(Fs)
plot(Fs)
}
\keyword{distribution}
\keyword{models}
