\name{aggregateDist}
\alias{aggregateDist}
\alias{print.aggregateDist}
\alias{plot.aggregateDist}
\alias{mean.aggregateDist}
\title{Compute the Distribution of Claim Amounts }
\description{
  Use one of five methods to generate the
  aggregate claims distribution of a portfolio over one year given a frequency and a
  severity model or the true moments of the distribution.  
}
\usage{
aggregateDist(method = c("normal", "np2", "simulation", "recursive", "exact"),
      model.sev, model.freq, moments = c(mean = 0, var = 1, skewness = NULL), 
      x.scale = 1, n, p0, TOL = 1e-06, echo = FALSE, ...)

\method{print}{aggregateDist}(x, \dots)

\method{plot}{aggregateDist}(x, xlim, \dots)

\method{mean}{aggregateDist}(x, \dots)
}

\arguments{
  \item{method}{ method to be used}
  \item{model.freq}{ vector of the (discretized) claim amount
    distribution or a named list containing the frequency model (see details
    below). Ignored when using \code{normal} and \code{np2} methods.}
  \item{model.sev}{ vector of the number of claims probabilities or a
    named list containing the severity model (see details
    below). Ignored when using \code{normal} and \code{np2} methods.}
  \item{moments}{ vector containing the moments required by the
    method chosen. Ignored unless \code{normal} or \code{np2} method is called.}
  \item{x.scale}{ value of an amount of 1 in the severity model.}
  \item{n}{ number simulations when using the simulation method} 
  \item{p0}{ arbitrary amount of probability at zero given to the
    frequency distribution. It creates zero-modified or
    zero-truncated distributions. Ignored unless \code{recursive}
    method is called.}
  \item{TOL}{ used to stop recursion when cumulative probability function is
    less than TOL away from 1 in the \code{recursive} method and passed
    on to a quantile function to limit the discretization of
    distributions when parameterized models (named lists) are given in argument for
    the \code{recursive} and \code{exact} methods.}
  \item{echo}{print the cumulative distribution of the total amount of
    claims as it is computed when approximating by recursion.}
  \item{x}{ an \code{aggregateDist} object}
  \item{xlim}{the x limits (x1,x2) of the plot.} 
  \item{...}{ arguments potentially required by additionnal methods and
    additionnal plotting options in the case of \code{plot.aggregateDist}}   	  
}
\details{

  \code{aggregateDist} essentiallly calls functions creating \code{aggregateDist}
  objects. \code{\link{normal}}, \code{\link{np2}}, and \code{\link{exact}} methods call
  functions of the same name.  \code{simulation} and \code{recursive}
  methods respectively call functions \code{simS} and \code{panjer}.

  \code{normal} and \code{np2} methods only require the \code{moments} argument. The
  latter must contain the correctly ordered moments required by
  \code{normal} or \code{np2}, that is the true mean and the true
  variance in the first case, and additionally the true skewness in the second.
  Other arguments will be ignored. \code{moments} will only be read if one of those two methods is called. 

  The simulation method directly calls simS(), with
  a number of observations equal to \code{n}. \code{model.freq} and
  \code{model.sev} must necessarily be named lists representing
  parameterized distributions (see \code{\link{simS}}). Only those three
  arguments will be read.

  The exact method gives the possibility to directly send
  user-vectorized severity and frequency models in argument as required
  by \code{exact}, but also has the ability to discretize models sent as
  named lists (as in the \code{simulation} method).

  The \code{recursive} method gives the same latitude, but only for the
  severity model. That is, \code{model.freq} must be given accordingly to the
  original usage of \code{\link{panjer}}. 
}
\value{
A function of class \code{"aggregateDist"}, inheriting from the
\code{"function"} class when using Normal and Normal Power
approximations and additionally inheriting from the \code{"ecdf"} and
\code{"stepfun"} classes when other methods are used.  
  
  
}
\seealso{
  \code{\link{panjer}} \code{\link{exact}}
  \code{\link{simS}} \code{\link{np2}} \code{\link{normal}}
  }

\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Louis-Philippe Pouliot}

\examples{
## Normal Power approximation
  res <- aggregateDist("np2", moments = c(200, 200, 0.5))
  res(210)

## Simulation method
  model.freq <- list(dist = "pois", par = list(lambda = 4))
  model.sev <- list(dist = "gamma", par = list(shape = 100, rate = 2))
  res <- aggregateDist("simulation", model.sev = model.sev,
                       model.freq = model.freq, n = 10000)
  mean(res)
  res(250)

## Approximation by the recursive method
  # with a parameterized distribution severity model

  model.sev <- list(dist = "gamma", par = list(shape = 100, rate = 2))
  (res <- aggregateDist("recursive", model.sev, model.freq,
                       x.scale = 1))

  # with a vectorized severity model

  model.sev <- c(0.1, 0.1, 0.2, 0.3, 0.2, 0.1)
  (res <- aggregateDist("recursive", model.sev, model.freq,
                       x.scale = 100))
  plot(res)

## Exact calculation

  model.freq <- c(0.1, 0.1, 0.2, 0.3, 0.2, 0.1)
  (res <- aggregateDist("exact", model.sev = model.sev, model.freq = model.freq, x.scale = 100))
  plot(res)
  res(1000)

}
\keyword{distribution}
\keyword{models}