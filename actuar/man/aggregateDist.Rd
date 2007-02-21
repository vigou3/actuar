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
aggregateDist(method = c("recursive", "convolution", "normal", "npower", "simulation"),
              model.freq = NULL, model.sev = NULL, p0 = NULL, x.scale = 1, 
              moments, nb.simul, \dots, TOL = 1e-06, echo = FALSE)

\method{print}{aggregateDist}(x, \dots)

\method{plot}{aggregateDist}(x, xlim, \dots)

\method{summary}{aggregateDist}(object, \dots)

\method{mean}{aggregateDist}(x, \dots)
}
\arguments{
  \item{method}{method to be used}
  \item{model.freq}{for \code{"recursive"} method: a character string
    giving the name of a distribution in the \eqn{(a, b, 0)} or \eqn{(a,
      b, 1)} families of distributions. For \code{"convolution"} method:
    a vector of claim number probabilities. For \code{"simulation"}
    method: a frequency simulation model (see \code{\link{simpf}} for
    details) or \code{NULL}. Ignored with \code{normal} and
    \code{npower} methods.}
  \item{model.sev}{for \code{"recursive"} and \code{"convolution"}
    methods: a vector of claim amount probabilities. For
    \code{"simulation"} method: a severity simulation model (see
    \code{\link{simpf}} for details) or \code{NULL}. Ignored with
    \code{normal} and \code{npower} methods.}
  \item{p0}{arbitrary probability at zero for the frequency
    distribution. Creates a zero-modified or zero-truncated
    distribution if not \code{NULL}. Used only with \code{"recursive"}
    method.}
  \item{x.scale}{value of an amount of 1 in the severity model (monetary
    unit). Used only with \code{"recursive"} and \code{"convolution"}
    methods.}
  \item{moments}{vector of the true moments of the aggregate claim
    amount distribution; required only by the \code{"normal"}
    or \code{"npower"} methods.}
  \item{nb.simul}{number of simulations for the \code{"simulation"} method.} 
  \item{\dots}{parameters of the frequency distribution for the
    \code{"recursive"} method, or additionnal plotting options in the
    case of \code{plot.aggregateDist}.}
  \item{TOL}{the recursion in the \code{"recursive"} method stops when the
    cumulative distribution function is less than \code{TOL} away from 1.}
  \item{echo}{echo the recursions to screen in the \code{"recursive"} method.}
  \item{x, object}{an object of class \code{"aggregateDist"}.}
  \item{xlim}{the x limits (x1, x2) of the plot.} 
}
\details{
  \code{aggregateDist} essentially checks the arguments and then calls
  functions creating \code{"aggregateDist"} objects. Each method returns
  a function to compute the aggregate claim amount cumulative
  distribution function (cdf).

  The \code{"recursive"} method computes the cdf using the Panjer
  algorithm; the \code{"convolution"} method using convolutions; the
  \code{"normal"} method using a normal approximation; the
  \code{"npower"} method using the Normal Power 2 approximation; the
  \code{"simulation"} method using simulations. More details follow.
}
\section{Recursive method}{
  The frequency distribution is a member of the \eqn{(a, b, 0)} family
  of discrete distributions if \code{p0} is \code{NULL} and a member of
  the \eqn{(a, b, 1)} family if \code{p0} is specified. The logarithmic
  distribution is a limiting case of the negative binomial distribution
  where the size parameter is equal to 0.

  \code{model.freq} must be one of \code{"binomial"},
  \code{"geometric"}, \code{"negative binomial"}, \code{"poisson"} or
  \code{"logarithmic"} (these can abbreviated). The parameters of the
  frequency distribution must be specified using names identical to the
  arguments of functions \code{\link{dbinom}}, \code{\link{dgeom}},
  \code{\link{dnbinom}}, \code{\link{dpois}} and \code{\link{dnbinom}},
  respecively.

  \code{model.sev} is a vector of the (discretized) claim amount
  distribution; the first element \strong{must} be \eqn{f_X(0)}{fx(0)}.
}
\section{Convolution method}{
  The cumulative distribution function (cdf) \eqn{F_S(x)}{Fs(x)} of the
  aggregate claim amount of a portfolio in the collective risk model is
  \deqn{F_S(x) = \sum_{n = 0}^{\infty} F_X^{*n}(x) p_n,}{%
    Fs(x) = sum(n; Fx^*n(x) * pn)}
  for \eqn{x = 0, 1, \dots}; \eqn{p_n = \Pr[N = n]}{pn = Pr[N = n]} is
  the frequency probability mass function and
  \eqn{F_X^{*n}(x)}{Fx^*n(x)} is the cdf of the \eqn{n}th convolution of
  the (discrete) claim amount random variable.

  \code{model.freq} is vector \eqn{p_n}{pn} of the number of claims
  probabilities; the first element \strong{must} be \eqn{Pr[N = 0]}.

  \code{model.sev} is vector \eqn{f_X(x)}{fx(x)} of the (discretized)
  claim amount distribution; the first element \strong{must} be
  \eqn{f_X(0)}{fx(0)}.
}
\section{Normal and Normal Power 2 methods}
{
  The Normal approximation of a cumulative distribution function (cdf)
  \eqn{F(x)} with mean \eqn{\mu}{m} and standard deviation
  \eqn{\sigma}{s} is
  \deqn{F(x) \approx \Phi\left( \frac{x - \mu}{\sigma} \right).}{%
    F(x) ~= pnorm((x - m)/s).}

  The Normal Power 2 approximation of a cumulative distribution function (cdf)
  \eqn{F(x)} with mean \eqn{\mu}{m}, standard deviation \eqn{\sigma}{s}
  and skewness \eqn{\gamma}{g} is
  \deqn{F(x) \approx \Phi \left(%
    -\frac{3}{\gamma} + \sqrt{\frac{9}{\gamma^2} + 1 %
      + \frac{6}{\gamma} \frac{x - \mu}{\sigma}} \right).}{%
    F(x) ~= pnorm(-3/g + sqrt(9/g^2 + 1 + (6/g) * (x - m)/s)).}
  This formula is valid only for the right-hand tail of the distribution
  and skewness should not exceed unity.
}
\section{Simulation method}{
  This methods returns the empirical distribution function of a sample
  of size code{nb.simul} of the aggregate claim amount distribution
  specified by \code{model.freq} and
  \code{model.sev}. \code{\link{simpf}} is used for the simulation of
  claim amounts, hence both the frequency and severity models can be
  mixtures of distributions.
}
\value{
  A function of class \code{"aggregateDist"}, inheriting from the
  \code{"function"} class when using normal and Normal Power
  approximations and additionally inheriting from the \code{"ecdf"} and
  \code{"stepfun"} classes when other methods are used.  
}
\seealso{
  \code{"\link{discretize}"} to discretize a severity distribution;
  \code{\link{simpf}}.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2004),
  \emph{Loss Models, From Data to Decisions, Second Edition}, Wiley.

  Daykin, C.D., Pentikäinen, T. and Pesonen, M. (1994), \emph{Practical 
    Risk Theory for Actuaries}, Chapman & Hall.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Louis-Philippe Pouliot
}
\examples{
## Convolution method (example 6.6 of Klugman et al. (2004))
fx <- c(0, 0.15, 0.2, 0.25, 0.125, 0.075, 0.05, 0.05, 0.05, 0.025, 0.025)
pn <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.15, 0.06, 0.03, 0.01)
Fs <- aggregateDist("convolution", model.freq = pn, model.sev = fx,
                    x.scale = 25)
summary(Fs)
c(Fs(0), diff(Fs(25 * 0:21))) # prob. mass function
plot(Fs)

## Recursive method
model.freq <- list(dist = "pois", par = list(lambda = 3))
Fs <- aggregateDist("recursive", model.freq = "poisson", model.sev = fx,
                    lambda = 3, x.scale = 25)
plot(Fs)

## Normal Power approximation
Fs <- aggregateDist("npower", moments = c(200, 200, 0.5))
Fs(210)

## Simulation method
model.freq <- expression(rpois(3))
model.sev <- expression(rgamma(100, 2))
Fs <- aggregateDist("simulation", model.freq, model.sev, nb.simul = 1000)
mean(Fs)
plot(Fs)
}
\keyword{distribution}
\keyword{models}
