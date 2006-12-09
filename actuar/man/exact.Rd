\name{exact}
\alias{exact}
\title{Exact Calculation of the Aggregate Claim Amount Distribution}
\description{
  Compute the aggregate claim amount cumulative distribution function of
  a portfolio over a period by convolutions.
}
\usage{
exact(fx, pn, x.scale = 1)
}
\arguments{
  \item{fx}{a vector of the (discretized) claim amount distribution;
    first element \strong{must} be \eqn{f_X(0)}.}
  \item{pn}{a vector of the number of claims probabilities; first
    element \strong{must} be \eqn{Pr[N = 0]}.}
  \item{x.scale}{real value of a claim amount of 1 (monetary unit).}  
}
\details{
  The cumulative distribution function (cdf) \eqn{F_S(x)}{Fs(x)} of the
  aggregate claim amount of a portfolio in the collective risk model is
  \deqn{F_S(x) = \sum_{n = 0}^{\infty} F_X^{*n}(x) p_n,}{%
    Fs(x) = sum(n; Fx^*n(x) * pn)}
  for \eqn{x = 0, 1, \dots}; \eqn{p_n = \Pr[N = n]}{pn = Pr[N = n]} is
  the frequency probability mass function and
  \eqn{F_X^{*n}(x)}{Fx^*n(x)} is the cdf of the \eqn{n}th convolution of
  the (discrete) claim amount random variable.
}
\value{
  An object of class \code{"aggregateDist"}. See
  \code{\link{aggregateDist}} for details on such objects.
}
\seealso{
  \code{\link{discretize}},
  \code{\link{aggregateDist}}
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
fx <- c(0, 0.15, 0.2, 0.25, 0.125, 0.075, 0.05, 0.05, 0.05, 0.025, 0.025)
pn <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.15, 0.06, 0.03, 0.01)
Fs <- exact(fx, pn, x.scale = 25)
c(Fs(0), diff(Fs(25 * 0:21)))
}
\keyword{distribution}
\keyword{models}
