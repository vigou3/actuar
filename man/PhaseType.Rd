\name{PhaseType}
\alias{dphtype}
\alias{pphtype}
\alias{rphtype}
\alias{mphtype}
\alias{mgfphtype}
\title{The Phase-type Distribution}
\description{
  Density, distribution function, random generation, raw moments and
  moment generating function for the (continuous) Phase-type
  distribution with parameters \code{prob} and \code{rates}.
}
\usage{
dphtype(x, prob, rates, log = FALSE)
pphtype(q, prob, rates, lower.tail = TRUE, log.p = FALSE)
rphtype(n, prob, rates)
mphtype(order, prob, rates)
mgfphtype(x, prob, rates, log = FALSE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{prob}{vector of initial probabilities for each of the transient
    states of of the underlying Markov chain. The initial probability of
    the absorbing state is \code{1 - sum(prob)}.} 
  \item{rates}{square matrix of the rates of transition among the states
    of the underlying Markov chain.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment.}
}
\details{
  The phase-type distribution with parameters \code{prob} \eqn{= \pi}{=
    pi} and \code{rates} \eqn{= T} has density:%
  \deqn{f(x) = \pi e^{T x} t}{f(x) = pi \%*\% exp(T * x) \%*\% t}% 
  for \eqn{x > 0} and \eqn{f(0) = 1 - \pi e}{f(0) = 1 - pi \%*\% e},
  where \eqn{e} is a column vector with all components equal to one,
  \eqn{t = -T e}{t = -T \%*\% e} is the exit rates vector and
  \eqn{e^{Tx}}{exp(T * x)} denotes the matrix exponential of \eqn{Tx}{T
    * x}. The parameters must satisfy \eqn{\pi e \leq 1}{pi \%*\% e <=
    1}, \eqn{T_{ii} < 0}{T[i, i] < 0}, \eqn{T_{ij} \geq 0}{T[i, j] >= 0}
  and \eqn{Te \leq 0}{T \%*\% e <= 0}.

  If the value of \code{prob} do no sum to one, the distribution is
  mixed with an atom of \code{1 - sum(prob)} at \code{x = 0}.

  The matrix exponential of a matrix \eqn{M} is defined as the
  Taylor series%
  \deqn{e^M = \sum_{n = 0}^{\infty} \frac{M^n}{n!}.}{%
    exp(M) = sum(n = 0:Inf; (M^n)/(n!)).}

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]} and the moment generating function is
  \eqn{E[e^{xX}]}.
}
\value{
  \code{dphasetype} gives the density,
  \code{pphasetype} gives the distribution function,
  \code{rphasetype} generates random deviates,
  \code{mphasetype} gives the \eqn{k}th raw moment, and
  \code{mgfphasetype} gives the moment generating function in \code{x}.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}

\references{
  \url{http://en.wikipedia.org/wiki/Phase-type_distribution}
  
  Neuts, M. F. (1981), \emph{Generating random variates from a
    distribution of phase type}, WSC '81: Proceedings of the 13th
  conference on Winter simulation, IEEE Press.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Christophe Dutang
}
\examples{
## Erlang(3, 2) distribution
T <- cbind(c(-2, 0, 0), c(2, -2, 0), c(0, 2, -2))
pi <- c(1,0,0)
x <- 0:10

dphtype(x, pi, T)		# density
dgamma(x, 3, 2)			# same
pphtype(x, pi, T)		# cdf
pgamma(x, 3, 2)			# same

rphtype(10, pi, T)		# random values

mphtype(1, pi, T)		# expected value

curve(mgfphtype(x, pi, T), from = -10, to = 1)
}
\keyword{distribution}
