\name{InverseGamma}
\alias{InverseGamma}
\alias{dinvgamma}
\alias{pinvgamma}
\alias{qinvgamma}
\alias{rinvgamma}
\alias{minvgamma}
\alias{levinvgamma}
\alias{mgfinvgamma}
\title{The Inverse Gamma Distribution}
\description{
  Density function, distribution function, quantile function, random generation,
  raw moments, and limited moments for the Inverse Gamma distribution
  with parameters \code{shape} and \code{scale}.
}
\usage{
dinvgamma(x, shape, rate = 1, scale = 1/rate, log = FALSE)
pinvgamma(q, shape, rate = 1, scale = 1/rate,
          lower.tail = TRUE, log.p = FALSE)
qinvgamma(p, shape, rate = 1, scale = 1/rate,
          lower.tail = TRUE, log.p = FALSE)
rinvgamma(n, shape, rate = 1, scale = 1/rate)
minvgamma(order, shape, rate = 1, scale = 1/rate)
levinvgamma(limit, shape, rate = 1, scale = 1/rate,
            order = 1)
mgfinvgamma(t, shape, rate =1, scale = 1/rate, log =FALSE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{shape, scale}{parameters. Must be strictly positive.}
  \item{rate}{an alternative way to specify the scale.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
  \item{t}{numeric vector.}
}
\details{
  The inverse gamma distribution with parameters \code{shape} \eqn{=
    \alpha}{= a} and \code{scale} \eqn{= \theta}{= s} has density:
  \deqn{f(x) = \frac{u^\alpha e^{-u}}{x \Gamma(\alpha)}, %
    \quad u = \theta/x}{%
    f(x) = u^a exp(-u)/(x Gamma(a)), u = s/x}
  for \eqn{x > 0}, \eqn{\alpha > 0}{a > 0} and \eqn{\theta > 0}{s > 0}.
  (Here \eqn{\Gamma(\alpha)}{Gamma(a)} is the function implemented
  by \R's \code{\link{gamma}()} and defined in its help.)

  The special case \code{shape == 1} is an
  \link[=dinvexp]{Inverse Exponential} distribution.

  The \eqn{k}th raw moment of the random variable \eqn{X} is
  \eqn{E[X^k]}{E[X^k]}, \eqn{k < \alpha}{k < shape}, and the \eqn{k}th
  limited moment at some limit \eqn{d} is \eqn{E[\min(X, d)^k]}{E[min(X,
  d)^k]}, all \eqn{k}.

  The moment generating function is given by \eqn{E[e^{tX}]}.
}
\value{
  \code{dinvgamma} gives the density,
  \code{pinvgamma} gives the distribution function,
  \code{qinvgamma} gives the quantile function,
  \code{rinvgamma} generates random deviates,
  \code{minvgamma} gives the \eqn{k}th raw moment,
  \code{levinvgamma} gives the \eqn{k}th moment of the limited loss
  variable, and
  \code{mgfinvgamma} gives the moment generating function in \code{t}.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  \code{levinvgamma} computes the limited expected value using
  \code{gammainc} from package \pkg{expint}.

  Also known as the Vinci distribution. See also Kleiber and
  Kotz (2003) for alternative names and parametrizations.
}
\references{
  Kleiber, C. and Kotz, S. (2003), \emph{Statistical Size Distributions
  in Economics and Actuarial Sciences}, Wiley.

  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
exp(dinvgamma(2, 3, 4, log = TRUE))
p <- (1:10)/10
pinvgamma(qinvgamma(p, 2, 3), 2, 3)
minvgamma(-1, 2, 2) ^ 2
levinvgamma(10, 2, 2, order = 1)
mgfinvgamma(-1, 3, 2)
}
\keyword{distribution}
