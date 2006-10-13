\name{GeneralizedPareto}
\alias{GeneralizedPareto}
\alias{dgenpareto}
\alias{pgenpareto}
\alias{qgenpareto}
\alias{rgenpareto}
\alias{mgenpareto}
\alias{levgenpareto}
\title{The Generalized Pareto Distribution}
\description{
  Density, distribution function, quantile function, random generation,
  raw moments and limited moments for the Generalized Pareto
  distribution with parameters \code{shape1}, \code{shape2} and
  \code{scale}.
}
\usage{
  dgenpareto(x, shape1, shape2, rate = 1, scale = 1/rate,
             log = FALSE)
  pgenpareto(q, shape1, shape2, rate = 1, scale = 1/rate,
             lower.tail = TRUE, log.p = FALSE)
  qgenpareto(p, shape1, shape2, rate = 1, scale = 1/rate,
             lower.tail = TRUE, log.p = FALSE)
  rgenpareto(n, shape1, shape2, rate = 1, scale = 1/rate)
  mgenpareto(order, shape1, shape2, rate = 1, scale = 1/rate)
  levgenpareto(limit, shape1, shape2, rate = 1, scale = 1/rate,
               order = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{shape1, shape2, scale}{parameters. Must be strictly positive.}
  \item{rate}{an alternative way to specify the scale.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities
    \eqn{p} are returned as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
}
\details{
  The Generalized Pareto distribution with parameters \code{shape1}
  \eqn{= \alpha}{= a}, \code{shape2} \eqn{= \tau}{= b} and \code{scale}
  \eqn{= \theta}{= s} has density:
  \deqn{f(x) = \frac{\Gamma(\alpha + \tau)}{\Gamma(\alpha)\Gamma(\tau)}
    \frac{\theta^\alpha x^{\tau - 1}}{%
      (x + \theta)^{\alpha + \tau}}}{%
    f(x) = Gamma(a + b)/(Gamma(a) * Gamma(b)) (s^a x^(b - 1))/%
    (x + s)^(a + b))}
  for \eqn{x > 0}, \eqn{\alpha > 0}{a > 0}, \eqn{\tau > 0}{b > 0} and
  \eqn{\theta > 0}{s > 0}.
  (Here \eqn{\Gamma(\alpha)}{Gamma(a)} is the function implemented
  by R's \code{\link{gamma}()} and defined in its help.)

  The Generalized Pareto is the distribution of the random variable
  \deqn{\theta \left(\frac{X}{1 - X}\right),}{%
    s (X/(1 - X)),}
  where \eqn{X} has a Beta distribution with parameters \eqn{\alpha}{a}
  and \eqn{\tau}{b}.

  The Generalized Pareto distribution has the following special cases:
  \itemize{
    \item A \link[=dpareto]{Pareto} distribution when \code{shape2 ==
      1};
    \item An \link[=dinvpareto]{Inverse Pareto} distribution when
    \code{shape1 == 1}.
  }
}
\value{
  \code{dgenpareto} gives the density,
  \code{pgenpareto} gives the distribution function,
  \code{qgenpareto} gives the quantile function,
  \code{rgenpareto} generates random deviates,
  \code{mgenpareto} gives the \eqn{k}th raw moment, and
  \code{levgenpareto} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  Distribution also known as the Beta of the Second Kind.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2004),
  \emph{Loss Models, From Data to Decisions, Second Edition}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
exp(dgenpareto(3, 3, 4, 4, log = TRUE))
p <- (1:10)/10
pgenpareto(qgenpareto(p, 3, 3, 1), 3, 3, 1)
qgenpareto(.3, 3, 4, 4, lower.tail = FALSE)
mgenpareto(1, 3, 2, 1) ^ 2
levgenpareto(10, 3, 3, 3, 2)
}
\keyword{distribution}
