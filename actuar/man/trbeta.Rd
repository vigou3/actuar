\name{Trbeta}
\alias{trbeta}
\alias{print.trbeta}
\alias{dtrbeta}
\alias{ptrbeta}
\alias{qtrbeta}
\alias{rtrbeta}
\alias{mtrbeta}
\alias{levtrbeta}
\alias{pearson6}
\alias{print.pearson6}
\alias{dpearson6}
\alias{ppearson6}
\alias{qpearson6}
\alias{rpearson6}
\alias{mpearson6}
\alias{levpearson6}
\title{The Transformed Beta Distribution}
\description{
  Density, distribution function, quantile function, random generation,
  raw moments, and limited moments for the Transformed Beta distribution with parameters \code{shape1},
  \code{scale}, \code{shape2}, and \code{shape3}.
}
\usage{
  dtrbeta(x, shape1, scale, shape2, shape3, log = FALSE)
  ptrbeta(q, shape1, scale, shape2, shape3, lower.tail = TRUE,
          log.p = FALSE)
  qtrbeta(p, shape1, scale, shape2, shape3, lower.tail = TRUE,
          log.p = FALSE)
  rtrbeta(n, shape1, scale, shape2, shape3)
  mtrbeta(k, shape1, scale, shape2, shape3, log = FALSE)
  levtrbeta(d, shape1, scale, shape2, shape3, order = 1, log = FALSE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{d}{limit of the loss variable.}
  \item{shape1, scale, shape2, shape3}{parameters. Must be strictly positive.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities p are
    returned as log(p).}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are P[X
    <= x], otherwise, Pr[X > x].}
  \item{k, order}{order of the moment.}
}
\details{ 
  The Transformed Beta distribution with parameters \code{shape1} = a,
  \code{scale} = b, \code{shape2} = c, and \code{shape3} = d has density:
  \deqn{p(x) = \frac{\Gamma(a + d)}{\Gamma(a)\Gamma(d)}
    \frac{c(\frac{x}{b})^{(c d)}}{x(1 + (\frac{x}{b})^c)^{(a + d)}}}{p(x)
    = (gamma(a + d) / gamma(a)*gamma(d))((c (x / b)^(a c)) / (x
    (1 + (x / b)^c)^(a + d)))} for \eqn{x > 0, a > 0, b > 0, c > 0} and
  \eqn{d > 0}. Here Gamma(a) is the function implemented by R's
  \code{gamma()} and defined in its help.

  Let Y be a Beta random variable with parameters
  \code{shape1} and \code{shape3}. Then, \deqn{scale\left(\frac{Y}{1 - Y}\right)^{\left(\frac{1}{shape2}\right)}}{
    scale(Y / (1 - Y))^(1 / shape2)} has the
  Transformed Beta probability density function with parameters
  \code{shape1}, \code{scale}, \code{shape2}, and \code{shape3}.

  The Transformed Beta distribution defines a family of distributions
  with the following special cases:
  \itemize{
    \item \code{shape3} == 1 is a Burr distribution,
    \item \code{shape3} == \code{shape2} == 1 is a Pareto distribution,
    \item \code{shape3} == \code{shape1} == 1 is a Loglogistic
    distribution,
    \item \code{shape3} == 1 and \code{shape2} == \code{shape1} is a
    Paralogistic distribution,
    \item \code{shape2} == 1 is a Generalized Pareto distribution,
    \item \code{shape2} == \code{shape1} == 1 is an Inverse Pareto
    distribution,
    \item \code{shape1} == 1 is an Inverse Burr distribution, and
    \item \code{shape1} == 1 and \code{shape2} == \code{shape3} is an
    Inverse Paralogistic distribution.
  }
}
\value{
\code{dtrbeta} gives the density, \code{ptrbeta} gives the distribution function,
\code{qtrbeta} gives the quantile function, \code{rtrbeta} generates random
deviates, \code{mtrbeta} gives the kth raw moment, and \code{levtrbeta} gives
the kth moment of the limited loss variable.

Invalid arguments will return warning.
}
\note{
  Also known as the generalized beta of the second kind and Pearson Type VI.
}
\references{
  Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
  Edition, Wiley, 2004.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
exp(dtrbeta(2, 2, 3, 4, 5, log = TRUE))
p <- (1:10)/10
ptrbeta(qtrbeta(p, 2, 3, 4, 5), 2, 3, 4, 5)
qpearson6(0.3, 2, 3, 4, 5, lower.tail = FALSE)
mtrbeta(2, 1, 2, 3, 4) - mtrbeta(1, 1, 2, 3, 4) ^ 2
levtrbeta(10, 1, 2, 3, 4, 2)
}
\keyword{distribution}
