\name{TransformedBeta}
\alias{TransformedBeta}
\alias{dtrbeta}
\alias{ptrbeta}
\alias{qtrbeta}
\alias{rtrbeta}
\alias{mtrbeta}
\alias{levtrbeta}
\alias{Pearson6}
\alias{dpearson6}
\alias{ppearson6}
\alias{qpearson6}
\alias{rpearson6}
\alias{mpearson6}
\alias{levpearson6}
\title{The Transformed Beta Distribution}
\description{
  Density, distribution function, quantile function, random generation,
  raw moments and limited moments for the Transformed Beta distribution
  with parameters \code{shape1}, \code{shape2}, \code{shape3} and
  \code{scale}.
}
\usage{
  dtrbeta(x, shape1, shape2, shape3, rate = 1, scale = 1/rate,
          log = FALSE)
  ptrbeta(q, shape1, shape2, shape3, rate = 1, scale = 1/rate,
          lower.tail = TRUE, log.p = FALSE)
  qtrbeta(p, shape1, shape2, shape3, rate = 1, scale = 1/rate,
          lower.tail = TRUE, log.p = FALSE)
  rtrbeta(n, shape1, shape2, shape3, rate = 1, scale = 1/rate)
  mtrbeta(order, shape1, shape2, shape3, rate = 1, scale = 1/rate,
          log = FALSE)
  levtrbeta(limit, shape1, shape2, shape3, rate = 1, scale = 1/rate,
            order = 1, log = FALSE)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length is
    taken to be the number required.}
  \item{shape1, shape2, shape3, scale}{parameters. Must be strictly
    positive.}
  \item{rate}{an alternative way to specify the scale.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities p are
    returned as \eqn{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise, \eqn{P[X > x]}.}
  \item{order}{order of the moment.}
  \item{limit}{limit of the loss variable.}
}
\details{ 
  The Transformed Beta distribution with parameters \code{shape1} \eqn{=
    \alpha}{= a}, \code{shape2} \eqn{= \gamma}{= b}, \code{shape3}
  \eqn{= \tau}{= c} and \code{scale} \eqn{= \theta}{= s}, has
  density:
  \deqn{f(x) = \frac{\Gamma(\alpha + \tau)}{\Gamma(\alpha)\Gamma(\tau)}
    \frac{\gamma (x/\theta)^{(\gamma \tau)}}{%
      x (1 + (x/\theta)^\gamma)^{(\alpha + \tau)}}}{%
    f(x) = (Gamma(a + c)/Gamma(a) * Gamma(c)) ((b (x/s)^(bc))/%
    (x (1 + (x/s)^b)^(a + c)))}
  for \eqn{x > 0}, \eqn{\alpha > 0}{a > 0}, \eqn{\gamma > 0}{b > 0},
  \eqn{\tau > 0}{c > 0} and \eqn{\theta > 0}{s > 0}.
  (Here \eqn{\Gamma(\alpha)}{Gamma(a)} ma(a) is the function implemented
  by R's \code{\link{gamma}()} and defined in its help.)

  The Transformed Beta is the distribution of the random variable
  \deqn{\theta \left(\frac{X}{1 - X}\right)^{(1/\gamma)},}{%
    s (X/(1 - X))^(1/b),}
  where \eqn{X} has a Beta distribution with parameters \eqn{\alpha}{a}
  and \eqn{\tau}{c}.

  The Transformed Beta distribution defines a family of distributions
  with the following special cases:
  \itemize{
    \item A \link[=dburr]{Burr} distribution when \code{shape3 == 1};
    \item A \link[=dllogis]{Loglogistic} distribution when \code{shape1
      == shape3 == 1};
    \item A \link[=dparalogis]{Paralogistic} distribution when
    \code{shape3 == 1} and \code{shape2 == shape1};
    \item A \link[=dgenpareto]{Generalized Pareto} distribution when
    \code{shape2 == 1};
    \item A \link[=dpareto]{Pareto} distribution when \code{shape2 ==
      shape3 == 1};
    \item An \link[=dinvburr]{Inverse Burr} distribution when
    \code{shape1 == 1};
    \item An \link[=dinvpareto]{Inverse Pareto} distribution when
    \code{shape2 == shape1 == 1};
    \item An \link[=dinvparalogis]{Inverse Paralogistic} distribution
    when \code{shape1 == 1} and \code{shape3 == shape2}.
  }
}
\value{
  \code{dtrbeta} gives the density,
  \code{ptrbeta} gives the distribution function,
  \code{qtrbeta} gives the quantile function,
  \code{rtrbeta} generates random deviates,
  \code{mtrbeta} gives the \eqn{k}th raw moment, and
  \code{levtrbeta} gives the \eqn{k}th moment of the limited loss
  variable.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  Distribution also known as the Generalized Beta of the Second
  Kind and Pearson Type VI.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (1998),
  \emph{Loss Models, From Data to Decisions}, Wiley.
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
