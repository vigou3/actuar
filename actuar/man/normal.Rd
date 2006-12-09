\name{normal}
\alias{normal}
\alias{np2}
\title{
  Normal and Normal Power Approximations
}
\description{
  Normal and Normal Power 2 approximations of a cumulative distribution
  function.
}
\usage{
normal(mean, var)

np2(mean, var, skewness)
}
\arguments{
  \item{mean}{true mean of the distribution}
  \item{var}{true variance of the distribution}
  \item{skewness}{true skewness of the distribution}
}
\details{
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
\value{
  An object of class \code{"aggregateDist"}. See
  \code{\link{aggregateDist}} for details on such objects.
}
\references{
  Daykin, C.D., Pentikäinen, T. and Pesonen, M. (1994), \emph{Practical 
    Risk Theory for Actuaries}, Chapman & Hall.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Louis-Philippe Pouliot}
\examples{
normal(mean = 200, var = 500)(c(220, 250, 300))
np2(mean = 200, var = 500, skewness = 0.25)(c(220, 250, 300))
}
\keyword{distribution}
\keyword{models}
