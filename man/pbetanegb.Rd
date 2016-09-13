\name{pbetanegb}
\alias{pbetanegb}
\title{The Beta Cumulative Distribution for Negative Second Parameter}
\description{
  (Scaled) Cumulative distribution function for the Beta distribution
  with negative, non integer, second parameter.
}
\usage{
pbetanegb(x, a, b)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{a, b}{parameters. See Details for admissible values.}
}
\details{
  The function computes the value of \eqn{\Gamma(a) \Gamma(b) B_x(a, b)}
  when \eqn{b < 0}, \eqn{b \neq -1, -2, \ldots}{b != -1, -2, \ldots} and
  \eqn{a > 1 + [-b]} as follows:
  \deqn{%
    \begin{array}{rcl}
      \Gamma(a) \Gamma(b) B_x(a, b)
      &=& \displaystyle
          -\Gamma(a + b) \left[ \frac{x^{a-1} (1-x)^b}{b}
          + \frac{(a-1) x^{a-2} (1-x)^{b+1}}{b (b+1)} \right. \\
      & & \displaystyle\left.
	  + \cdots + \frac{(a-1) \cdots (a-r) x^{a-r-1}
	  (1-x)^{b+r}}{b (b+1) \cdots (b+r)} \right] \\
      & &  \displaystyle
  	  + \frac{(a-1) \cdots (a-r-1)}{b (b+1) \cdots (b+r)}
	  \Gamma(a-r-1) \\
      & & \times \Gamma(b+r+1) B_x(a-r-1, b+r+1),
      \end{array}}{%
    Gamma(a) Gamma(b) B_x(a, b)
    = -Gamma(a + b) \{(x^(a-1) (1-x)^b)/b
        + [(a-1) x^(a-2) (1-x)^(b+1)]/[b(b+1)]
        + \dots
	+ [(a-1)\dots(a-r) x^(a-r-1) (1-x)^(b+r)]/[b(b+1)\dots(b+r)]\}
      + [(a-1)\dots(a-r-1)]/[b(b+1)\dots(b+r)] Gamma(a-r-1)
        * Gamma(b+r+1) B_x(a-r-1, b+r+1),}
  where \eqn{r = [-b]}{r = floor(-b)} and \eqn{B_x(a, b)} is the
  \code{\link{pbeta}} function.

  The function is mostly mostly used at the C level to compute the
  limited expected value for distributions of the transformed beta
  family; see, for example, \code{\link{levtrbeta}}.
}
\value{
  The (scaled) distribution function.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
a <- 7
b <- -2.2
r <- floor(-b)        # r = 2
x <- 0.3

pbetanegb(x, a, b)

## "manual" calculation
s <- (x^(a-1) * (1-x)^b)/b +
    ((a-1) * x^(a-2) * (1-x)^(b+1))/(b * (b+1)) +
    ((a-1) * (a-2) * x^(a-3) * (1-x)^(b+2))/(b * (b+1) * (b+2))
-gamma(a+b) * s +
    (a-1)*(a-2)*(a-3) * gamma(a-r-1)/(b*(b+1)*(b+2)) *
    gamma(b+r+1)*pbeta(x, a-r-1, b+r+1)
}
\keyword{distribution}
