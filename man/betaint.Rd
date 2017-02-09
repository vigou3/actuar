\name{betaint}
\alias{betaint}
\title{The \dQuote{Beta Integral}}
\description{
  The \dQuote{beta integral} which is just a multiple of the non
  regularized incomplete beta function. This function merely
  provides an R interface to the C level routine. It is not exported
  by the package.
}
\usage{
betaint(x, a, b)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{a, b}{parameters. See Details for admissible values.}
}
\details{
  Function \code{betaint} computes the \dQuote{beta integral}
  \deqn{
    B(a, b; x) = \Gamma(a + b) \int_0^x t^{a-1} (1-t)^{b-1} dt}{%
    B(a, b; x) = Gamma(a + b) int_0^x t^(a-1) (1-t)^(b-1) dt}
  for \eqn{a > 0}, \eqn{b \neq -1, -2, \ldots}{b != -1, -2, \ldots} and
  \eqn{0 < x < 1}.
  (Here \eqn{\Gamma(\alpha)}{Gamma(a)} is the function implemented
  by \R's \code{\link{gamma}()} and defined in its help.)
  When \eqn{b > 0},
  \deqn{
    B(a, b; x) = \Gamma(a) \Gamma(b) I_x(a, b),}
  where \eqn{I_x(a, b)} is \code{pbeta(x, a, b)}. When \eqn{b < 0},
  \eqn{b \neq -1, -2, \ldots}{b != -1, -2, \ldots}, and \eqn{a > 1 +
  [-b]}{a > 1 + floor(-b)},
  \deqn{%
    \begin{array}{rcl}
      B(a, b; x)
      &=& \displaystyle
          -\Gamma(a + b) \left[ \frac{x^{a-1} (1-x)^b}{b}
          + \frac{(a-1) x^{a-2} (1-x)^{b+1}}{b (b+1)} \right. \\
      & & \displaystyle\left.
	  + \cdots + \frac{(a-1) \cdots (a-r) x^{a-r-1}
	  (1-x)^{b+r}}{b (b+1) \cdots (b+r)} \right] \\
      & &  \displaystyle
  	  + \frac{(a-1) \cdots (a-r-1)}{b (b+1) \cdots (b+r)}
	  \Gamma(a-r-1) \\
      & & \times \Gamma(b+r+1) I_x(a-r-1, b+r+1),
      \end{array}}{%
    B(a, b; x)
    = -Gamma(a+b) \{(x^(a-1) (1-x)^b)/b
        + [(a-1) x^(a-2) (1-x)^(b+1)]/[b(b+1)]
        + \dots
	+ [(a-1)\dots(a-r) x^(a-r-1) (1-x)^(b+r)]/[b(b+1)\dots(b+r)]\}
      + [(a-1)\dots(a-r-1)]/[b(b+1)\dots(b+r)] Gamma(a-r-1)
        * Gamma(b+r+1) I_x(a-r-1, b+r+1),}
  where \eqn{r = [-b]}{r = floor(-b)}.

  This function is used (at the C level) to compute the
  limited expected value for distributions of the transformed beta
  family; see, for example, \code{\link{levtrbeta}}.
}
\value{
  The value of the integral.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  The need for this function in the package is well explained in the
  introduction of Appendix A of Klugman et al. (2012). See also chapter
  6 and 15 of Abramowitz and Stegun (1972) for definitions and relations
  to the hypergeometric series.
}
\references{
  Abramowitz, M. and Stegun, I. A. (1972), \emph{Handbook of Mathematical
  Functions}, Dover.

  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2012),
  \emph{Loss Models, From Data to Decisions, Fourth Edition}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca}
}
\examples{
x <- 0.3
a <- 7

## case with b > 0
b <- 2
actuar:::betaint(x, a, b)
gamma(a) * gamma(b) * pbeta(x, a, b)    # same

## case with b < 0
b <- -2.2
r <- floor(-b)        # r = 2
actuar:::betaint(x, a, b)

## "manual" calculation
s <- (x^(a-1) * (1-x)^b)/b +
    ((a-1) * x^(a-2) * (1-x)^(b+1))/(b * (b+1)) +
    ((a-1) * (a-2) * x^(a-3) * (1-x)^(b+2))/(b * (b+1) * (b+2))
-gamma(a+b) * s +
    (a-1)*(a-2)*(a-3) * gamma(a-r-1)/(b*(b+1)*(b+2)) *
    gamma(b+r+1)*pbeta(x, a-r-1, b+r+1)
}
\keyword{math}
\keyword{distribution}
