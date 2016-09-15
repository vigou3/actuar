\name{SpecialIntegrals}
\alias{SpecialIntegrals}
\alias{expint}
\alias{betaint}
\alias{gammaint}
\title{The Exponential Integral and Two Other Special Integrals}
\description{
  The exponential integral \eqn{E_1(x)} and two integrals related to the
  incomplete beta function and the incomplete gamma function. These
  functions merely provide R interfaces to the C level functions. They
  are not exported by the package.
}
\usage{
expint(x)
betaint(x, a, b)
gammaint(x, a)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{a, b}{parameters. See Details for admissible values.}
}
\details{
  Function \code{expint} computes the exponential integral
  \deqn{
    G(0; x) = E_1(x) = \int_x^\infty \frac{e^{-t}}{t}\, dt}{%
    G(0; x) = E_1(x) = int_x^Inf \exp(-t)/t dt}
  for \eqn{x \neq}{x != 0} a real number.

  Function \code{betaint} computes the integral
  \deqn{
    B(a, b; x) = \Gamma(a + b) \int_0^x t^(a-1) (1-t)^(b-1) dt}{%
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
      & & \times \Gamma(b+r+1) B_x(a-r-1, b+r+1),
      \end{array}}{%
    Gamma(a) Gamma(b) B_x(a, b)
    = -Gamma(a + b) \{(x^(a-1) (1-x)^b)/b
        + [(a-1) x^(a-2) (1-x)^(b+1)]/[b(b+1)]
        + \dots
	+ [(a-1)\dots(a-r) x^(a-r-1) (1-x)^(b+r)]/[b(b+1)\dots(b+r)]\}
      + [(a-1)\dots(a-r-1)]/[b(b+1)\dots(b+r)] Gamma(a-r-1)
        * Gamma(b+r+1) I_x(a-r-1, b+r+1),}
  where \eqn{r = [-b]}{r = floor(-b)}.

  Function \code{gammaint} computes the integral
  \deqn{
    G(\alpha; x) = \int_x^\infty t^{\alpha - 1} e^{-t}\, dt}{%
    G(a; x) = int_x^Inf t^(a - 1) exp(-t) dt}
  for \eqn{\alpha}{a} real and \eqn{x > 0}. When \eqn{\alpha > 0}{a >
  0},
  \deqn{
    G(\alpha; x) = \Gamma(a) [1 - \Gamma(\alpha; x)],}{%
    G(a; x) = Gamma(a) [1 - Gamma(a; x)],}
  where \eqn{\Gamma(\alpha; x)}{Gamma(a; x)} is \code{pgamma(x, a)}.
  When \eqn{\alpha \le 0}{a <= 0},
  \deqn{
    G(\alpha; x) = -\frac{x^\alpha * e^{-x}}{\alpha}
      + \frac{1}{\alpha} G(\alpha + 1; x)}{%
    G(a; x) = -(x^a * exp(-x))/a + G(a + 1; x)/a}
  which can be repeated until \eqn{\alpha + k}{a + k} is a positive
  number. If \eqn{\alpha = 0, -1, -2, \ldots}{a = 0, -1, -2, \ldots}
  this requires \eqn{G(0; x)}, the exponential integral above.

  These function are mostly mostly used at the C level to compute the
  limited expected value for distributions of the transformed beta
  (\code{betaint}) and transformed gamma (\code{gammaint} and
  \code{expint}) families; see, for example, \code{\link{levtrbeta}} and
  \code{\link{levinvtrgamma}}.
}
\value{
  The value of the respective integral.

  Invalid arguments will result in return value \code{NaN}, with a warning.
}
\note{
  The C implementation of \code{expint} is based on code from the GNU
  Software Library \url{https://www.gnu.org/software/gsl/}.

  The rationale for these function in the package is explained in the
  introduction of Appendix A of Klugman et al. (2012). See also chapter
  5 and 6 of Abramowitz and Stegun (1972).
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
###
### Exponential integral
###
actuar:::expint(1.275)           # 0.1408099 in Abramowitz and Stegun
actuar:::expint(10)              # 4.15697e-6 in Abramowitz and Stegun

###
### Incomplete beta related integral
###
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

###
### Incomplete gamma related integral
###
x <- 1.2

## case with a > 0
a <- 1.5
actuar:::gammaint(x, a)
gamma(a) * pgamma(x, a, lower = FALSE)  # same

## case with a < 0 non integer
a <- -1.5
actuar:::gammaint(x, a)

## "manual" calculation
(((gamma(a+2)*pgamma(x, a+2, 1, lower = FALSE)
    -(x^(a+1) * exp(-x)))/(a+1))
    - (x^a * exp(-x)))/a

## case with a < 0 integer
a <- -2
actuar:::gammaint(x, a)

## "manual" calculation
(((actuar:::expint(x) -(x^(a+1) * exp(-x)))/(a+1))
    - (x^a * exp(-x)))/a
}
\keyword{distribution}
