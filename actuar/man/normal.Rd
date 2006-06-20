\name{normal}
\alias{normal}
\alias{np2}
\title{Normal and Normal Power Approximation of the Total Amount of Claims Ditribution }
\description{
  Approximates the distribution of the total amount of claims. 
}
\usage{
normal(x, mean, var)

np2(x, mean, var, skewness)

}

\arguments{
  \item{x}{ value at which the function is evaluated}
  \item{mean}{ true mean of the distribution}
  \item{var}{ true variance of the distribution}
  \item{skewness}{true skewness of the distribution}
}
\details{
  \code{normal} computes the very common normal standardization, having:

  \deqn{F(X) \approx N(v(X))}{F(X) ~= N(v(X))}

  where \emph{N} is the normal d.f. and

  \deqn{v(X) = \frac{(X - \mu_x)}{\sigma_X}}{
    v(X) = (X - mu_x)/sigma_x}

  standardizes the variable.

  \code{np2} approximation has formula:

  \deqn{F(x) \approx N(v(X)) =
    N\left(-\frac{3}{\gamma_X}+\sqrt{\frac{9}{\gamma_X^2}+1+\frac{6}{\gamma_X}*\frac{X-\mu_X}{\sigma_X}}\right)}{
      F(X) ~= N(v(X)) = N(-3/gamma_x + sqrt(9/gamma_x^2 + 1 + 6/gamma_x
      *(X-mu_x)/sigma_x))}

  This formula is valid only for the right-hand tail of the distribution and skewness should not exceed unity.
}
\value{
 An object of class \code{aggregateDist}. For details on the elements of
 such objects, refer to the more general function \code{\link{aggregateDist}}. 
}
\seealso{
  \code{\link{aggregateDist}} \code{\link{panjer}} \code{\link{exact}}
  \code{\link{simS}}
  }
\references{ Daykin, C.D. and Pentikäinen, T. and Pesonen, M. (1994), \emph{Practical
Risk Theory for Actuaries}, Chapman & Hall, Sections 4.2.2 - 4.2.4}

\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Louis-Philippe Pouliot}

\examples{
x <- seq(0, qnorm(1-1E-5, mean = 200, sd = sqrt(500)), by = 0.001)

distribution <- normal(x, mean = 200, var = 500)
plot(distribution)


distribution <- np2(x, mean = 200, var = 500, skewness = 0.25)
plot(distribution)
}
\keyword{distribution}
\keyword{models}
