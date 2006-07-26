\name{mde}
\alias{Mde}
\alias{mde}
\title{Minimum distance estimation}
\description{
  Minimum distance estimators for individual and grouped data.
}
\usage{
  mde(x, cdf, start, methods  = c("CvM", "chi-square", "LAS"), w = NULL, ...)
}
\arguments{
  \item{x}{individual or grouped data.}
  \item{cdf}{Either a character string or a function returning a
    cumulative distribution evaluated at its first argument.}
  \item{start}{A named list giving the parameters to be optimized with
          initial values }
  \item{methods}{Either \code{CvM} for Cramer von Mises' method,
    \code{chi-square} for modified chi-square's method, and \code{LAS} for ???.}
  \item{w}{weights}
  \item{\dots}{potential further arguments required by generic.}
}
\details{ 
  \code{CvM}'s method minimizes distance between theorical and empirical
  cumulative distribution function, modified \code{chi-square}'s method
  minimizes chi-square statistic, and \code{LAS}'s method minimizes
  distance between ???.
}
\value{
 A list with three components,

 estimate: the parameter estimates,

 sd: the estimated standard errors, and

 loglik: the log-likelihood.

 Invalid arguments will return warning.
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
x <- rpareto(10000, 8, 10)
d <- list(shape = 2, scale = 2)
mde(x, ppareto, d, "CvM")
t <- grouped(x = c(2, 4, 6, 10, 20), y = c(1, 4, 2, 1))
mde(t, pgamma, d, "chi-square")
}
\keyword{models}
