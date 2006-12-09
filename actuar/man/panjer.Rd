\name{panjer}
\alias{panjer}
\title{Panjer's recursion formula}
\description{
  Panjer recursion formula to compute the approximate aggregate claim
  amount distribution of a portfolio. 
}
\usage{
panjer(fx, model.freq, p0 = NULL, x.scale = 1, TOL = 1e-8, echo = FALSE)
}

\arguments{
  \item{fx}{a vector of the (discretized) claim amount distribution;
    first element \strong{must} be \eqn{f_X(0)}.}
  \item{model.freq}{named list containing the frequency model (see details
    below).}
  \item{x.scale}{real value of a claim amount of 1 (monetary unit).}
  \item{p0}{arbitrary amount of probability at zero given to the
    frequency distribution. If not \code{NULL}, creates a zero-modified
    or zero-truncated distribution.}
  \item{TOL}{stop recursion when cumulative distribution function is
    less than \code{TOL} away from 1.}
  \item{echo}{print the cumulative distribution function as it is computed.}
}
\details{
  The frequency distribution is a member of the \eqn{(a, b, 0)} family
  of discrete distributions if \code{p0} is \code{NULL} and a member of
  the \eqn{(a, b, 1)} family if \code{p0} is specified. The logarithmic
  distribution is a limiting case of the negative binomial distribution
  where the size parameter is equal to 0.

  \code{model.freq} is a named list composed of:
  \describe{
    \item{\code{dist}}{one of
      \code{"binom"},
      \code{"geom"},
      \code{"nbinom"},
      \code{"pois"} or
      \code{"logarithmic"};}
    \item{\code{par}}{named list of the parameters of \code{dist} as
      they are defined in \code{ddist}.}
  }
}
\value{
  An object of class \code{"aggregateDist"}. See
  \code{\link{aggregateDist}} for details on such objects.
}
\seealso{
  \code{\link{discretize}},
  \code{\link{aggregateDist}}
}
\references{
  Klugman, S. A., Panjer, H. H. and Willmot, G. E. (2004),
  \emph{Loss Models, From Data to Decisions, Second Edition}, Wiley.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  Sébastien Auclair and Louis-Philippe Pouliot
}
\examples{
### (a, b, 0) family member (binomial)
model.freq <- list(dist = "binom", par = list(size = 3, prob = 0.5))
panjer(fx = rep(0.5, 2), model.freq = model.freq, echo = TRUE)

### (a, b, 1) family member (extended truncated negative binomial)
### [First part of Example 6.18 of Klugman et al. (2004)]
model.freq <- list(dist = "nbinom", par = list(size = 0.2, prob = 0.25))  
Fs <- panjer(fx = c(0.3, 0.5, 0.2), model.freq = model.freq, p0 = 0,
             x.scale = 10, TOL = 1E-5)
Fs(c(0, 10, 20, 30, 40))
}
\keyword{distribution}
\keyword{models}
