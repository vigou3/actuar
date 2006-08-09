\name{distrmod}
\alias{distrmod}
\alias{Distrmod}
\title{Density and Cumulative Distribution Function for Modified Data}
\description{
  Compute density and cumulative distribution function for data with
  deductible \code{d} and limit \code{u}. Distribution can be created
  for per-loss variable and per-payment variable.
}
\usage{
distrmod(dist, param, d, u, density = TRUE, per.payment = TRUE, ...)
}
\arguments{
  \item{dist}{named list which contains two elements: name of
    distribution and names of parameters. The distribution must exist in
    R.}
   \item{param}{values of parameters.}
  \item{d}{deductible. It must be less than limit.}
  \item{u}{limit}
  \item{density}{logical; if \code{TRUE} (default) probability/density
    is returned.}
  \item{per.payment}{logical; if \code{TRUE} (default) per-payment
    function is returned, if \code{FALSE} per-loss function is returned.}
  \item{\dots}{potential further arguments required by generic.}
}
\details{
  Per-payment variable is defined as follow:
  \itemize{
    \item undefined for \eqn{X <= d},
    \item \eqn{(X - d)} for \eqn{d < X <= u}, and
    \item \eqn{(u - d)} for \eqn{X > u}.
  }
  Per-loss variable is define as follow:
  \itemize{
    \item \eqn{0} for \eqn{X <= d},
    \item \eqn{(X - d)} for \eqn{d < X <= u}, and
    \item \eqn{(u - d)} for \eqn{X > u}.
    }
}
\value{
  Returns an object of \code{class} \code{function}. 
}
\references{Klugman, Panjer & Willmot, Loss Models, Second
  Edition, Wiley, 2004.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Mathieu Pigeon
}
\examples{
dmod <- distrmod(list(dist="gamma", par=c("shape", "rate")), c(2, 3), 5, 10,
TRUE, TRUE)
curve(dmod)
}
\keyword{models}
