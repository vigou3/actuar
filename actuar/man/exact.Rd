\name{exact}
\alias{exact}
\title{Exact Calculation of the Total Amount of Claims Ditribution}
\description{
  Directly compute the exact distribution of the total amount of claims
  by convolution.
}
\usage{
exact(x.scale = 1, fx, pn)
}

\arguments{
  \item{x.scale}{ integer, the value of an amount of 1 in the severity
  of model, expressed through the vector of probabilities \code{fx}}  
  \item{fx}{ a vector of the (discretized) claim amount distribution;
    first element *must* be \eqn{f_X(0)}}
  \item{pn}{ a vector of the number of claims probabilities; first
    element must be Pr[N = 0]}
}
\details{
  The distribution function of the total amount of claims \code{S} has
  formula:

  \deqn{Fs(x) = \sum_{n=0}^{\infty} p_n F_X^{*n}(x)}{
    Fs(x) = sum(n=0 to inf;pn Fx^*n(x))}
  
  where

  \deqn{p_n = Pr(X \le x)}{pn = Pr(X <= x)}

  and

  \deqn{F_X^{*k}(x) = \sum_{y=0}^{x} F_X^{*(k-1)}(x-y)f_X(y)}{
    Fx^*k(x) = sum(y=0 to x;Fx^*(k-1)(x-y) fx(y))}
  
  for \eqn{x = 0, 1,..., k = 2, 3,....}

  Since this convolution product for discrete distributions is defined
  for unitary amounts of claims, the vector of probabilities \code{fx}
  \emph{must} have elements 1, 2, 3, ... respectively equal to \eqn{f_X(0)},
  \eqn{f_X(1)}, \eqn{f_X(2)}, and so on. It is to the user's discretion
  to choose a scale for the amounts that will give a relatively good
  precision without demanding considerable computer time. It is well
  known that the n-fold convolutions are computer intensive, thereupon
  the \code{actuar} package proposes a series of approximation methods
  (see below) to compute CDFs.
}
\value{
 An empirical distribution function of class \code{aggregateDist}. For details on
 such objects, refer to the more general function \code{\link{aggregateDist}}. 
}
\seealso{
  \code{\link{aggregateDist}} \code{\link{panjer}} \code{\link{normal}}
  \code{\link{np2}} \code{\link{simS}}
}
\references{  Klugman, S.A and Panjer,H.H and Willmot, G.E (2004),
  \emph{Loss Models: from data to decision, Second Edition}, Wiley, Sections
  4.6.5 - 4.6.7, 6.6 and Appendix B}

\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Louis-Philippe Pouliot
  }

  \examples{
x.scale <- 100
limsev <- qbeta(1 - 1E-5, 2, 5)                         ### Discretization of
Fx <- pbeta(seq(from=0, to=limsev, by=1/x.scale), 2, 5) ### the severity distribution
fx <- c(0, diff(Fx))


limfreq <- qpois(1-1E-5, 3)                      ### Discretization of
Pn <- ppois(seq(from=0, to=limfreq), 3)          ### the frequency distribution
pn <- c(0, diff(Pn))

(d <- exact(x.scale, fx, pn)) 
}
\keyword{distribution}
\keyword{models}
