\name{panjer}
\alias{panjer}
\title{Panjer's recursion formula}
\description{
  Panjer recursion formula to compute the approximate distribution of the total amount of claims
  of a portfolio. 
}
\usage{
panjer(fx, freq.dist = c("pois", "nbinom", "binom",
"geom", "logarithmic"),par, p0, TOL = 1e-08, echo= FALSE)

}

\arguments{
  \item{fx}{ a vector of the (discretized) claim amount distribution;
    first element *must* be \eqn{f_X(0)}}
  \item{freq.dist}{name of the counting distribution}
  \item{par}{ named list of the parameter(s) of the counting
    distribution as they are defined in "rdist".}
  \item{p0}{arbitrary amount of probability at zero given to the
    frequency distribution. It creates zero-modified or zero-truncated distributions} 
  \item{TOL}{ stop recursion when cumulative probability function is
    less than TOL away from 1}
  \item{echo}{ print the cumulative distribution of the total amount of
    claims as it is computed.}
}
\details{
  The formula of the (a, b, 1) class is:
  \deqn{f_S(x)=\frac{[p_1-(a+b)p_0]f_X(x)+\sum_{y=1}^{\min(x, m)}(a+\frac{by}{x})f_X(y)f_S(x-y)}{1-af_X(0)}%
  }{([p1-(a+b)*p0]*fx[x+1] + sum(y=1 to min(x,m);(a+b*y/x)fx[y+1]*fs(x-y)))/(1-a*fx[1])}
  For the (a, b, 0) class, the result reduces to
  \deqn{f_S(x)=\frac{\sum_{y=1}^{\min(x,m)}(a+\frac{by}{x})f_X(y)f_S(x-y)}{1-af_X(0)}%
  }{sum(y=1 to min(x,m);(a+b*y/x)fx[y+1]*fs(x-y))/(1-a*fx[1])}
  
  The counting variable is a member of the (a, b, 0) family of
  discrete distributions if p0 is not specified and a member of the
  (a,b,1) family if p0 is specified. The logarithmic distribution is
  a limiting case of the negative binomial distribution where the size parameter
  is equal to 0. 
}
\value{
 An object of class \code{AggregateDist}. For details on the elements of
 such objects, refer to the more general function \code{\link{ADgen}}. 
}
\references{  Klugman, S.A and Panjer,H.H and Willmot, G.E (2004),
  \emph{Loss Models: from data to decision, Second Edition}, Wiley, Sections
  4.6.5 - 4.6.7, 6.6 and Appendix B}
\seealso{
  \code{\link{ADgen}} \code{\link{exact}} \code{\link{normal}}
  \code{\link{np2}} \code{\link{simS}}
}

\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca} and
  Sébastien Auclair
  }

\examples{
### (a,b,0) class with the binomial distribution
(a <- panjer(fx=rep(0.5,2), freq.dist="binom",par=list(size=3, prob=0.5),echo=TRUE))
 

### Example 6.18 of "Loss Models" (Second Edition).
### (a,b,1) class with the Extended Truncated Negative Binomial distribution.  
(a <-
panjer(fx=c(0.3, 0.5, 0.2), "nbinom", par=list(size=0.2, prob=0.25), p0=0, TOL=1E-5, echo=TRUE))
plot(a)
}
\keyword{distribution}
\keyword{models}
