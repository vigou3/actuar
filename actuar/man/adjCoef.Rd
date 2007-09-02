\name{adjCoef}
\alias{adjCoef}
\title{Ruin Theory}
\description{
	Compute the adjustment coefficient R, the (strictly) positive
	root of the Lundberg equation :
	\deqn{ h(r) = E[ e^(r X - r c W) ] = 1}
}
\usage{
adjCoef (mgfClaim, mgfWaitTime, premRate, upper, h, retRate=1, retLimit=Inf, toplot=FALSE,...)
}
\arguments{
	TODO
}
\details{
We want to solve the equation 
\deqn{ \underbrace{E[ e^{r X - r c W} ]}_{h(r)} = 1}
on \eqn{[0, upper]} where \eqn{X} stands for the claim size distribution, \eqn{W} the inter-occurence times
and \eqn{c} the premium rate, which must satisfy the positive safety loading constraint \eqn{E[X-cW]<0}. The
moment generating function of \eqn{X} on \eqn{[0, upper]}.


}
\value{
Function \code{adjustCoef} returns the adjustment coefficient. If the return value is \eqn{0}, then the premium
rate \eqn{c} supplied does not satisfy the positive safety loading constraint. 
}
\references{
	Actuarial Mathematics????	
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  and Christophe Dutang
}
\examples{

## no reinsurance, exponentially claim sizes and wait times
res <- adjCoef(function(x) mgfexp(x), function(x) mgfexp(x), 1.2, 1)

## proportional reinsurance,  exponentially claim sizes and wait times
res <- adjCoef(function(x,y)mgfexp(x*y) , function(x)mgfexp(x), function(x)-.1+x*1.3, 1, r="proportional", from=0, to=1, n=6)

##
## R(1) i.e. no reinsurance
res(1)

## plot result
plot(res)

## increase the point number of the retention rate grid
res <- adjCoef(function(x,y)mgfexp(x*y) , function(x)mgfexp(x), function(x)-.1+x*1.3, 1, r="proportional", from=0, to=1, n=101)
lines(res, type="l", col="blue")

## use the other parametrization with argument 'h' 
premium <-function(x)-.1+x*1.3
theHfunction <- function(x,a) mgfexp(a*x)*mgfexp(-x*premium(a))
res <- adjCoef(h=theHfunction, upper=1, p=premium, r="proportional", from=0, to =1, n=11)

## plot the result
## see the impact of point number 'n'
lines(res, type="l", col="green")

##
## excess of loss reinsurance

premium <- function(L) -.1+levgamma(L,2,2)*1.3
mgfClaimSize <- function(x,L) mgfgamma(x,2,2)*pgamma(L,2,2-x)+exp(x*L)*(1-pgamma(L,2,2))
theHfunction <- function(x,L) mgfClaimSize(x,L)*mgfexp(-x*premium(L))

## different grid of retention limits
res1 <- adjCoef(h=theHfunction, upper=1, p=premium, r="excess-of-loss", from=0, to=10, n=11)
res2 <- adjCoef(h=theHfunction, upper=1, p=premium, r="excess-of-loss", from=0, to=10, n=101)

L<-seq(0,10,length.out=100)

## plot the result
plot(res1, type="l", ylim=c(0,max(res1(L),res2(L))))
lines(res2, type="l", col="blue")

}
\keyword{models}
