\name{adjustCoef}
\alias{adjustCoef}
\title{Ruin Theory}
\description{
	Compute the adjustment coefficient R, the (strictly) positive
	root of the Lundberg equation :
	\deqn{ h(r) = E[ e^(r X - r c W) ] = 1}
}
\usage{
adjustCoef (mgfClaim, mgfWaitTime, premRate, upper, h, retRate=1, retLimit=Inf, toplot=FALSE,...)
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

res <- adjustCoef(function(x,y)mgfexp(x*y) , function(x)mgfexp(x), function(x)-.1+x*1.3, 1, retRate=0:5/5)
res(1)

plot(a<-seq(0,1,length.out=100),res(a),type="l")
res <- adjustCoef(function(x,y)mgfexp(x*y) , function(x)mgfexp(x), function(x)-.1+x*1.3, 1, retRate=0:100/100)
lines(a,res(a),type="l",col="blue")
 
premium <-function(x)-.1+x*1.3
theHfunction <- function(x,a) mgfexp(a*x)*mgfexp(-x*premium(a))

res <- adjustCoef(h=theHfunction, upper=1, premRate=premium, retRate=0:10/10)

lines(a,res(a),type="l",col="green")

premium <- function(L) -.1+levgamma(L,2,2)*1.3
mgfClaimSize <- function(x,L) mgfgamma(x,2,2)*pgamma(L,2,2-x)+exp(x*L)*(1-pgamma(L,2,2))

theHfunction <- function(x,L) mgfClaimSize(x,L)*mgfexp(-x*premium(L))

res1 <- adjustCoef(h=theHfunction, upper=1, premRate=premium, retLimit=0:10)


res2 <- adjustCoef(h=theHfunction, upper=1, premRate=premium, retLimit=0:100/10)
plot(L<-seq(0,10,length.out=100),res1(L),type="l",ylim=c(0,max(res1(L),res2(L))))
lines(L,res2(L),type="l",col="blue")

}
\keyword{models}
