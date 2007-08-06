\name{adjustCoef}
\alias{adjustCoef}
\title{Ruin Theory}
\description{
	Compute the adjustment coefficient R, the (strictly) positive
	root of the Lundberg equation :
	\deqn{ h(r) = E[ e^(r X - r c W) ] = 1}
}
\usage{
adjustCoef(model="CramerLundberg",marginal="exp",param,premRate,toplot=FALSE)
}
\arguments{
	\item{model}{ the function 'h' or a string, either  CramerLundberg (default value) or
	SparreAndersen, which use pre-coded 'h' function (see Details)}
	\item{marginal}{a string used, if \code{model} is either CramerLundberg or SparreAndersen, which
	indicates the marginal distribution of claim sizes and arrival times. \code{marginal}
	could be \code{exp}, \code{gamma}, \code{mixexp}, \code{mixgamma}, \code{phasetype},
	\code{expgamma} or \code{gammagamma} (see Details). }
	\item{param}{if \code{model} is a function, \code{param} is a list with the component \code{RBound}, 
	such that \code{h} is defined on [0,\code{RBound}].
		
		if \code{model} is the string \code{CramerLundberg}, \code{param} is a list with the following components: 
	either \code{lambda} and \code{beta} when \code{marginal} equals to \code{exp} ; either \code{lambda},
	\code{weight} and \code{beta} when \code{marginal} equals to \code{mixexp} ; either \code{lambda}, \code{alpha}
	and \code{beta} when \code{marginal} equals to \code{gamma} ; either \code{lambda}, \code{alpha},
	\code{beta} and \code{weight} when \code{marginal} equals to \code{mixgamma} or \code{lambda}, \code{pi}, 
	\code{T} and \code{m} when \code{marginal} equals to \code{phasetype}.
	
		if \code{model} is the string \code{SparreAndersen}, \code{param} is a list with the following components:
	either \code{eta}, \code{lambda} and \code{beta} when \code{marginal} equals to \code{expgamma} ; either
	\code{eta}, \code{lambda}, \code{alpha} and \code{beta} when \code{marginal} equals to \code{gammagamma}
	or \code{nu}, \code{S}, \code{n}, \code{pi}, \code{T}, \code{m} when \code{marginal} equals to \code{phasetype} }
	 \item{premRate}{the premium rate (which must respect the positive safety loading constraint)}
	 \item{toplot}{TO DELETE}
}
\details{
We want to solve the equation 
\deqn{ \underbrace{E[ e^{r X - r c W} ]}_{h(r)} = 1}
on \eqn{[0, R_{bound}]} where \eqn{X} stands for the claim size distribution, \eqn{W} the inter-occurence times
and \eqn{c} the premium rate, which must satisfy the positive safety loading constraint \eqn{E[X-cW]<0}. The
moment generating function of \eqn{X} on \eqn{[0, R_{bound}]}.


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
\dontrun{
## Cramer Lundberg - exponential claim sizes
    adjustCoef(param=list(beta=1.5,lambda=1),premRate=1.1*2/3)

## Cramer Lundberg - mixture of 4 exponential or hyper-exponential claim sizes
## E(1), E(2), E(3) and E(4) with respective weights 0.1, 0.2, 0.3 and 0.4.
    premium <- sum(c(.1,.2,.3,.4)/1:4)*1*1.1
    adjustCoef("CramerLundberg","mixexp",list(lambda=1,weight=c(.1,.2,.3,.4),beta=1:4),premium)

## Cramer Lundberg - phase-type distributed claim sizes
## claim sizes is a mixture of 2 erlang distribution Er(3,1) and Er(3,2)
    matrixT <- array(0,c(6,6))
    weight <- c(1/3,0,0,2/3,0,0)
    diag(matrixT) <- -c(rep(1,3),rep(2,3))
    matrixT[1,2] <- 1
    matrixT[2,3] <- 1
    matrixT[4,5] <- 2
    matrixT[5,6] <- 2
    
    premium <- -sum( weight %*% solve(matrixT, diag(6) ) )*2*1.1

    print(premium)

    adjustCoef("CramerLundberg", "phasetype", list(lambda=2, pi=weight, T=matrixT, m=6), premium)

## NB. we could use the other parametrization with 'marginal' argument equals to 'mixgamma'

## Sparre Andersen - both claim sizes and arrival times are gamma distributed
## respectively G(2,2) and G(3,1)
    adjustCoef("SparreAndersen", "gammagamma", list(lambda=1,eta=3,beta=2,alpha=2), 2/2*1/3*1.1)

## Sparre Andersen - claim size is a mixture of exponential and inter-occurence time is erlang
## respectively a mixture of exponential E(1), E(2) and E(3) with weight 1/3 and an erlang distribution
## with Er(6,1).
    matS <- array(0,c(6,6))
    probnu <- c(1,0,0,0,0,0)
    diag(matS) <- -c(rep(1,6))
    matS[1,2] <- 1
    matS[2,3] <- 1
    matS[3,4] <- 1
    matS[4,5] <- 1
    matS[5,6] <- 1

    matT <- array(0,c(3,3))
    probpi <- c(1/3,1/3,1/3)
    diag(matT) <- -(1:3)
        

    expX <- -sum(probpi%*%solve(matT,diag(3)))
    expW <- -sum(probnu%*%solve(matS,diag(6)))

    adjustCoef("SparreAndersen", "phasetype", 
    	list(nu=probnu,S=matS, n=6, pi=probpi, T=matT, m=3), 1.1*expX/expW)


## the "conditional" model based Boudreault et al (2006)
## as usual, if X are claim sizes, and W are inter-occurence times,
## we suppose in this model, that X and W are dependent:
## M_X/W=t (x) = exp(-alpha t) M_Y1(x) + (1-exp(-alpha t)) M_Y2(x)
## where Y1, Y2 are independent
   
    lambda <- 1.5
    beta1 <- 1
    beta2 <- 2
    alpha <- 1/2
    expX <- lambda/(lambda+alpha)/beta1 + (1-lambda/(lambda+alpha))/beta2
    prem <- expX*lambda*1.1    

    theHfunction <- function(r)
    		 beta1/(beta1-r)*lambda/(lambda+r*prem) + beta2/(beta2-r)*(lambda/(lambda+r*prem) -lambda/(lambda+r*prem+alpha))
		 
    adjustCoef(theHfunction, "", list(RBound=min(beta1,beta2)), prem)

## see Boudreault et al (2006) page 269, complete reference: Boudreault, M., Cossette, H., Landriault, D. and Marceau, E. (2006).
## On a risk model with dependence between interclaim arrivals and claim sizes. Scandinavian Actuarial Journal, 2006(5):265-285
}
}
\keyword{models}
