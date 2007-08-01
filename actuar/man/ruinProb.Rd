\name{ruinProb}
\title{Ruin Theory}
\description{
	Compute the infinite time ruin probability in the model
	of Cramér-Lundberg or Sparre Andersen, using the results
	of Gerber-Dufresnes (1988) and Asmussen-Rolski (1991).
}
\usage{
ruinProb(model="CramerLundberg",param,premRate)
}
\arguments{
	\item{model}{a string indicating the model used, either Cramér-Lundberg (default value) or
	Sparre Andersen}
	\item{param}{a list with the following components: either \code{lambda}, \code{pi}, \code{T} and
	 \code{m} in the Cramér-Lundberg model or \code{nu}, \code{S}, \code{n}, \code{pi}, \code{T}, 
	 \code{m} in the Sparre Andersen model; where \code{lambda} is the Poisson process parameter
	 (of the claim arrival process), \code{pi,T,m} the parameters of the claim size phase-type distribution and
	 \code{nu,S,n} the parameters of the inter-occurence times phase-type distribution}
	 \item{premRate}{the premium rate (which must respect the positive safety loading constraint)}
}
\details{
TODO
}
\value{
Function \code{ruinProb} computes the ruin probability and returns the ruin probability
as a function of one parameter, the initial capital.
}
\references{
Asmussen, S. and Rolski, T. (1991). \emph{Computational methods in risk theory:
A matrix algorithmic approach}. Insurance: Mathematics and Economics, 10:259-274

Gerber, H. U. and Dufresnes, F. (1991). \emph{Three methods to calculate the probability
of ruin}. Astin Bull., 19(1):71-90	
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  and Christophe Dutang
}
\examples{
## Cramer Lundberg - exponential claim sizes
    resExpPsi <- ruinProb("CramerLundberg", list( lambda=3 ,pi=1, T=2, m=1),  1.1*3/2)

    gridU <- seq(0,5/2,length.out=100) #vector of initial capitals
    resExp <- resExpPsi(gridU)

## Cramer Lundberg - mixture of 3 exponential or hyper-exponential claim sizes
## E(1), E(3) and E(6) with respective weights 1/3, 1/2, 1/6
    Lambda <- 3
    matT <- array(c(-1,0,0,0,-3,0,0,0,-6),c(3,3))
    weight <- c(1/3,1/2,1/6)
    beta <- c(1,3,6)
    premiumRate <- sum(weight/beta)*Lambda*1.1
        
    resMixExpPsi <- ruinProb("CramerLundberg", list( lambda=Lambda ,pi=poids, T=matT, m=3), premiumRate)

    gridU <- seq(0,5*sum(weight/beta),length.out=100)#vector of initial capitals
    resMixExp <- resMixExpPsi(gridU)

## Sparre Andersen - Dickson (1992) numerical applications
## Dickson, D. (1992). On the distribution of surplus prior to ruin. North American Actuarial Journal
    S <-array(c(-1,0,0,1,-1,0,0,1,-1),c(3,3))
    probnu <- c(.4,.2,.4)
    
    matT <- array(c(-1,0,0,0,1,-3,0,0,0,3,-2,0,0,0,2,-4),c(4,4))
    probpi <- c(.2,.3,.4,.1)    

    resErlangPsi <- ruinProb("SparreAndersen", 
    		 list( pi=probpi,T=matT,m=length(probpi),nu=probnu,S=S,n=length(probnu) ), 1)
    
    gridU <- seq(0,5*16/15,length.out=100)
    resErlang <- resErlangPsi(gridU)
}
\keywords{models}
