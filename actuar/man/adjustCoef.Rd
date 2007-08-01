\name{adjustCoef}
\title{Ruin Theory}
\description{
	Compute the adjustment coefficient R, the (strictly) positive
	root of the Lundberg equation :
	h(r) = E[ e^(r X - r c W) ] = 1
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
}
\details{
We want to solve the equation 
\deqn{ \underbrace{E[ e^{r X - r c W} ]}_{h(r)} = 1}
on \eqn{[0, R_{bound}]} where \eqn{X} stands for the claim size distribution, \eqn{W} the inter-occurence times
and \eqn{c} the premium rate, which must satisfy the positive safety loading constraint \eqn{E[X-cW]<0}. The
moment generating function of \eqn{X} on \eqn{[0, R_{bound}]}.

There are 2 ways to use
the function \code{adjustCoef}, either to supply the \eqn{h} function in the argument \code{model} or to use the
pre-coded \eqn{h} function. If you give the \eqn{h} function, then \code{param} must be a list with the
component \code{RBound}, for \eqn{R_{bound}}. Otherwise, you must give the value of the parameters of the
pre-coded \eqn{h} function.

The pre-coded \eqn{h} function are :
    	if \code{model} is the string \code{CramerLundberg}, 
	\deqn{h(r) = \frac{\lambda}{\lambda+rc} M_X(r)}
	where \eqn{M_X} denotes the moment generating function of \eqn{X}. \eqn{M_X(r)} equals to 
	\itemize{
		\item{\deqn{M_X(r) = \frac{\beta}{\beta-r}} when \eqn{X} is exponentially distributed. \code{marginal} must be 
			\code{exp} and \code{param} must be a list with the following components: \code{lambda} and \code{beta}.
			In this case, \eqn{R_{bound}=\beta}. }
		\item{\deqn{M_X(r) = \sum_{i=1}^m \frac{w_i \beta_i}{\beta_i-r}} when \eqn{X} is a mixture of exponential.
			\code{marginal} must be \code{mixexp} and \code{param} must be a list with the following components:
			\code{lambda}, \code{weight} and \code{beta}. Note that \code{weight} and \code{beta} must be a vector
			of the same length. In this case, \eqn{R_{bound}=\underset{i}{\min}\beta_i}.}
		\item{ \deqn{M_X(r) = \left(\frac{ \beta}{\beta-r}\right)^{\alpha}} when \eqn{X} is gamma distributed. 
		       \code{marginal} must be \code{gamma} and \code{param} must be a list with the following components:
		       \code{lambda}, \code{alpha} and \code{beta}. In this case, \eqn{R_{bound}=\beta}.}	
		\item{ \deqn{M_X(r) = \sum_{i=1}^m w_i \left(\frac{ \beta_i}{\beta_i-r}\right)^{\alpha_i}} when \eqn{X} is
		       a mixture of gamma distribution. \code{marginal} must be \code{mixgamma} and \code{param} must
		       be a list with the following components: \code{lambda}, \code{alpha}, \code{beta} and \code{weight}.
		       In this case, \eqn{R_{bound}=\underset{i}{\min}\beta_i}.}
		\item{ \deqn{M_X(r) = \pi (-rI_m -T)^{-1} t_0 } with \eqn{I_m} stands for the identity matrix of dimension
		       \eqn{m}, and \eqn{t_0} denotes the exit rate vector \eqn{-T \mathbf{1}_m} where \eqn{\mathbf{1}_m}
		       	is the \eqn{1} vector of dimension \eqn{m}, when \eqn{X} is phase-type distributed \eqn{PH(\pi,T,m)}.
			\code{marginal} must be \code{phasetype} and \code{param} must be a list with the following 
			components:  \code{lambda}, \code{pi}, \code{T} and \code{m}.  In this case, 
			\eqn{R_{bound}=-\underset{i}{\min} diag(T)}, where \eqn{eqn} denotes the diagonal of a matrix.}	
			}	
	
	if \code{model} is the string \code{SparreAndersen}, 
	\deqn{h(r) = M_W{-rc} M_X(r)}
	where \eqn{M_X} denotes the moment generating function of \eqn{X} and \eqn{M_W} the moment generating
	function \eqn{W}. Theses 2 moment generating function are equals to
	\itemize{
		\item{ \deqn{M_X(r)=\frac{\beta}{\beta-r}~\textrm{and}~M_W(r)=\left(\frac{\lambda}{\lambda-r}\right)^\eta}
		       when \eqn{X} is exponentially distributed, and \eqn{W} is gamma distributed.\code{marginal} must be 
			\code{expgamma} and \code{param} must be a list with the following components: 
			\code{eta}, \code{lambda} and \code{beta}.  In this case, \eqn{R_{bound}=\beta.}
		\item{ \deqn{M_X(r)=\left(\frac{\beta}{\beta-r}\right)^\alpha~\textrm{and}~
		       M_W(r)=\left(\frac{\lambda}{\lambda-r}\right)^\eta} when \eqn{X} and \eqn{W} are gamma distributed. 
		       \code{marginal} must be \code{gammagamma} and \code{param} must be a list with the following 
		       components: \code{eta}, \code{lambda}, \code{alpha} and \code{beta}. 
		       In this case, \eqn{R_{bound}=\beta.}
		\item{ \deqn{M_X(r)=\pi (-rI_m -T)^{-1} t_0 ~\textrm{and}~ M_W(r)= \nu (-rI_n -S)^{-1} s_0} when 
		       \eqn{X} and \eqn{W} are phase-type distributed respectively \eqn{PH(\pi,T,m)} and \eqn{PH(\nu,S,n)}.
		       Note that \eqn{I_m} and \eqn{I_n} stands for the identity matrix of dimension with the respective
		       dimensions \eqn{m} and \eqn{n}; and \eqn{t_0=-T\mathbf{1}_m} and \eqn{s_0=-S\mathbf{1}_n} 
		       are the exit rate vector for \eqn{X} and \eqn{W}. \code{marginal} must be \code{phasetype} and 
		       \code{param} must be a list with the following components:
		        \code{nu}, \code{S}, \code{n}, \code{pi}, \code{T}, \code{m}. In this case, 
			\eqn{R_{bound}=-\underset{i}{\min} diag(T)}, where \eqn{eqn} denotes the diagonal of a matrix.}      		
			}

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
    
    premium <- -sum(weight %*% solve(matrixT,diag(6)))*2*1.1

    adjustCoef("CramerLundberg","phasetype",list(lambda=2,pi=weight, T=matrixT,m=6),premium)

## NB. we could have just the other parametrization with 'marginal' argument equals to 'mixgamma'

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
\keywords{models}
