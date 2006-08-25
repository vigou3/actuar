\name{predict.ghcm}
\alias{predict.ghcm}
\title{Predict method for Generalized Hierarchical Credibility Models}
\description{Predicted values of claim ratios (credibility premiums)
  based on credibility model object}
\usage{
\method{predict}{ghcm}(object, \dots)
}
\arguments{
  \item{object}{Object of class inheriting from \code{"ghcm"}}
  \item{\dots}{potential further arguments required by generic.}  
}
\details{
  \code{ghcm} objects contain all the estimated values needed to
  calculate the credibility premiums. To give an example of how the
  premiums are computed, let us consider a portfolio divided in sectors,
  which are further divided in classes, which themselves contain
  homogeneous insureds. Let subsctripts \eqn{p, k, i} denote
  respectively the different sectors, classes, and insureds. Then, the
  credibility premium, \eqn{P_i^{(kp)}}{P[i]}, is determined as follows:

  \deqn{P_i^{(kp)} = z_i^{(kp)} B_i^{(kp)} + (1 - z_i^{(kp)})
    \hat{M}_{kp},}{P[i] = z[i] B[i] + (1 - z[i]) M[kp],}
  \deqn{\hat{M}_{kp} = z_k^{(p)} B_k^{(p)} + (1 - z_k^{(p)})
    \hat{M}_{k},}{M[kp] = z[k] B[k] + (1 - z[k]) M[k],}
  \deqn{\hat{M}_{k} = z_p B_p + (1 - z_p) \hat{M},}{M[k] = z[p] B[p] +
    (1 - z[p]) M,}
  where \eqn{B_i^{(kp)}, B_k^{(p)}, B_p}{B[i], B[k] , B[p]} are the
  contract, class, and sector linearly sufficient statictics,
  respectively. \eqn{\hat{M}_{kp}, \hat{M}_{k}, \hat{M}}{M[kp], M[k], M}
  are respectively the class, sector, and portfolio risk premiums.
}
\value{
  Returns a vector of of credibility premiums for each contract.
}
\references{
  Goulet, V. (1998), \emph{Principles and Application of Credibility
    Theory}, Journal of Actuarial Practice, Volume 6, ISSN 1064-6647.
  
  Goovaerts, M. J. and Hoogstad, W. J. (1987), \emph{Credibility Theory},
  Surveys of Actuarial Studies, NO.4, Nationale-Nederlanden N.V.
}
\author{
  Vincent Goulet \email{vincent.goulet@act.ulaval.ca},
  and Louis-Philippe Pouliot
}
\examples{
data(hachemeister)
contract <- c("C1", "C2" ,"C3", "C4", "C5")
sector <- LETTERS[c(1,2,1,2,2)]
x <- data.frame(sector, contract, hachemeister$claims,
hachemeister$weights)

(z <- ghcm(~sector + sector:contract, x, X1:X12, X1.1:X12.1))

predict(z)
}
\keyword{models}
  
  

