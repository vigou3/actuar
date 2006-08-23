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
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  To see how premiums are calculated, let us take a portfolio with
  two hierarchical levels: sector subsctripted j, and contracts
  subsctripted i.
  Let \eqn{B} be 
  
  The risk premium of sector \eqn{j} is 

  \deqn{\hat{M}_j = z_p^{(2)} + (1 - z_i) X_{zw},}{z[i] X[iw] + (1 - z[i]) X[zw],}


  ah putain jsais pas koi mettre
  ..............................
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
  
  

