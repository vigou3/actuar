### ===== actuar: an R package for Actuarial Science =====
###
### Phase-type Theory
###
### Compute the exponential matrix of exp( T*x ) %*% v 
###
### function used by ruinProb
###
### AUTHORS: Christophe Dutang,
### Vincent Goulet <vincent.goulet@act.ulaval.ca>,

calcMatrixExpGen <- function(x, u, T, v)
{
    .Call("calcMatExpGen", x, u, T, v)
}
 
