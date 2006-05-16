## "cdpareto" <-
## function (x,alpha,lambda)
## {
##     y <- .C("dpareto", as.double(x), as.double(alpha), as.double(lambda), y = double(length(x)), as.integer(length(x)), PACKAGE = "actuar")
## return(y[[4]])
## }

## "cppareto" <-
## function (x,alpha,lambda)
## {
## y<-.C("ppareto",as.double(x),as.double(alpha),as.double(lambda),y = double(length(x)),as.integer(length(x)),PACKAGE = "pareto")
## return(y[[4]])
## }

## "cqpareto" <-
## function (x,alpha,lambda)
## {
## y<-.C("qpareto",as.double(x),as.double(alpha),as.double(lambda),y = double(length(x)),as.integer(length(x)),PACKAGE = "pareto")
## return(y[[4]])
## }

rpareto <- function(n, alpha, lambda)
{
    .C("R_rpareto", as.integer(n), as.double(alpha), as.double(lambda),
       y = double(n))
}
