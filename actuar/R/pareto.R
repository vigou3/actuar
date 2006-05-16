## dpareto <- function (x, alpha, lambda)
##     .C("R_dpareto", as.double(x), as.double(alpha), as.double(lambda),
##        y = double(length(x)), as.integer(length(x)))$y

## ppareto <- function(x, alpha, lambda)
##     .C("R_ppareto", as.double(x), as.double(alpha), as.double(lambda),
##        y = double(length(x)), as.integer(length(x)))$y

## qpareto <- function(q, alpha, lambda)
##     .C("R_qpareto", as.double(q), as.double(alpha), as.double(lambda),
##        y = double(length(x)), as.integer(length(x)))$y

rpareto <- function(n, alpha, lambda)
    .C("rpareto_sym", as.integer(n), as.double(alpha), as.double(lambda),
       y = double(n))$y
