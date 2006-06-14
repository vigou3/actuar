WHgen <- function(x, mean, var, skewness)
{
    call <- match.call()
    b1 <- skewness*skewness/108
    b2 <- skewness/6 - 6/skewness
    b3 <- 2/skewness
    FF <- mean + (b1 * (rnorm(x, 0, 1) - b2)^3 - b3) *  sqrt(var)
    FF
    
}
