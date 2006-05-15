crpareto <- function (n, alpha, lambda)
{
    .C("rpareto", as.integer(n), as.double(alpha), as.double(lambda),
            y = double(n), PACKAGE = "actuar")
}
