"cdpareto" <-
function (x,alpha,lambda)
{
    y <- .C("dpareto", as.double(x), as.double(alpha), as.double(lambda), y = double(length(x)), as.integer(length(x)), PACKAGE = "actuar")
return(y[[4]])
}
