"cppareto" <-
function (x,alpha,lambda) 
{
y<-.C("ppareto",as.double(x),as.double(alpha),as.double(lambda),y = double(length(x)),as.integer(length(x)),PACKAGE = "pareto")
return(y[[4]])
}

