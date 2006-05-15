"crpareto" <-
function (x,alpha,lambda) 
{
y<-.C("rpareto",as.double(alpha),as.double(lambda),y = double(x),as.integer(x),package="pareto")
return(y[[3]])
}

