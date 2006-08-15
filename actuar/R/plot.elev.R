### ===== actuar: an R package for Actuarial Science =====
###
### Create a graphic for sample empirical limited value functions for individual and
### grouped data.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca> and
###          Mathieu Pigeon

plot.elev <- function(x, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL, col = 1, ...)
{
    if (attr(x, "grouped"))  
    {
        plot(knots(x), x(knots(x)), main = "Empirical Limited Function", xlim = xlim, ylim = ylim, xlab = "Limits", ylab = "Empirical limited values", col = col, type = "o", pch = 20)
    }
    else
    {
        xval <- eval(expression(x), env = environment(x))
        plot(xval, x(xval), main = "Empirical Limited Function", xlim = xlim, ylim = ylim, xlab = "Limits", ylab = "Empirical limited values", col = col, type = "o", pch = 20)
    }
}

