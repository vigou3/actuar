
### ===== actuar: an R package for Actuarial Science =====
###
### Ogive and histogram for grouped data
###
### See Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Mathieu Pigeon

## Method to create graphic of empirical distribution function.
plot.ogive <- function(x, y = NULL, xlim = NULL, ylim = NULL, xlab = "boundaries", ylab = "F(x)", col = 1, ...)
{
    xval <- eval(expression(x), env = environment(x))
    plot(xval, x(xval),  main = "Ogive", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = col, type = "o", pch = 20)
}

