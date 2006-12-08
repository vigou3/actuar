### ===== actuar: an R package for Actuarial Science =====
###
### Use one of five methods to compute the aggregate claim amount
### distribution of a portfolio over a period given a frequency and a
### severity model or the true moments of the distribution.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Louis-Philippe Pouliot

discretize <- function (expr, from, to, step = 1, method = c("upper", "lower"),
                        by = step, xlim = NULL)
{
    ## If 'expr' is only the name of a function (say f), build a call
    ## 'f(x)'. Otherwise, check that the expression is a function call
    ## containing an 'x'. Taken from 'curve'.
    sexpr <- substitute(expr)
    if (is.name(sexpr))
    {
        fcall <- paste(sexpr, "(x)")
        expr <- parse(text = fcall)
    }
    else
    {
        if (!(is.call(sexpr) && match("x", all.vars(sexpr), nomatch = 0)))
            stop("'expr' must be a function or an expression containing 'x'")
        expr <- sexpr
    }

    ## If 'from' and/or 'to' are not specified, take their values in 'xlim'.
    if (missing(from))
        from <- xlim[1]
    if (missing(to))
        to <- xlim[2]

    ## Evaluate specified cdf at every 'step' between 'from' and 'to'.
    x <- seq.int(from, to, by)
    Fx <- eval(expr, envir = list(x = x), enclos = parent.frame())

    ## The "upper" discretization is the forward difference of Fx:
    ##
    ##   Fx(step) - Fx(from), ..., Fx(to) - Fx(to - step),
    ##
    ## whereas the "lower" discretization is the backward difference
    ##
    ##   Fx(from), Fx(step) - Fx(from), ..., Fx(to) - Fx(to - step).
    ##
    ## Hence, the latter case simply has one more element.
    method <- match.arg(method)
    c(if(method == "lower") Fx[1], diff(Fx))
}

discretise <- discretize
