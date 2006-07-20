### ===== actuar: an R package for Actuarial Science =====
###
### Definition of the {m,lev}exp functions. The Exponential
### distribution used in these functions has cumulative distribution
### function
###
###   Pr[X <= x] = 1 - exp(-x / scale), x > 0.
###
### See Appendix A of Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mde <-function (x, cdf, start, w, ...) 
{
    myfn <- function(parm, ...) sum(w * (dens(parm, ...) - e(k)) ^ 2)
    Call <- match.call(expand.dots = FALSE)
    if (missing(start)) 
        start <- NULL
    if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
        stop("'x' must be a non-empty numeric vector")
    if (missing(cdf) || !(is.function(cdf) || is.character(cdf))) 
        stop("'cdf' must be supplied as a function or name")
    n <- length(x)
    if (is.null(start) || !is.list(start)) 
        stop("'start' must be a named list")
    nm <- names(start)
    f <- formals(cdf)
    args <- names(f)
    m <- match(nm, args)
    if (any(is.na(m))) 
        stop("'start' specifies names which are not arguments to 'cdf'")
    formals(cdf) <- c(f[c(1, m)], f[-c(1, m)])
    dens <- function(parm, x, w, ...) cdf(x, parm, w, ...)
    if ((l <- length(nm)) > 1) 
        body(dens) <- parse(text = paste("cdf(x,", paste("parm[", 
            1:l, "]", collapse = ", "), ", w, ...)"))
    e <- ecdf(x)
    k <- knots(e)
    Call[[1]] <- as.name("optim")
    Call$cdf <- Call$start <- NULL
    Call$x <- k
    Call$par <- start
    Call$fn <- myfn
    Call$hessian <- TRUE
    if (is.null(Call$method)) {
        if (any(c("lower", "upper") %in% names(Call))) 
            Call$method <- "L-BFGS-B"
        else if (length(start) > 1) 
            Call$method <- "BFGS"
        else Call$method <- "Nelder-Mead"
    }
    res <- eval(Call)
    if (res$convergence > 0) 
        stop("optimization failed")
    res
}


