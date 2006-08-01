### ===== actuar: an R package for Actuarial Science =====
###
### [...]
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mde <- function(x, fun, start, measure = c("CvM", "chi-square", "LAS"), weights = NULL, ...)
{
    ## General form of the function to minimize.
    myfn <- function(parm, x, weights, ...)
    {
        y <- G(parm, x, ...) - Gn(x)
        drop(crossprod(weights * y, y))
    }

    ## Extract call; used to build the call to optim().
    Call <- match.call(expand.dots = TRUE)

    ## Argument checking.
    if (missing(start) || !is.list(start))
        stop("'start' must be a named list")
    if (missing(fun) || !(is.function(fun)))
        stop("'fun' must be supplied as a function")
    grouped <- inherits(x, "grouped.data")
    if (!(is.numeric(x) || grouped))
        stop("'x' must be a numeric vector or an object of class 'grouped.data'")

    ## Make sure that any argument of 'fun' specified in '...' is held
    ## fixed.
    dots <- names(list(...))
    dots <- dots[!is.element(dots, c("upper", "lower"))]
    start <- start[!is.element(names(start), dots)]

    ## Adapt 'fun' to our needs; taken from MASS::fitdistr.
    nm <- names(start)
    f <- formals(fun)
    args <- names(f)
    m <- match(nm, args)
    if (any(is.na(m)))
        stop("'start' specifies names which are not arguments to 'fun'")
    formals(fun) <- c(f[c(1, m)], f[-c(1, m)]) # reorder arguments
    fn <- function(parm, x, ...) fun(x, parm, ...)
    if ((l <- length(nm)) > 1)
        body(fn) <- parse(text = paste("fun(x,", paste("parm[", 1:l, "]", collapse = ", "), ")"))

    measure <- match.arg(measure)

    ## Cramer-von Mises. Use the true and empirical cdf for individual
    ## data, or the true cdf and the ogive for grouped data.
    if (measure == "CvM")
    {
        G <- fn
        Gn <- if (grouped) ogive(x) else ecdf(x)
        if (is.null(weights))
            weights <- 1
        Call$x <- knots(Gn)
        Call$par <- start
    }

    ## Modified Chi-square.
    if (measure == "chi-square")
    {
        if (!grouped)
            stop("'chi-square' measure requires an object of class 'grouped.data'")
        if (any(x$nj[-1] == 0))
            stop("frequency must be larger than 0 in all classes")
        og <- ogive(x)
        x <- knots(og)
        G <- function(...) diff(fn(...))
        Gn <- function(...) diff(og(...))
        if (is.null(weights))
            weights <- 1/og(x)[-1]
        Call$x <- x
        Call$par <- start
    }

    ## Layer average severity.
    if (measure == "LAS")
    {
        if (!grouped)
            stop("'LAS' measure requires an object of class 'grouped.data'")
        e <- elev(x)
        x <- knots(e)
        G <- function(...) diff(fn(...))
        Gn <- function(...) diff(e(...))
        if (is.null(weights))
            weights <- 1
        Call$x <- x
        Call$par <- start
    }

    ## optim() call
    Call[[1]] <- as.name("optim")
    Call$fun <- Call$start <- Call$measure <- NULL
    Call$fn <- myfn
    Call$weights <- weights
    Call$hessian <- FALSE
    if (is.null(Call$method))
    {
        if (any(c("lower", "upper") %in% names(Call)))
            Call$method <- "L-BFGS-B"
        else if (length(start) > 1)
            Call$method <- "BFGS"
        else
            Call$method <- "Nelder-Mead"
    }
    res <- eval(Call)

    ## Return result
    if (res$convergence > 0)
        stop("optimization failed")
    structure(list(estimate = res$par, distance = res$value), class = c("mde","list"))
}

print.mde <- function(x, digits = getOption("digits"), ...)
{
    ans <- format(rbind(x$estimate, "      distance", x$distance), digits = digits)
    ans[1, ] <- sapply(ans[1, ], function(x) paste("", x))
    ans[3, ] <- sapply(ans[3, ], function(x) paste("(", x, ")", sep=""))
    dn <- dimnames(ans)
    dn[[1]] <- rep("", 3)
    dn[[2]] <- paste(substring("      ", 1, (nchar(ans[3,]) - nchar(dn[[2]])) %/% 2), dn[[2]])
    dn[[2]] <- paste(dn[[2]], substring("      ", 1, (nchar(ans[3,]) - nchar(dn[[2]])) %/% 2))
    ans[2,(2:ncol(ans))] <- ""
    ans[3,(2:ncol(ans))] <- ""
    dimnames(ans) <- dn
    print(ans, quote = FALSE)
    x
}


########### Junkyard #############

### Données "grouped dental"  individualisée pour contrôler avec
### l'exemple 2.21 de Loss Models (1ere édition).
gd <- grouped.data(x = c(0, 25, 50, 100, 150, 250, 500, 1000, 1500, 2500, 4000),
                   y = c(30, 31, 57, 42, 65, 84, 45, 10, 11, 3))
#
id <- c(141, 16, 46, 40, 351, 259, 317, 1511, 107, 567)
