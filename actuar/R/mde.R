### ===== actuar: an R package for Actuarial Science =====
###
### [...]
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mde <- function(x, fun, start, measure = c("CvM", "chi-square", "LAS"), weights = NULL, ...)
{

    myfn <- function(parm, ...) sum(weights * (G(parm, x, ...) - Gn(x))^2)

    Call <- match.call(expand.dots = FALSE)

    if (missing(start))
        start <- NULL
    dots <- names(list(...))
    dots <- dots[!is.element(dots, c("upper", "lower"))]

    if (missing(fun) || !(is.function(fun) || is.character(fun)))
        stop("'fun' must be supplied as a function or name")

    if (is.null(start) || !is.list(start))
        stop("'start' must be a named list")
    nm <- names(start)
    f <- formals(fun)
    args <- names(f)
    m <- match(nm, args)
    if (any(is.na(m)))
        stop("'start' specifies names which are not arguments to 'fun'")
    formals(fun) <- c(f[c(1, m)], f[-c(1, m)])
    fn <- function(parm, x, ...) fun(x, parm, ...)
    if ((l <- length(nm)) > 1)
        body(fn) <- parse(text = paste("fun(x,", paste("parm[", 1:l, "]", collapse = ", "), ")"))

    if (measure == "CvM")
    {
        x <- sort(x)
        G <- fn
        Gn <- ecdf(x)
        weights <- if (is.null(weights)) 1
        Call$x <- knots(Gn)
        Call$par <- start
    }

    if (measure == "chi-square")
    {
        if (class(x) != "grouped.data")
            stop("'x' must be an object of class 'grouped.data'")
        myfn <- function(parm, ...)
            sum((sum(nj)*diff(fn(parm, cj)) - nj) ^ 2 / weights)
        nj <- x$nj[-1]
        cj <- x$cj
        if(is.null(weights))
            weights <- nj
        if(any(nj == 0))
            stop("all class must contain at least one element")
        Call$par <- start
    }
    if (measure == "LAS")
    {
        if(class(x) != "grouped.data")
            stop("'class' of 'x' must be 'grouped.data'")
        myfn <- function(parm, ...) sum(weights[-c(1, n)] * (diff(c(0, fn(parm, cj[-c(1, n)]))) - diff(c(0, empLEV))) ^ 2)
        cj <- x$cj
        nj <- x$nj[-1]
        n <- length(nj)
        empLEV <- cumsum(nj / n * (cj[-1] + cj[-n]) / 2)[-(n - 1)] + cj[-c(1, n)] * (1 - cumsum(nj / sum(nj)))[-c(n)]
        if(is.null(weights))
            weights <- rep(1, n + 1)
        Call$par <- start
      }
    ## optim() call
    Call[[1]] <- as.name("optim")
    Call$fun <- Call$start <- Call$measure <- Call$weights <- NULL
    Call$fn <- myfn
    Call$hessian <- TRUE
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
    return(res)
    
    if (res$convergence > 0)
        stop("optimization failed")
    sds <- sqrt(diag(solve(res$hessian)))
    structure(list(estimate = res$par, sd = sds, loglik = -res$value))
}

### Données "grouped dental" individualisée pour contrôler avec
### l'exemple 2.21 de Loss Models (1ere édition).
gd <- c(rep(25, 30), rep(50, 31), rep(100, 57), rep(150, 42),
        rep(250, 65), rep(500, 84), rep(1000, 45), rep(1500, 10),
        rep(2500, 11), rep(4000, 3))

id <- c(141, 16, 46, 40, 351, 259, 317, 1511, 107, 567)

distanceLAS <- function (p,dons,poids,dist)
{
    f1<-match.fun(paste("lev",dist$dist,sep=""))
    formals(f1)[dist$par] <- p
    n<-length(dons$nj)
    sum(poids[-c(1,n)]*(diff(c(0,f1(dons$cj[-c(1,n)])))-diff(c(0,emp.lev.moments(dons))))^2)
}

distanceChi2 <- function(p,dons,dist)
  {
    f1 <- match.fun(paste("p", dist$dist, sep=""))
    formals(f1)[dist$par] <- p
    nj <- dons$nj[-1]
    cj <- dons$cj
    n <- length(nj) 
    sum((sum(nj)*diff(f1(cj))-nj)^2/nj)
  }

distanceCvM <- function(p,dons,poids,dist)
  {
    f1 <- match.fun(paste("p",dist$dist,sep=""))
    formals(f1)[dist$par]<-p
    n <- length(dons)
    repart <- (1:n)/n
    sum(poids*(f1(dons)-repart)^2)
  }

