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

mde <-function (x, cdf, start, methods, w, ...) 
{
  if(methods == "CvM")
    {
      myfn <- function(parm, ...) sum(w * (fun(parm, k) - e(k)) ^ 2)
      Call <- match.call(expand.dots = FALSE)
      if (missing(start)) 
        start <- NULL
      if (missing(x) || length(x) == 0 || mode(x) != "numeric") 
        stop("'x' must be a non-empty numeric vector")
      if (missing(cdf) || !(is.function(cdf) || is.character(cdf))) 
        stop("'cdf' must be supplied as a function or name")
      if (is.null(start) || !is.list(start)) 
        stop("'start' must be a named list")
      nm <- names(start)
      f <- formals(cdf)
      args <- names(f)
      m <- match(nm, args)
      if (any(is.na(m))) 
        stop("'start' specifies names which are not arguments to 'cdf'")
      formals(cdf) <- c(f[c(1, m)], f[-c(1, m)])
      fun <- function(parm, x, w, ...) cdf(x, parm, w, ...)
      if ((l <- length(nm)) > 1) 
        body(fun) <- parse(text = paste("cdf(x,", paste("parm[", 
                             1:l, "]", collapse = ", "), ")"))
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
    }
  if(methods == "chi-square")
    {
      myfn <- function(parm, ...) sum((sum(nj) * diff(fun(parm, cj)) - nj) ^ 2 / w )
      nj <- x$nj[-1]
      cj <- x$cj
      if(any(nj == 0))
        stop("all class must contain at least one element")
      Call <- match.call(expand.dots = FALSE)
      if (missing(start)) 
        start <- NULL
      if (missing(cdf) || !(is.function(cdf) || is.character(cdf))) 
        stop("'cdf' must be supplied as a function or name")
      if (is.null(start) || !is.list(start)) 
        stop("'start' must be a named list")
      nm <- names(start)
      f <- formals(cdf)
      args <- names(f)
      m <- match(nm, args)
      if (any(is.na(m))) 
        stop("'start' specifies names which are not arguments to 'cdf'")
      formals(cdf) <- c(f[c(1, m)], f[-c(1, m)])
      fun <- function(parm, x, ...) cdf(x, parm, ...)
      if ((l <- length(nm)) > 1) 
        body(fun) <- parse(text = paste("cdf(x,", paste("parm[", 
                             1:l, "]", collapse = ", "), ")"))
      Call[[1]] <- as.name("optim")
      Call$cdf <- Call$start <- NULL
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
    }
  if(methods == "LAS")
    {
      myfn <- function(parm, ...) sum(w[-c(1,n)] * (diff(c(0,fun(parm, cj[-c(1,n)]))) - diff(c(0,emp.lev.moments(x)))) ^ 2)
      cj <- x$cj
      nj <- x$nj
      n <- length(nj)
      Call <- match.call(expand.dots = FALSE)
      if (missing(start)) 
        start <- NULL
      if (missing(cdf) || !(is.function(cdf) || is.character(cdf))) 
        stop("'cdf' must be supplied as a function or name")
      if (is.null(start) || !is.list(start)) 
        stop("'start' must be a named list")
      nm <- names(start)
      f <- formals(cdf)
      args <- names(f)
      m <- match(nm, args)
      if (any(is.na(m))) 
        stop("'start' specifies names which are not arguments to 'cdf'")
      formals(cdf) <- c(f[c(1, m)], f[-c(1, m)])
      fun <- function(parm, x, ...) cdf(x, parm, ...)
      if ((l <- length(nm)) > 1) 
        body(fun) <- parse(text = paste("cdf(x,", paste("parm[", 
                             1:l, "]", collapse = ", "), ")"))
      Call[[1]] <- as.name("optim")
      Call$cdf <- Call$start <- NULL
      Call$par <- start
      Call$x <- x
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
    }
  res
}

test <- function (p,dons,dist) 
{
f1<-match.fun(paste("lev",dist$dist,sep=""))
formals(f1)[dist$par]<-p
n<-length(dons$nj)
sum(diff(c(0,f1(dons$cj[-c(1,n)])))-diff(c(0,emp.lev.moments(dons))))
}


