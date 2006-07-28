
### ===== actuar: an R package for Actuarial Science =====
###
### Ogive and histogram for grouped data
###
### See Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Mathieu Pigeon


grouped.data <- function(x, y = NULL, digits = 2)
{
    ## 'data.frame' must contain boundaries in first column and number of
    ## data by class in second column.
    if (class(x) == "data.frame")
    {
        x <- x[, 1]
        y <- x[, 2]
    }
    nx <- length(x)
    ny <- length(y)

    ## First data in 'y' won't be used.
    if (nx - ny > 1 || nx - ny < 0)
        stop("length(x) incorrect")
    if (nx - ny == 1)
        y = c(0, y)

    ## Create an object of class 'grouped.data'.
    res = list(cj = x, nj = y, digits = digits, j = 3)
    class(res) <- "grouped.data"
    attr(res, "call") <- sys.call()
    res
}

## Calculate an empirical distribution function.
ogive <- function(x, y = NULL)
{
  ## Use object created by 'grouped' function.
  if (class(x) == "grouped.data"){
    y <- x$nj
    x <- x$cj
  }
  ## An error message is issued if 'x' is empty.
  if (length(x) < 1)
    stop("'x' must have 1 or more non-missing values")

  ##Create an object of class 'ogive'.
  Fnt <- approxfun(x, cumsum(y) / sum(y), yleft = 0, yright = 1, method = "linear", ties = "ordered")
  class(Fnt) <- c("ogive", class(Fnt))
  attr(Fnt, "call") <- sys.call()
  Fnt
}

### Method of knots() for objects of class 'ogive'. Identical to
### stats::knots.stepfun().
knots.ogive <- stats:::knots.stepfun

## Calculate an empirical density function and create the histogram.
hist.grouped.data <- function (x, y = NULL, main = "Histogram", xlim = NULL, ylim = NULL, xlab = "boundaries", ylab = "f(x)", plot = TRUE, ...)
{
  ## Use object created by 'grouped' function.
  if (class(x) == "grouped.data"){
    y <- x$nj
    x <- x$cj
  }

  ## Create an object of class 'histogram'.
  fnt <- approxfun(x, c(0, y[-1] / (sum(y) * diff(x))), yleft = 0, yright = 0, f = 1, method = "constant")
  r <- structure(list(cj = x, nj = y[-1], density = fnt(x)), class = "histogram")

  ## If 'plot' is true, histogram is created, else, boudaries, number of data by class and density are returned.
  if (plot){
    plot(x , fnt(x), main = main, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, type = "S", frame = FALSE)
    segments(x, 0, x, fnt(x))
    segments(0, 0, max(x), 0)
    invisible(r)
  }
  else{
    r
  }
}

## Method to create graphic of empirical distribution function.
plot.ogive <- function (x, y = NULL, xlim = NULL, ylim = NULL, xlab = "boundaries", ylab = "F(x)", col = 1, ...)
{
  xval <- eval(expression(x), env = environment(x))
  plot(xval, x(xval),  main = "Ogive", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = col, type = "o", pch = 20)
}

print.ogive <- function (x, digits = getOption("digits") - 2, ...)
{
  ## To formate numbers.
  numform <- function(x) paste(formatC(x, dig = digits), collapse = ", ")

  ## Create a structure for the presentation of empirical distribution function.
  cat("Empirical CDF for grouped data \nCall: ")
  print(attr(x, "call"), ...)
  nc <- length(xxc <- eval(expression(x), env = environment(x)))
  nn <- length(xxn <- eval(expression(y), env = environment(x)))
  i1 <- 1:min(3, nc)
  i2 <- if(nc >= 4) max(4, nc - 1):nc else integer(0)
  i3 <- 1:min(3, nn)
  i4 <- if(nn >= 4) max(4, nn - 1):nn else integer(0)
  cat(" cj = ", numform(xxc[i1]), if(nc > 3) ", ", if(nc > 5) " ..., ", numform(xxc[i2]), "\n", sep = "")
  cat(" nj = ", numform(xxn[i1]), if(nn > 3) ", ", if(nn > 5) " ..., ", numform(xxn[i2]), "\n", sep = "")
  invisible(x)
}

print.grouped.data <- function(x, ...)
{
  ## To formate numbers.
  numform <- function(x) formatC(x, digits = digits, format = "e")
  numformy <- function(x) formatC(x)

  ## Use object created by 'grouped' function
  j <- x$j
  digits <- x$digits

  ## Choose which column(s) is(are) presented.
  if (j == 1){
    x <- x$cj
    x1 <- length(x)
    cat("          cj     ", "\n", paste("[",numform(x[-x1]),", ",numform(x[-1]),")", "\n", sep = ""))
  }
  if (j == 2){
    y <- unclass(x$nj)

  }
  if (j != 1 && j != 2){
    y <- x$nj
    x <- x$cj
    x1 <- length(x)
    cat("          cj     ", "          nj       ", "\n",paste("[",numform(x[-x1]),", ",numform(x[-1]),")", "       ", numformy(y[-1]), "\n", sep = ""))
  }
}

"[.grouped.data" <- function(x, i, j)
{
  ## If 'i' is missing, all rows are presented.
  if (missing(i)){
    i <- (1:(length(x$cj)-1))
  }
  ## If 'j' is missing, all columns are presented.
  if (missing(j)){
    j <- 3
  }

  ## Create an object of class 'grouped.data'.
  if (j == 1){
    cj <- x$cj
    cj <- cj[c((i), max(i) + 1)]
    res = list(cj = cj, digits = x$digits, j = j)
    class(res) <- c("grouped.data")
    attr(res, "call") <- sys.call()
  }
  if (j == 2){
    nj <- x$nj[-1]
    nj <- nj[(i)]
    res = nj
    class(res) <- c("numeric")
  }
  if (j != 1 && j != 2){
  nj <- x$nj[-1]
  nj <- c(0, nj[(i)])
  cj <- x$cj
  cj <- cj[c((i), max(i) + 1)]
  res = list(cj = cj, nj =  nj, digits = x$digits, j = j)
  class(res) <- c("grouped.data")
  attr(res, "call") <- sys.call()
  }
  res
}
