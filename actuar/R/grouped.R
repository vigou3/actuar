
### ===== actuar: an R package for Actuarial Science =====
###
### Ogive and histogram for grouped data
###
### See Klugman, Panjer & Willmot, Loss Models, Second
### Edition, Wiley, 2004.
###
### AUTHORS: Vincent Goulet <vincent.gouletaact.ulaval.ca> Mathieu Pigeon
 

grouped <- function(x, y = NULL, digits = 2)
{
  if(class(x) == "data.frame"){
    y <- x[, 2]
    x <- x[, 1]
  }
  nx <- length(x)
  ny <- length(y)
  if(nx - ny > 1 || nx - ny < 0)
    stop("length(x) incorrect")
  if(nx - ny == 1)
      y = c(0, y)
  if(nx == ny && y[1] != 0)
    stop("length(y) incorrect")
   
  res = list(cj = x, nj = y, digits = digits)
  class(res) <- c("grouped.data")
  attr(res, "call") <- sys.call()
  res
  }

ogive <- function(x, y = NULL)
{
    if(class(x) == "grouped.data"){
      y <- x$nj
      x <- x$cj
      }
    if (length(x) < 1) 
        stop("'x' must have 1 or more non-missing values")
    Fnt <- approxfun(x, cumsum(y) / sum(y), yleft = 0, yright = 1, method = "linear", ties = "ordered")
    class(Fnt) <- c("ogive", class(Fnt))
    attr(Fnt, "call") <- sys.call()
    Fnt
}

hist.grouped.data <- function (x, y = NULL, main = "Histogram", xlim = NULL, ylim = NULL, xlab = "boundaries", ylab = "f(x)", plot = TRUE, ...)
{
  if(class(x) == "grouped.data"){
     y <- x$nj
     x <- x$cj
   }
  fnt <- approxfun(x, c(0, y[-1] / (sum(y) * diff(x))), yleft = 0, yright = 0, f = 1, method = "constant")
  r <- structure(list(cj = x, nj = y, density = fnt(x)), class = "histogram")
  if(plot){
    plot(x , fnt(x), main = main, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, type = "S", frame = FALSE)
    segments(x, 0, x, fnt(x))
    segments(0, 0, max(x), 0)
    invisible(r)
  }
  else{
    r
  }
}

plot.ogive <- function (x, y = NULL, xlim = NULL, ylim = NULL, xlab = "boundaries", ylab = "F(x)", col = 1, ...)
{
  xval <- eval(expression(x), env = environment(x))
  plot(xval, x(xval),  main = "Ogive", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = col, type = "o", pch = 20)
}

print.ogive <- function (x, digits = getOption("digits") - 2, ...)
{
  numform <- function(x) paste(formatC(x, dig = digits), collapse = ", ")
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
  digits <- x$digits
  numform <- function(x) formatC(x, digits = digits, format = "e")
  numformy <- function(x) formatC(x)
  y <- x$nj
  x <- x$cj
  x1 <- length(x)
 
    cat("          cj     ", "          nj       ", "\n",paste("[",numform(x[-x1]),", ",numform(x[-1]),")", "       ", numformy(y[-1]), "\n", sep = ""))
}

"[.grouped.data" <- function(x, is = 1, ie = NULL)
{
  if(is.null(ie)){
    ie <- length(x$cj) - 1
  }
  is <- is + 1
  ie <- ie + 1
  nj <- x$nj
  nj <- nj[(is:ie)]
  cj <- x$cj
  is <- is - 1
  cj <- cj[(is:ie)]
  digits <- x$digits
  res = list(cj = cj, nj =  c(0, nj), digits = digits)
  class(res) <- c("grouped.data")
  attr(res, "call") <- sys.call()
  res   
}
  

  
  


