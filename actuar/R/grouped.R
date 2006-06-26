
### ===== actuar: an R package for Actuarial Science =====
###
### Ogive et histogramme pour données groupées

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
  class(res) <- c("data.grouped")
  attr(res, "call") <- sys.call()
  res
  }

ogive <- function(x, y = NULL)
{
    if(class(x) == "data.grouped"){
      y <- x$nj
      x <- x$cj
      }
    if (length(x) < 1) 
        stop("'x' must have 1 or more non-missing values")
    Fnt <- approxfun(x, cumsum(y) / sum(y), yleft = 0, yright = 1, method = "linear", ties = "ordered")
    class(Fnt) <- c("groupedData", class(Fnt))
    attr(Fnt, "call") <- sys.call()
    Fnt
}

hist.data.grouped <- function (x, y = NULL, freq = NULL, main = "Histogram", xlim = NULL, ylim = NULL, xlab = "grouped boundaries", ylab = "f(x)", plot = TRUE, ...)
{
  if(class(x) == "data.grouped"){
     y <- x$nj
     x <- x$cj
   }
  fnt <- approxfun(x, y / (sum(y) * diff(c(0, x))), yleft = 0, yright = 0, f = 1, method = "constant")
  r <- structure(list(cj = x, nj = y, density = fnt(x)), class = "histogram")
  if(plot){
    plot(c(0, x) , fnt(c(0, x)), main = main, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, type = "S", frame = FALSE)
    segments(x, 0, x, fnt(x))
    segments(0, 0, max(x), 0)
    invisible(r)
  }
  else{
    r
  }
}

plot.groupedData <- function (x, y = NULL, xlim = NULL, ylim = NULL, xlab = "group boundaries", ylab = "F(x)", col = 1, ...)
{
  xval <- eval(expression(x), env = environment(x))
  plot(xval, x(xval),  main = "Ogive", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = col, type = "o", pch = 20)
}

print.groupedData <- function (x, digits= getOption("digits") - 2, ...)
{
  numform <- function(x) paste(formatC(x, dig=digits), collapse=", ")
  cat("Empirical CDF for grouped data \nCall: ")
  print(attr(x, "call"), ...)
  nc <- length(xxc <- eval(expression(x),env = environment(x)))
  nn <- length(xxn <- eval(expression(y),env = environment(x)))
  i1 <- 1:min(3, nc)
  i2 <- if(nc >= 4) max(4, nc - 1):nc else integer(0)
  i3 <- 1:min(3, nn)
  i4 <- if(nn >= 4) max(4, nn - 1):nn else integer(0)
  cat(" cj = ", numform(xxc[i1]), if(nc > 3) ", ", if(nc > 5) " ..., ", numform(xxc[i2]), "\n", sep = "")
  cat(" nj = ", numform(xxn[i1]), if(nn > 3) ", ", if(nn > 5) " ..., ", numform(xxn[i2]), "\n", sep = "")
  invisible(x)
}

print.data.grouped <- function(x, ...)
{
  digits <- x$digits
  numform <- function(x) formatC(x, digits = digits, format = "e")
  numformy <- function(x) formatC(x)
  y <- x$nj
  x <- x$cj
  x1 <- length(x)
  
  if(length(x$cj) == length(x$nj))
    cat("          cj     ", "          nj       ", "\n",paste("[",numform(c(0, x[-x1])),", ",numform(x),")", "       ", numformy(y), "\n", sep = ""))
  else
    cat("          cj     ", "          nj       ", "\n",paste("[",numform(x[-x1]),", ",numform(x[-1]),")", "       ", numform(y), "\n", sep = ""))
}

  
  


