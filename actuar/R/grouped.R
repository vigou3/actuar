
### ===== actuar: an R package for Actuarial Science =====
###
### Ogive et histogramme pour données groupées

grouped <- function(x, y = NULL)
{
  if(class(x) == "data.frame"){
    y <- x$nj
    x <- x$cj
  }

  Fnt <- approxfun(x, cumsum(y) / sum(y), yleft = 0, yright = 1, method = "linear")
  fnt <- approxfun(x, c(0, y[-1] / (sum(y) * diff(x))), yleft = 0, yright = 0, f = 1, method = "constant")

  result <- list(cj = x, Fnt = Fnt(x), fnt = fnt(x))
  class(result) <- "groupedData"
  result
  }

ogive <- function(x, y = NULL, xlim = NULL, ylim = NULL, xlab = "group boundaries", ylab = "F(x)", col = 1)
{
  data <- grouped(x, y)
  plot(data, xlim = xlim, ylim = ylim, xlab = "group boundaries", ylab = "F(x)", col = 1)
}

hist.groupedData <- function (x, xlim = NULL, ylim = NULL, xlab = "group boundaries", ylab = "f(x)", ...)
{
  cj <- x$cj
  fnt <- x$fnt

  plot(cj, fnt, main = "Histogram", xlim = xlim, ylim = ylim, xlab = xlab, ylab =ylab, type = "S", frame = FALSE) 
  segments(cj, 0, cj, fnt)             
  segments(0, 0, max(cj), 0)
}


plot.groupedData <- function (x, y = NULL, xlim = NULL, ylim = NULL, xlab = "group boundaries", ylab = "F(x)", col = 1, ...)
{
  cj <- x$cj
  Fnt <- x$Fnt

  plot(cj, Fnt, main = "Ogive", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, type = "o", pch = 20, col = col)
}
