### ===== actuar: an R package for Actuarial Science =====
###
### Ogive pour données groupées

ogive <- function(x, lines = FALSE, xlim = NULL, ylim = NULL, xlab = "group boundaries", ylab = "F(x)", col = 1)
{
  class(x) <- "ogive"
  
  if (lines == FALSE){
  plot(x, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, col = col)
  }
  else{
  lines(x, xlim = xlim, ylim = ylim, col = col)
  }
}


plot.ogive <- function (x, y = NULL, xlim = NULL, ylim = NULL, xlab = "group boundaries", ylab = "F(x)", col = 1, ...)
{
  cj <- x$cj
  nj <- x$nj

  Fnt <- approxfun(cj, cumsum(nj) / sum(nj),
                 yleft = 0, yright = 1, method = "linear")

  plot(cj, Fnt(cj), main = "Ogive", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, type = "o", pch = 20, col = col)
  
  p<-sum(is.na(cj*0))
  if(p >= 1){
	warning("Inf in group boundaries")
      }
}

lines.ogive <- function (x, xlim = NULL, ylim = NULL, col = 1, ...)
{
  cj <- x$cj
  nj <- x$nj
  
  Fnt <- approxfun(cj, cumsum(nj) / sum(nj),
                 yleft = 0, yright = 1, method = "linear")

  lines(cj, Fnt(cj), xlim = xlim, ylim = ylim, col = col, type = "o")
}

  
