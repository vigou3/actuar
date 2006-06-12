### ===== actuar: an R package for Actuarial Science =====
###
### Ogive pour données groupées

ogive <- function(data, lines = FALSE, xlim = NULL, ylim = NULL, xlab = "group boundaries", ylab = "F(x)", col = 1)
{
  class(data) <- "ogive"
  
  if (lines == FALSE){
  plot(data, xlim, ylim, xlab, ylab, col)
  }
  else{
  lines(data, xlim = xlim, ylim = ylim, col = col)
  }
}


plot.ogive <- function (data, xlim, ylim, xlab, ylab, col)
{
  cj <- data$cj
  nj <- data$nj

  Fnt <- approxfun(cj, cumsum(nj) / sum(nj),
                 yleft = 0, yright = 1, method = "linear")

  plot(cj, Fnt(cj), main = "Ogive", xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab, type = "o", pch = 20, col = col)
  
  p<-sum(is.na(cj*0))
  if(p >= 1){
	warning("Group boundaries contain Inf")
      }
}

lines.ogive <- function (data, xlim, ylim, col)
{
  cj <- data$cj
  nj <- data$nj
  
  Fnt <- approxfun(cj, cumsum(nj) / sum(nj),
                 yleft = 0, yright = 1, method = "linear")

  lines(cj, Fnt(cj), xlim = xlim, ylim = ylim, col = col, type = "o")
}

  
