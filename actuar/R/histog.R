### ===== actuar: an R package for Actuarial Science =====
###
### Histogramme pour données groupées

histog <- function(x, xlim = NULL, ylim = NULL, xlab = "group boundaries", ylab = "f(x)")
{
  class(x) <- "histog"
  plot(x, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
}


plot.histog <- function (x, y = NULL, xlim = NULL, ylim = NULL, xlab = "group boundaries", ylab = "f(x)", ...)
{
  cj <- x$cj
  nj <- x$nj

  fnt <- approxfun(cj, c(0, nj[-1] / (sum(nj) * diff(cj))),
                 yleft = 0, yright = 0, f = 1, method = "constant")

  plot(cj, fnt(cj), main = "Histogram", xlim = xlim, ylim = ylim, xlab = xlab, ylab =ylab, type = "S", frame = FALSE) 
  segments(cj, 0, cj, fnt(cj))             
  segments(0, 0, max(cj), 0)
  
  p<-sum(is.na(cj*0))
  if(p >= 1){
	warning("Inf in group boundaries")
     } 
}
  
