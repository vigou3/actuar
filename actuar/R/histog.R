### ===== actuar: an R package for Actuarial Science =====
###
### Histogramme pour données groupées

histog <- function(data, ...)
{
  class(data) <- "histog"
  plot(data)
}


plot.histog <- function (data, ...)
{
  cj <- data$cj
  nj <- data$nj

  fnt <- approxfun(cj, c(0, nj[-1]/(sum(nj) * diff(cj))),
                 yleft=0, yright=0, f=1, method="constant")

  plot(cj, fnt(cj), type="S", frame=FALSE) 
  segments(cj, 0, cj, fnt(cj))             
  segments(0, 0, max(cj), 0)    
}
  
