### ===== actuar: an R package for Actuarial Science =====
###
### Ogive pour données groupées

ogive <- function(data, ...)
{
  class(data) <- "ogive"
  plot(data)
}


plot.ogive <- function (data, ...)
{
  cj <- data$cj
  nj <- data$nj

  Fnt <- approxfun(cj, cumsum(nj)/sum(nj),
                 yleft=0, yright=1, method="linear")

  plot(cj, Fnt(cj), type="o", pch=20)
}
  
