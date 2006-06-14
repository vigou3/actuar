mean.AggregateDist <- function(x) sum(x$fs * (1:length(x$fs)))

var.AggregateDist <- function(x) sum(x$fs * (1:length(x$fs))^2) - (mean.AggregateDist(x))^2
