### ===== actuar: an R package for Actuarial Science =====
###
### 'summary' method for 'aggregateDist' objects 
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>, Louis-Philippe Pouliot

summary.aggregateDist <- function(object, ...)
    {structure(object, class = c("summary.aggregateDist", class(object)))}

print.summary.aggregateDist <- function(x, ...)
{
    cat(ifelse(comment(x) %in% c("Normal approximation","Normal Power approximation"),
               "Aggregate Claim Amounts CDF:\n",
               "Aggregate Claim Amounts Empirical CDF:\n"))
    q <- quantile(x, p = c(0.25, 0.5, 0.75))
    expectation <- mean(x)
    
    if (comment(x) %in% c("Normal approximation","Normal Power approximation"))
    {
        min <- 0; max <- NA
    }
    else
    {
        max <- tail(eval(expression(x), environment(x)), 1)
        min <- head(eval(expression(x), environment(x)), 1)
    }
    res <- c(min, q[1:2], expectation, q[3], max)
    names(res) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
    print(res)
}
