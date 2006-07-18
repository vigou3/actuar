form1 <- ~ x1 + x1:x2 | a1/a2/a3
form <- form1[[length(form1)]] ## class 'call', devient de longueur '3'

asOneSidedFormula(form[[3]]) ## côté droit

asOneSidedFormula(form[[2]]) ## côté gauche

sep <- "/"
val <- splitFormula(asOneSidedFormula(form[[2]]), sep = sep)
names(val) <- unlist(lapply(val, function(el) deparse(el[[2]])))
as.formula(eval(parse(text = paste("~", paste(names(val), collapse = "/")))))

getGroupsFormula(form1)
getGroupsFormula(form1, asList = TRUE)
getCovariateFormula(form1)

library(actuar)
data(hachemeister) ; hach <- hachemeister$claims

dframepf <- function(x, contracts, subs = NULL)
{
    d <- dim(x)
    nrows <- d[1]
    ir <- seq(length = nrows)
    ncols <- d[2] + 1 + length(subs)
    ic <- seq(length = d[2])    
    collabs <- c(names(subs), "contract", paste("Y", ic, sep = ""))
    value <- vector("list", ncols)
    for (i in (ic)) value[[i+1+length(subs)]] <- as.vector(x[, i])
    if (!is.null(subs))
        value[1:length(subs)] <- subs
    value[[length(subs)+1]] <- as.factor(contracts)
    row.names <- as.character(ir)
    names(value) <- collabs
    attr(value, "row.names") <- row.names
    class(value) <- "data.frame"
    value
}
x <- dframepf(hach, contract = c("4B12", "66FZ", "76H5", "F44R", "FGG5"), list(unit = c(1,1,2,3,2)))
cm(formula = ~ contract | Y1 + Y2 + Y3, data = x)


cm(formula = ~ contract | ., data = x)
cm(formula = ~ contract | Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11+Y12, data = x)
cm(formula = ~ contract | ., data = x, subset = (x$unit == 1))
